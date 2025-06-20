#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <time.h>

#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#include "args.h"
#include "vtk.h"
#include "data.h"
#include "setup.h"

// Global timers
double total_time_build_rhs = 0.0;
double total_time_solve_poissons = 0.0;
double total_time_update_velocities = 0.0;
double total_time_apply_boundary = 0.0;

#define IDX(i,j) ((i)*(ny) + (j))

// Timer helper function
double get_time_seconds() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}


/**
 * @brief Build the right hand side of the Poisson equation on GPU.
 */
__global__ void kernel_build_rhs(double *u_d, double *v_d, double *b_d,
                                 double rho, double dt,
                                 double dx, double dy,
                                 int nx, int ny)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    
    if (i >= 1 && i < nx-1 && j >= 1 && j < ny-1) {
        double dudx = (u_d[IDX(i, j+1)] - u_d[IDX(i, j-1)]) / (2.0 * dx);
        double dvdy = (v_d[IDX(i+1, j)] - v_d[IDX(i-1, j)]) / (2.0 * dy);
        double dudx_sq = dudx * dudx;
        double dvdx = (v_d[IDX(i, j+1)] - v_d[IDX(i, j-1)]) / (2.0 * dx);
        double dudy = (u_d[IDX(i+1, j)] - u_d[IDX(i-1, j)]) / (2.0 * dy);
        double dvdy_sq = dvdy * dvdy;
        
        b_d[IDX(i,j)] = rho * (1.0 / dt) * (dudx + dvdy)
                        - dudx_sq - 2.0 * (dudy * dvdx) - dvdy_sq;
    }
}

/**
 * @brief Apply boundary conditions on GPU.
 */
__global__ void kernel_apply_boundary(double *u_d, double *v_d,
                                      int nx, int ny)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    // bottom row (j=0)
    if (i < nx) {
        u_d[IDX(i,0)] = 0.0;
        v_d[IDX(i,0)] = 0.0;
    }
    // top row (j=ny-1)
    if (i < nx && ny > 1) {
        u_d[IDX(i,ny-1)] = 0.0;
        v_d[IDX(i,ny-1)] = 0.0;
    }

    int j = blockIdx.y * blockDim.y + threadIdx.y;
    // left column (i=0)
    if (j < ny) {
        u_d[IDX(0,j)] = 0.0;
        v_d[IDX(0,j)] = 0.0;
    }
    // right column (i=nx-1)
    if (j < ny && nx > 1) {
        u_d[IDX(nx-1,j)] = 1.0; // // set the velocity on the cavity lid to 1.0
        v_d[IDX(nx-1,j)] = 0.0;
    }
}

/**
 * @brief One step of Poisson solver on GPU.
 */
__global__ void kernel_solve_poissons_step(double *p_d, double *pn_d, double *b_d,
                                           double dx, double dy,
                                           int nx, int ny)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < nx && j < ny) {
        // copy old data
        pn_d[IDX(i,j)] = p_d[IDX(i,j)];
    }

    __syncthreads();  

    if (i >= 1 && i < nx-1 && j >= 1 && j < ny-1) {
        p_d[IDX(i,j)] =
            (((pn_d[IDX(i,  j+1)] + pn_d[IDX(i,  j-1)]) * (dy*dy) +
              (pn_d[IDX(i+1,j)] + pn_d[IDX(i-1,j)]) * (dx*dx))
             / (2.0 * (dx*dx + dy*dy))
             - (dx*dx) * (dy*dy) / (2.0 * (dx*dx + dy*dy))
             * b_d[IDX(i,j)]);
    }
}

/**
 * @brief Apply poisson boundary conditions to p on GPU.
 */
__global__ void kernel_poisson_boundary(double *p_d, int nx, int ny)
{
    // for i in [0..nx-1], fix top/bottom
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < nx) {
        p_d[IDX(i, ny-1)] = p_d[IDX(i, ny-2)];
        p_d[IDX(i, 0)]    = p_d[IDX(i, 1)];
    }

    // for j in [0..ny-1], fix left/right
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (j < ny) {
        p_d[IDX(0, j)] = p_d[IDX(1, j)];
        p_d[IDX(nx-1, j)] = 0.0;
    }
}

/**
 * @brief Update velocities on GPU.
 */
__global__ void kernel_update_velocities(double *u_d, double *v_d,
                                         double *un_d, double *vn_d,
                                         double *p_d,
                                         double rho, double nu, double dt,
                                         double dx, double dy,
                                         int nx, int ny)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < nx && j < ny) {
        // store old velocities
        un_d[IDX(i,j)] = u_d[IDX(i,j)];
        vn_d[IDX(i,j)] = v_d[IDX(i,j)];
    }

    __syncthreads();

    if (i >= 1 && i < nx-1 && j >= 1 && j < ny-1) {
        // u-velocity
        double unij   = un_d[IDX(i,j)];
        double uni_jm = un_d[IDX(i,j-1)];
        double uni_mj = un_d[IDX(i-1,j)];
        double p_jp   = p_d[IDX(i,j+1)];
        double p_jm   = p_d[IDX(i,j-1)];

        u_d[IDX(i,j)] = unij
            - unij * dt / dx * (unij - uni_jm)
            - vn_d[IDX(i,j)] * dt / dy * (unij - uni_mj)
            - dt / (2.0 * rho * dx) * (p_jp - p_jm)
            + nu * (dt / (dx * dx) *
                    (un_d[IDX(i,j+1)] - 2.0 * unij + uni_jm)
                  + dt / (dy * dy) *
                    (un_d[IDX(i+1,j)] - 2.0 * unij + uni_mj));

        // v-velocity
        double vnij   = vn_d[IDX(i,j)];
        double vni_jm = vn_d[IDX(i,j-1)];
        double vni_mj = vn_d[IDX(i-1,j)];
        double p_ip   = p_d[IDX(i+1,j)];
        double p_im   = p_d[IDX(i-1,j)];

        v_d[IDX(i,j)] = vnij
            - un_d[IDX(i,j)] * dt / dx * (vnij - vni_jm)
            - vnij          * dt / dy * (vnij - vni_mj)
            - dt / (2.0 * rho * dy) * (p_ip - p_im)
            + nu * (dt / (dx * dx) *
                    (vn_d[IDX(i,j+1)] - 2.0 * vnij + vni_jm)
                  + dt / (dy * dy) *
                    (vn_d[IDX(i+1,j)] - 2.0 * vnij + vni_mj));
    }
}

/**
 * @brief Build the right hand side of the Poisson equation.
 */
void build_rhs() {
    double start_time = get_time_seconds();

    // Copy from host -> device
    cudaMemcpy(d_u,  u[0], nx*ny*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_v,  v[0], nx*ny*sizeof(double), cudaMemcpyHostToDevice);

    // Launch kernel
    dim3 blockSize(16,16);
    dim3 gridSize((nx+blockSize.x-1)/blockSize.x, (ny+blockSize.y-1)/blockSize.y);
    kernel_build_rhs<<<gridSize, blockSize>>>(d_u, d_v, d_b, rho, dt, dx, dy, nx, ny);

    cudaDeviceSynchronize();

    // Copy back to host
    cudaMemcpy(b[0], d_b, nx*ny*sizeof(double), cudaMemcpyDeviceToHost);

    total_time_build_rhs += get_time_seconds() - start_time;
}
/**
 * @brief Apply the boundary conditions to the u and v arrays.
 */
void apply_boundary() {
    double start_time = get_time_seconds();

    // Copy current u, v to device
    cudaMemcpy(d_u, u[0], nx*ny*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_v, v[0], nx*ny*sizeof(double), cudaMemcpyHostToDevice);

    dim3 blockSize(16,16);
    dim3 gridSize((nx+blockSize.x-1)/blockSize.x, (ny+blockSize.y-1)/blockSize.y);
    kernel_apply_boundary<<<gridSize, blockSize>>>(d_u, d_v, nx, ny);
    
    cudaDeviceSynchronize();

    // Copy back to host
    cudaMemcpy(u[0], d_u, nx*ny*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(v[0], d_v, nx*ny*sizeof(double), cudaMemcpyDeviceToHost);

    total_time_apply_boundary += get_time_seconds() - start_time;
}
/**
 * @brief Solve the Poisson equation to calculate the pressure across the domain.
 */
void solve_poissons() {
    double start_time = get_time_seconds();

    // Copy p, b to device
    cudaMemcpy(d_p, p[0], nx*ny*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, b[0], nx*ny*sizeof(double), cudaMemcpyHostToDevice);

    //iteration on GPU
    dim3 blockSize(16,16);
    dim3 gridSize((nx+blockSize.x-1)/blockSize.x, (ny+blockSize.y-1)/blockSize.y);

    int poisson_iters = 0;
    while (poisson_iters < nit) {
        kernel_solve_poissons_step<<<gridSize, blockSize>>>(d_p, d_u, d_b, dx, dy, nx, ny);

        cudaDeviceSynchronize();

        // apply boundary to p
        kernel_poisson_boundary<<<gridSize, blockSize>>>(d_p, nx, ny);

        cudaDeviceSynchronize();

        poisson_iters++;
    }

    // Copy back to host
    cudaMemcpy(p[0], d_p, nx*ny*sizeof(double), cudaMemcpyDeviceToHost);

    total_time_solve_poissons += get_time_seconds() - start_time;
}
/**
 * @brief Update the u and v velocity arrays using the calculated pressure.
 */
void update_velocities() {
    double start_time = get_time_seconds();

    // Copy relevant arrays to device
    cudaMemcpy(d_u,  u[0], nx*ny*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_v,  v[0], nx*ny*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_p,  p[0], nx*ny*sizeof(double), cudaMemcpyHostToDevice);

    dim3 blockSize(16,16);
    dim3 gridSize((nx+blockSize.x-1)/blockSize.x, (ny+blockSize.y-1)/blockSize.y);
    kernel_update_velocities<<<gridSize, blockSize>>>(d_u, d_v, d_un, d_vn, d_p, rho, nu, dt, dx, dy, nx, ny);

    cudaDeviceSynchronize();

    // Copy updated velocities back to host
    cudaMemcpy(u[0], d_u, nx*ny*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(v[0], d_v, nx*ny*sizeof(double), cudaMemcpyDeviceToHost);

    total_time_update_velocities += get_time_seconds() - start_time;
}
/**
 * @brief The main routine that sets up the problem and executes the timestepping routines
 * 
 * @param argc The number of arguments passed to the program
 * @param argv An array of the arguments passed to the program
 * @return int The return value of the application
 */
int main(int argc, char *argv[]) {
    double program_start_time = get_time_seconds();

    set_defaults();
    parse_args(argc, argv);
    setup();

    printf("Running problem size %f x %f on a %d x %d grid.\n", X, Y, nx, ny);

    if (verbose) print_opts();

    allocate_arrays();

    x = (double*)calloc(nx, sizeof(double));
    for (int i = 0; i < nx; i++) x[i] = (2.0 / (nx-1)) * i;
    y = (double*)calloc(ny, sizeof(double));
    for (int i = 0; i < ny; i++) y[i] = (2.0 / (ny-1)) * i;

    int iters = 0;
    double t = 0.0;

    while (iters < n_iters) {
        build_rhs();
        solve_poissons();
        update_velocities();
        apply_boundary();

        if (iters % output_freq == 0) {
            printf("Step %8d, Time: %14.8e (dt: %14.8e)\n", iters, t, dt);

            if ((!no_output) && (enable_checkpoints))
                write_checkpoint(iters);
        }

        iters++;
        t += dt;
    }

    printf("Step %8d, Time: %14.8e (dt: %14.8e)\n", iters, t, dt);
    printf("Simulation complete.\n");

    if (!no_output)
        write_result();

    free_arrays();

    double program_end_time = get_time_seconds();
    printf("\nExecution Summary:\n");
    printf("Total Program Time: %f seconds\n", program_end_time - program_start_time);
    printf("Total Time in build_rhs: %f seconds\n", total_time_build_rhs);
    printf("Total Time in solve_poissons: %f seconds\n", total_time_solve_poissons);
    printf("Total Time in update_velocities: %f seconds\n", total_time_update_velocities);
    printf("Total Time in apply_boundary: %f seconds\n", total_time_apply_boundary);

    exit(0);
}
