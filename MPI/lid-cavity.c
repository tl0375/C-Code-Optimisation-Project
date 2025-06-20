#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <mpi.h>

#include "args.h"
#include "data.h"
#include "setup.h"
#include "vtk.h"

// Global timers
double total_time_build_rhs = 0.0;
double total_time_solve_poissons  = 0.0;
double total_time_update_velocities = 0.0;
double total_time_apply_boundary = 0.0;

//Timer helper function
static double get_time_seconds() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/**
 * @brief Exchange ghost columns in the x direction
 * 
 * @param field 2D array representing computational grid including ghost cells
 */
static void exchange_ghost_cells(double **field) {
    if (size == 1) return;  // No exchange if only 1 rank

    MPI_Status status;
    int left_rank  = (rank == 0)        ? MPI_PROC_NULL : rank - 1;
    int right_rank = (rank == size - 1) ? MPI_PROC_NULL : rank + 1;

    // Send right, receive from left
    if(right_rank != MPI_PROC_NULL) {
        MPI_Send(field[local_nx], ny, MPI_DOUBLE, right_rank, 0, MPI_COMM_WORLD);
    }
    if(left_rank != MPI_PROC_NULL) {
        MPI_Recv(field[0], ny, MPI_DOUBLE, left_rank, 0, MPI_COMM_WORLD, &status);
    }
    // Send left, receive from right
    if(left_rank != MPI_PROC_NULL) {
        MPI_Send(field[1], ny, MPI_DOUBLE, left_rank, 1, MPI_COMM_WORLD);
    }
    if(right_rank != MPI_PROC_NULL) {
        MPI_Recv(field[local_nx+1], ny, MPI_DOUBLE, right_rank, 1, MPI_COMM_WORLD, &status);
    }
}

/**
 * @brief Apply the boundary conditions to the u and v arrays.
 */
void apply_boundary() {
    double start_time = get_time_seconds();

    // top & bottom walls: j=0, j=ny-1 = zero velocity
    for(int i = 0; i < local_nx_intern; i++) {
        u[i][0] = 0;
        u[i][ny-1] = 0;
        v[i][0] = 0;
        v[i][ny-1] = 0;
    }

    // left boundary (rank=0)
    if(rank == 0) {
        for(int j = 0; j < ny; j++) {
            u[1][j] = 0;  
            v[1][j] = 0;
        }
    }

    // right boundary => lid velocity
    if(rank == size - 1) {
        for(int j = 0; j < ny; j++) {
            u[local_nx][j] = 1;
            v[local_nx][j] = 0;
        }
    }

    total_time_apply_boundary += get_time_seconds() - start_time;
}

/**
 * @brief Build the right hand side of the Poisson equation.
 */
void build_rhs() {
    double start_time = get_time_seconds();

    // build the "b" matrix (RHS)
    for(int i = 1; i <= local_nx; i++) {
        for (int j = 1; j < ny-1; j++) {
            b[i][j] = (rho * (1 / dt) *
                ((u[i][j+1] - u[i][j-1]) /
                (2 * dx) + (v[i+1][j] - v[i-1][j]) / (2 * dy)) -
                pow((u[i][j+1] - u[i][j-1]) / (2 * dx), 2) -
                2 * ((u[i+1][j] - u[i-1][j]) / (2 * dy) *
                (v[i][j+1] - v[i][j-1]) / (2 * dx)) -
                pow((v[i+1][j] - v[i-1][j]) / (2 * dy), 2));
        }
    }

    total_time_build_rhs += get_time_seconds() - start_time;
}

/**
 * @brief Solve the Poisson equation to calculate the pressure across the domain.
 */
static void solve_poissons() {
    double start_time = get_time_seconds();

    double **pn = alloc_2d_array(local_nx_intern, ny);

    for(int iter = 0; iter < nit; iter++) {
        exchange_ghost_cells(p);

        for(int i = 0; i < local_nx_intern; i++) {
            for(int j = 0; j < ny; j++) {
                pn[i][j] = p[i][j];
            }
        }

        for(int i = 1; i <= local_nx; i++) {
            for(int j = 1; j < ny - 1; j++) {
                p[i][j] = (((pn[i][j+1] + pn[i][j-1]) * (dy*dy) +
                   (pn[i+1][j] + pn[i-1][j]) * (dx*dx))/ 
                   (2*(dx*dx + dy*dy))- 
                   (dx*dx * dy*dy)/(2*(dx*dx + dy*dy)) * 
                   b[i][j]);
            }
        }

        // boundary conditions for p
        for(int i = 0; i < local_nx_intern; i++) {
            p[i][ny-1] = p[i][ny-2]; // dp/dx = 0 at x = 2
            p[i][0] = p[i][1];    // dp/dx = 0 at x = 0
        }
        if(rank == 0) {
            for(int j = 0; j < ny; j++) {
                p[1][j] = p[2][j];
            }
        }
        if(rank == size - 1) {
            for(int j = 0; j < ny; j++) {
                p[local_nx][j] = 0.0;
            }
        }
    }

    free_2d_array(pn);
    total_time_solve_poissons += get_time_seconds() - start_time;
}

/**
 * @brief Update the u and v velocity arrays using the calculated pressure.
 */
static void update_velocities() {
    double start_time = get_time_seconds();


    for(int i = 0; i < local_nx_intern; i++) {
        for(int j = 0; j < ny; j++) {
            un[i][j] = u[i][j];
            vn[i][j] = v[i][j];
        }
    }

    for(int i = 1; i <= local_nx; i++) {
        for (int j = 1; j < ny-1; j++) {
            u[i][j] = (un[i][j]-
                un[i][j] * dt / dx *
                (un[i][j] - un[i][j-1]) -
                vn[i][j] * dt / dy *
                (un[i][j] - un[i-1][j]) -
                dt / (2 * rho * dx) * (p[i][j+1] - p[i][j-1]) +
                nu * (dt / (dx*dx) *
                (un[i][j+1] - 2 * un[i][j] + un[i][j-1]) +
                dt / (dy*dy) *
                (un[i+1][j] - 2 * un[i][j] + un[i-1][j])));

            v[i][j] = (vn[i][j] -
                un[i][j] * dt / dx *
                (vn[i][j] - vn[i][j-1]) -
                vn[i][j] * dt / dy *
                (vn[i][j] - vn[i-1][j]) -
                dt / (2 * rho * dy) * (p[i+1][j] - p[i-1][j]) +
                nu * (dt / (dx*dx) *
                (vn[i][j+1] - 2 * vn[i][j] + vn[i][j-1]) +
                dt / (dy*dy) *
                (vn[i+1][j] - 2 * vn[i][j] + vn[i-1][j])));
        }
    }

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
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double program_start_time = get_time_seconds();

    set_defaults();
    parse_args(argc, argv);
    setup();  

    if(rank == 0) {
        printf("Running problem size %f x %f on a %d x %d grid, with %d ranks.\n", X, Y, nx, ny, size);
    }
    if(verbose) {
        print_opts();
    }

    allocate_arrays();

    
    x = (double*)calloc(local_nx_intern, sizeof(double));
    for(int i = 1; i <= local_nx; i++) {
        int global_i = local_nx_start + (i - 1);
        x[i] = (2.0 / (nx - 1)) * global_i;
    }
    y = (double*)calloc(ny, sizeof(double));
    for(int j = 0; j < ny; j++) {
        y[j] = (2.0 / (ny - 1)) * j;
    }

    int iters = 0;
    double t  = 0.0;

    while(iters < n_iters) {
        build_rhs();
        solve_poissons();
        update_velocities();
        exchange_ghost_cells(u);
        exchange_ghost_cells(v);
        apply_boundary();

        if(iters % output_freq == 0) {
            if(rank == 0) {
                printf("Step %8d, Time: %.8e (dt: %.8e)\n", iters, t, dt);
            }
            if((!no_output) && enable_checkpoints) {
                write_checkpoint(iters);
            }
        }
        iters++;
        t += dt;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if(rank == 0) {
        printf("Step %8d, Time: %.8e (dt: %.8e)\n", iters, t, dt);
        printf("Simulation complete.\n");
    }

    if(!no_output) {
        write_result();  
    }

    free_arrays();
    free(x);
    free(y);

    double program_end_time = get_time_seconds();
    if(rank == 0) {
        printf("\nExecution Summary:\n");
        printf("Total Program Time: %f seconds\n", program_end_time - program_start_time);
        printf("Total Time in build_rhs: %f seconds\n", total_time_build_rhs);
        printf("Total Time in solve_poissons: %f seconds\n", total_time_solve_poissons);
        printf("Total Time in update_velocities: %f seconds\n", total_time_update_velocities);
        printf("Total Time in apply_boundary: %f seconds\n", total_time_apply_boundary);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
    
    exit(0);
}
