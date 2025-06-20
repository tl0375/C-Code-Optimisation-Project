#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <mpi.h>

#include "vtk.h"
#include "data.h"


char checkpoint_basename[1024];
char result_filename[1024];

/**
 * @brief If user does not set a custom base, we use out/lid-cavity
 */
void set_default_base() {
    set_basename("out/lid-cavity");
}

/**
 * @brief Set the base name for file output
 * 
 * @param base Basename string
 */
void set_basename(char *base) {
    sprintf(checkpoint_basename, "%s-%%d.vtr", base);
    sprintf(result_filename, "%s.vtr", base);
}

/**
 * @brief Return the checkpoint basename pattern
 * 
 * @return char* Basename string
 */
char *get_basename() {
    return checkpoint_basename;
}

/**
 * @brief Write a checkpoint file (with iteration number in the filename).
 *        
 * @param iteration The current iteration number
 * @return int Return whether the write was successfu
 */
int write_checkpoint(int iteration) {
    char filename[1024];
    sprintf(filename, checkpoint_basename, iteration);
    return write_vtk(filename, iteration, 0.0);
}

/**
 * @brief Write the final result to a .vtr file.
 *        
 * @return int Return whether the write was successful
 */
int write_result() {
    return write_vtk(result_filename, 0, 0.0);
}

/**
 * @brief Gather data from all ranks (skipping ghost columns) and write to .vtr file
 *
 * @param filename The filename to write out
 * @return int Return whether the write was successful
 */
int write_vtk(char* filename, int iters, double t) {
    // local_nx is how many interior columns this rank owns
    int local_interior_size = local_nx * ny;

    // Flatten local arrays
    double *u_local = (double*)malloc(local_interior_size * sizeof(double));
    double *v_local = (double*)malloc(local_interior_size * sizeof(double));
    double *p_local = (double*)malloc(local_interior_size * sizeof(double));
    double *x_local = (double*)malloc(local_nx * sizeof(double));

    {
        int idx = 0;
        // Flatten u, v, p from [1..local_nx], skipping i=0 and i=local_nx+1 (ghosts)
        for(int i = 1; i <= local_nx; i++) {
            for(int j = 0; j < ny; j++) {
                u_local[idx] = u[i][j];
                v_local[idx] = v[i][j];
                p_local[idx] = p[i][j];
                idx++;
            }
        }
        // Flatten x from [1..local_nx] similarly
        for(int i = 1; i <= local_nx; i++) {
            x_local[i-1] = x[i];
        }
    }

    double *u_global = NULL;
    double *v_global = NULL;
    double *p_global = NULL;
    double *x_global = NULL;

    if(rank == 0) {
        u_global = (double*)malloc(nx * ny * sizeof(double));
        v_global = (double*)malloc(nx * ny * sizeof(double));
        p_global = (double*)malloc(nx * ny * sizeof(double));
        x_global = (double*)malloc(nx * sizeof(double));
    }

    MPI_Gather(x_local, local_nx, MPI_DOUBLE,
               x_global, local_nx, MPI_DOUBLE,
               0, MPI_COMM_WORLD);

    MPI_Gather(u_local, local_interior_size, MPI_DOUBLE,
               u_global, local_interior_size, MPI_DOUBLE,
               0, MPI_COMM_WORLD);

    MPI_Gather(v_local, local_interior_size, MPI_DOUBLE,
               v_global, local_interior_size, MPI_DOUBLE,
               0, MPI_COMM_WORLD);

    MPI_Gather(p_local, local_interior_size, MPI_DOUBLE,
               p_global, local_interior_size, MPI_DOUBLE,
               0, MPI_COMM_WORLD);

    free(u_local);
    free(v_local);
    free(p_local);
    free(x_local);

    if(rank == 0) {
        FILE *f = fopen(filename, "w");
        if(!f) {
            perror("Error opening VTK file");
            free(u_global);
            free(v_global);
            free(p_global);
            free(x_global);
            return -1;
        }

        fprintf(f, "<?xml version=\"1.0\"?>\n");
        fprintf(f, "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
        fprintf(f, "<RectilinearGrid WholeExtent=\"0 %d 0 %d 0 0\">\n", nx-1, ny-1);
        fprintf(f, "<FieldData>\n");
        fprintf(f, "<DataArray type=\"Float64\" Name=\"TIME\" NumberOfTuples=\"1\" format=\"ascii\">\n");
        fprintf(f, "%.12e\n", t);
        fprintf(f, "</DataArray>\n");
        fprintf(f, "<DataArray type=\"Int32\" Name=\"CYCLE\" NumberOfTuples=\"1\" format=\"ascii\">\n");
        fprintf(f, "%d\n", iters);
        fprintf(f, "</DataArray>\n");
        fprintf(f, "</FieldData>\n");
        fprintf(f, "<Piece Extent=\"0 %d 0 %d 0 0\">\n", nx-1, ny-1);

        fprintf(f, "<Coordinates>\n");
        fprintf(f, "<DataArray type=\"Float64\" Name=\"X\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"%lf\">\n", X);
        for(int i = 0; i < nx; i++) fprintf(f, "%lf ", x_global[i]);
        fprintf(f, "\n</DataArray>\n");
        fprintf(f, "<DataArray type=\"Float64\" Name=\"Y\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"%lf\">\n", Y);
        for(int j = 0; j < ny; j++) fprintf(f, "%lf ", y[j]);
        fprintf(f, "\n</DataArray>\n");
        fprintf(f, "<DataArray type=\"Float64\" Name=\"Z\" format=\"ascii\">\n");
        fprintf(f, "0.0\n");
        fprintf(f, "</DataArray>\n");
        fprintf(f, "</Coordinates>\n");

        fprintf(f, "<PointData Vectors=\"uv\">\n");
        fprintf(f, "<DataArray type=\"Float64\" Name=\"uv\" NumberOfComponents=\"3\" format=\"ascii\">\n");

        int offset = 0;
        for(int r = 0; r < size; r++) {
            int chunk_x = nx / size; 
            int chunk_size = chunk_x * ny;
            for(int k = 0; k < chunk_size; k++) {
                double uu = u_global[offset + k];
                double vv = v_global[offset + k];
                fprintf(f, "%.12e %.12e 0\n", uu, vv);
            }
            offset += chunk_size;
        }
        fprintf(f, "</DataArray>\n");
        fprintf(f, "</PointData>\n");

        fprintf(f, "<CellData Scalars=\"p\">\n");
        fprintf(f, "<DataArray type=\"Float64\" Name=\"p\" format=\"ascii\">\n");

        offset = 0;
        for(int r = 0; r < size; r++) {
            int chunk_x = nx / size;
            for(int i = 0; i < chunk_x; i++) {
                for(int j = 0; j < ny-1; j++) {
                    int idx = offset + i*ny + j;
                    fprintf(f, "%.12e ", p_global[idx]);
                }
                fprintf(f, "\n");
            }
            offset += (chunk_x * ny);
        }
        fprintf(f, "</DataArray>\n");
        fprintf(f, "</CellData>\n");
        
        fprintf(f, "</Piece>\n");
        fprintf(f, "</RectilinearGrid>\n");
        fprintf(f, "</VTKFile>\n");

        fclose(f);
        free(u_global);
        free(v_global);
        free(p_global);
        free(x_global);
    }

    return 0;  
}
