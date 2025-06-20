#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include "vtk.h"
#include "data.h"

char checkpoint_basename[1024];
char result_filename[1024];

/**
 * @brief Set the default basename for file output to out/lid-cavity
 * 
 */
void set_default_base() {
    set_basename("out/lid-cavity");
}

/**
 * @brief Set the basename for file output
 * 
 * @param base Basename string
 */
void set_basename(char *base) {
    checkpoint_basename[0] = '\0';
    result_filename[0] = '\0';
    sprintf(checkpoint_basename, "%s-%%d.vtr", base);
    sprintf(result_filename, "%s.vtr", base);
}

/**
 * @brief Get the basename for file output
 * 
 * @return char* Basename string
 */
char *get_basename() {
    return checkpoint_basename;
}

/**
 * @brief Write a checkpoint VTK file (with the iteration number in the filename)
 * 
 * @param iteration The current iteration number
 * @return int Return whether the write was successful
 */
int write_checkpoint(int iteration) { 
    char filename[1024];
    sprintf(filename, checkpoint_basename, iteration);
    return write_vtk(filename, iteration, 0.0);
}

/**
 * @brief Write the final output to a VTK file
 * 
 * @return int Return whether the write was successful
 */
int write_result() {
    return write_vtk(result_filename, 0, 0.0);
}

/**
 * @brief Write a VTK file with the current state of the simulation
 * 
 * @param filename The filename to write out
 * @return int Return whether the write was successful
 */
int write_vtk(char* filename, int iters, double t) {
    FILE * f = fopen(filename, "w");
    if (f == NULL) {
        perror("Error");
        return -1;
    }

	fprintf(f, "<?xml version=\"1.0\"?>\n");
	fprintf(f, "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
	fprintf(f, "<RectilinearGrid WholeExtent=\"0 %d 0 %d 0 0\">\n", (nx-1), (ny-1));
    fprintf(f, "<FieldData>\n");
    fprintf(f, "<DataArray type=\"Float64\" Name=\"TIME\" NumberOfTuples=\"1\" format=\"ascii\">\n");
    fprintf(f, "%.12e\n", t);
    fprintf(f, "</DataArray>\n");
	fprintf(f, "<DataArray type=\"Int32\" Name=\"CYCLE\" NumberOfTuples=\"1\" format=\"ascii\">\n");
    fprintf(f, "%d\n", iters);
    fprintf(f, "</DataArray>\n");
    fprintf(f, "</FieldData>\n");
    fprintf(f, "<Piece Extent=\"0 %d 0 %d 0 0\">\n", (nx-1), (ny-1));

    fprintf(f, "<Coordinates>\n");
    fprintf(f, "<DataArray type=\"Float64\" Name=\"X\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"%lf\">\n", X);
    for (int i = 0; i < nx; i++) fprintf(f, "%lf ", x[i]);
    fprintf(f, "\n</DataArray>\n");
    fprintf(f, "<DataArray type=\"Float64\" Name=\"Y\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"%lf\">\n", Y);
    for (int i = 0; i < ny; i++) fprintf(f, "%lf ", y[i]);
    fprintf(f, "\n</DataArray>\n");
    fprintf(f, "<DataArray type=\"Float64\" Name=\"Z\" format=\"ascii\">\n");
    fprintf(f, "0.0\n");
    fprintf(f, "</DataArray>\n");
    fprintf(f, "</Coordinates>\n");

    fprintf(f, "<PointData Vectors=\"uv\">\n");
    fprintf(f, "<DataArray type=\"Float64\" Name=\"uv\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            fprintf(f, "%.12e %.12e 0\n", u[i][j], v[i][j]);
    fprintf(f, "</DataArray>\n");
    fprintf(f, "</PointData>\n");
    
    fprintf(f, "<CellData Scalars=\"p\">\n");
    fprintf(f, "<DataArray type=\"Float64\" Name=\"p\" format=\"ascii\">\n");
    for (int i = 0; i < nx-1; i++) {
        for (int j = 0; j < ny-1; j++)
            fprintf(f, "%.12e ", p[i][j]);
        fprintf(f, "\n");
    }
    fprintf(f, "</DataArray>\n");
    fprintf(f, "</CellData>\n");
    
    fprintf(f, "</Piece>\n");
    fprintf(f, "</RectilinearGrid>\n");
    fprintf(f, "</VTKFile>\n");
    
	fclose(f);
	return 0;
}