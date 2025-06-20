#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include "vtk.h"
#include "data.h"

double X;
double Y;

int nx;
int ny;
int n_iters;
int nit;

double dx;
double dy;

double rho = 1.0;
double nu = 0.1;
double dt = 0.001;

double *x;
double *y;

double **u;
double **v;
double **p;
double **b;

double **un;
double **vn;

/**
 * @brief Allocate a 2D array that is addressable using square brackets
 * 
 * @param m The first dimension of the array
 * @param n The second dimension of the array
 * @return double** A 2D array
 */
double **alloc_2d_array(int m, int n) {
  	double **x;
  	int i;

  	x = (double **)malloc(m*sizeof(double *));
  	x[0] = (double *)calloc(m*n,sizeof(double));
  	for ( i = 1; i < m; i++ )
    	x[i] = &x[0][i*n];
	return x;
}

/**
 * @brief Free a 2D array
 * 
 * @param array The 2D array to free
 */
void free_2d_array(double ** array) {
	free(array[0]);
	free(array);
}
