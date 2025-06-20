#ifndef DATA_H
#define DATA_H

#ifdef __cplusplus
extern "C" {
#endif

// size of the domain
extern double X;
extern double Y;

// dimension of the discretisation
extern int nx;
extern int ny;

// number of iterations
extern int n_iters;

// poisson iterations
extern int nit;

extern double dx;
extern double dy;

// experimental parameters
extern double rho;
extern double nu;
extern double dt;

// arrays for storing the simulation state (on host)
extern double *x;
extern double *y;

extern double **u;
extern double **v;
extern double **p;
extern double **b;

extern double **un;
extern double **vn;

// DEVICE pointers (on GPU)
extern double *d_u;
extern double *d_v;
extern double *d_p;
extern double *d_b;
extern double *d_un;
extern double *d_vn;

// helper functions for 2D arrays
double **alloc_2d_array(int m, int n);
void free_2d_array(double ** array);

#ifdef __cplusplus
}
#endif

#endif
