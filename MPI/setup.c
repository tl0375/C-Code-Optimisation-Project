#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include "vtk.h"
#include "data.h"
#include "setup.h"

/**
 * @brief Set up some default values before arguments have been loaded
 * 
 */
void set_defaults() {
    X = 2.0;
    Y = 2.0;

    nx = 504;
    ny = 504;

    n_iters = 5000;
    nit = 100;

    set_default_base();
}
/**
 * @brief Set up some of the values required for computation after arguments have been loaded
 * 
 */
void setup() {
    // Simple 1D decomposition in x
    local_nx = nx / size;  
    local_nx_start = rank * local_nx;

    // Add 2 for ghost layers on each rank
    local_nx_intern = local_nx + 2;

    dx = X / (nx - 1);
    dy = Y / (ny - 1);

    // set dt to a reasonable value based on the properties of the grid
    dt = 1.0 / (1.0 / (dx*dx) + 1.0 / (dy*dy)) * 0.9 / 2.0;
}

void allocate_arrays() {
    // Now allocate local_nx_intern in x, and ny in y
    u  = alloc_2d_array(local_nx_intern, ny);
    v  = alloc_2d_array(local_nx_intern, ny);
    p  = alloc_2d_array(local_nx_intern, ny);
    b  = alloc_2d_array(local_nx_intern, ny);

    un = alloc_2d_array(local_nx_intern, ny);
    vn = alloc_2d_array(local_nx_intern, ny);
}

void free_arrays() {
    free_2d_array(u);
    free_2d_array(v);
    free_2d_array(p);
    free_2d_array(b);
    free_2d_array(un);
    free_2d_array(vn);
}
