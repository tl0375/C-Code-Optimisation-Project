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

	nx = 500;
	ny = 500;
	
	n_iters = 5000;
	nit = 100;

	set_default_base();
}

/**
 * @brief Set up some of the values required for computation after arguments have been loaded
 * 
 */
void setup() {
	dx = X / (nx - 1);
	dy = Y / (ny - 1);

	// set dt to a reasonable value based on the properties of the grid
	dt = 1.0 / (1.0 / (dx * dx) + 1.0 / (dy * dy)) * 0.9 / 2.0;
}

/**
 * @brief Allocate all of the arrays used for computation
 * 
 */
void allocate_arrays() {
	u = alloc_2d_array(nx, ny);
	v = alloc_2d_array(nx, ny);
	p = alloc_2d_array(nx, ny);
	b = alloc_2d_array(nx, ny);

	// temp storage
	un = alloc_2d_array(nx, ny);
	vn = alloc_2d_array(nx, ny);
}

/**
 * @brief Free all of the arrays used for the computation
 * 
 */
void free_arrays() {
	free_2d_array(u);
	free_2d_array(v);
	free_2d_array(p);
	free_2d_array(b);
	free_2d_array(un);
	free_2d_array(vn);
}
