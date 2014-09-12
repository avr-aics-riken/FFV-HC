/// @file  bcx.h
/// @brief Functions to compute the right-hand side values at the outer boundaries of the simulation box for Cartesian grid data structure (for symmetric/asymmetric 7-band matrices)

#ifndef BCX_H
#define BCX_H

#include "real.h"

extern "C" {
	void bc_x1_d_(real *x, real *xc, int *sz, int *g);
	void bc_x3_d_(real *x, real *xc, int *sz, int *g);
	void bc_x2_d_(real *x, real *xc, int *sz, int *g);
	void bc_x4_d_(real *x, real *xc, int *sz, int *g);
	void bc_x5_d_(real *x, real *xc, int *sz, int *g);
	void bc_x6_d_(real *x, real *xc, int *sz, int *g);

	void bc_x1_d_f_(real *x, real *xc, int *sz, int *g);
	void bc_x3_d_f_(real *x, real *xc, int *sz, int *g);
	void bc_x2_d_f_(real *x, real *xc, int *sz, int *g);
	void bc_x4_d_f_(real *x, real *xc, int *sz, int *g);
	void bc_x5_d_f_(real *x, real *xc, int *sz, int *g);
	void bc_x6_d_f_(real *x, real *xc, int *sz, int *g);

	void bc_x1_n_(real *x, real *xc, int *sz, int *g);
	void bc_x3_n_(real *x, real *xc, int *sz, int *g);
	void bc_x2_n_(real *x, real *xc, int *sz, int *g);
	void bc_x4_n_(real *x, real *xc, int *sz, int *g);
	void bc_x5_n_(real *x, real *xc, int *sz, int *g);
	void bc_x6_n_(real *x, real *xc, int *sz, int *g);

	void bc_x1_p_(real *x, real *xc, int *sz, int *g);
	void bc_x2_p_(real *x, real *xc, int *sz, int *g);
	void bc_x5_p_(real *x, real *xc, int *sz, int *g);
}

#endif


