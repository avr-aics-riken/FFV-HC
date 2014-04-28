/// @file  bca4.h
/// @brief Functions to compute the coefficient matrix (BCA) at the outer boundaries of the simulation box for Cartesian grid data structure (for symmetric matrices)

#ifndef BCA4_H
#define BCA4_H

#include "real.h"

extern "C" {
	void bc_a1_d_4_(real *A, real *b, real *xc, int *sz, int *g);
	void bc_a3_d_4_(real *A, real *b, real *xc, int *sz, int *g);
	void bc_a2_d_4_(real *A, real *b, real *xc, int *sz, int *g);
	void bc_a4_d_4_(real *A, real *b, real *xc, int *sz, int *g);
	void bc_a5_d_4_(real *A, real *b, real *xc, int *sz, int *g);
	void bc_a6_d_4_(real *A, real *b, real *xc, int *sz, int *g);

	void bc_a1_n_4_(real *A, real *b, real *xc, int *sz, int *g);
	void bc_a3_n_4_(real *A, real *b, real *xc, int *sz, int *g);
	void bc_a2_n_4_(real *A, real *b, real *xc, int *sz, int *g);
	void bc_a4_n_4_(real *A, real *b, real *xc, int *sz, int *g);
	void bc_a5_n_4_(real *A, real *b, real *xc, int *sz, int *g);
	void bc_a6_n_4_(real *A, real *b, real *xc, int *sz, int *g);
}

#endif //

