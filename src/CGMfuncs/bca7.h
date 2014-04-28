/// @file  bca7.h
/// @brief Functions to compute the coefficient matrix (BCA) at the outer boundaries of the simulation box for Cartesian grid data structure (for symmetric/asymmetric matrices)

#ifndef BCA7_H
#define BCA7_H

#include "real.h"

extern "C" {
	void bc_a1_d_7_(real *A, real *b, real *xc, int *sz, int *g);
	void bc_a3_d_7_(real *A, real *b, real *xc, int *sz, int *g);
	void bc_a2_d_7_(real *A, real *b, real *xc, int *sz, int *g);
	void bc_a4_d_7_(real *A, real *b, real *xc, int *sz, int *g);
	void bc_a5_d_7_(real *A, real *b, real *xc, int *sz, int *g);
	void bc_a6_d_7_(real *A, real *b, real *xc, int *sz, int *g); 

	void bc_a1_n_7_(real *A, real *b, real *xc, int *sz, int *g);
	void bc_a3_n_7_(real *A, real *b, real *xc, int *sz, int *g);
	void bc_a2_n_7_(real *A, real *b, real *xc, int *sz, int *g);
	void bc_a4_n_7_(real *A, real *b, real *xc, int *sz, int *g);
	void bc_a5_n_7_(real *A, real *b, real *xc, int *sz, int *g);
	void bc_a6_n_7_(real *A, real *b, real *xc, int *sz, int *g);
}

#endif

