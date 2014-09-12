/// @file  blas.h
/// @brief Basic subprograms for linear algebra for Cartesian grid data structure (for symmetric/asymmetric 7-band matrices)

#ifndef BLAS_H
#define BLAS_H

#include "real.h"

extern "C" {
	void fill_(real *x, real *a, int *sz, int *g);
	void fill_2_(real *x, real *a, int *sz, int *g);
	void fill_vf3d_(real *x, real *a, int *sz, int *g, int *ne);

	void copy_(real *y, real *x, int *sz, int *g);
	void copy_integer_(int *y, int *x, int *sz, int *g);
	void add_(real *C, real *A, real *B, int *sz, int *g);
	void ave_(real *C, real *A, real *B, int *sz, int *g);
	void triad_(real *C, real *A, real *B, real *d, int *sz, int *g);
	void scal_(real *y, real *a, int *sz, int *g);
	void axpy_(real *y, real *x, real *a, int *sz, int *g);
	void xpay_(real *y, real *x, real *a, int *sz, int *g);
	void axpyz_(real *z, real *x, real *y, real *a, int *sz, int *g);
	void axpbypz_(real *z, real *x, real *y, real *a, real *b, int *sz, int *g);
	void dot_(real *xy, real *y, real *x, int *sz, int *g);
	void dotx2_(real *xy, real *xz, real *x, real *y, real *z, int *sz, int *g);
}

#endif

