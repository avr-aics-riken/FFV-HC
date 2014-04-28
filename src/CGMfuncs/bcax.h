#ifndef BCAX_H
#define BCAX_H

#include "real.h"

extern "C" {
	void bc_x3_poiseuille_u_(
					real *x, real *xc
				, int *sz, int *g
				, double *center
				, double *org, double *blockSize, double *cellSize);

	void bc_aw_poiseuille_u_(
					real* Ap,
					real* Aw,
					real* b,
					real* xc,
					int* sz, int* g,
					double *center,
					double *org, double *blockSize, double *cellSize);

	void bc_x3_poiseuille_p_(
					real *x, real *xc
				, int *sz, int *g
				, double *center
				, double *org, double *blockSize, double *cellSize);

	void bc_aw_poiseuille_p_(
					real* Ap,
					real* Aw,
					real* b,
					real* xc,
					int* sz, int* g,
					double *center,
					double *org, double *blockSize, double *cellSize);
}

#endif

