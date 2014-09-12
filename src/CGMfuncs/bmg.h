/// @file  bmg.h
/// @brief Basic subprograms for the multigrid method for Cartesian grid data structure

#ifndef MG_H
#define MG_H

#include "real.h"

extern "C" {
	void mg_restriction_1_(real *bc, real *r, int *sz, int *g);
	void mg_prolongation_1_(real *x, real *xc, int *sz, int *g);
	void mg_restriction_2_(real *bc, real *r, int *sz, int *g);
	void mg_prolongation_2_(real *x, real *xc, int *sz, int *g);

	void mg_restriction_1_mp_(real *bc, real *r, real *phi, real *psi, real *rhol, real *rhog, int *sz, int *g);
	void mg_prolongation_1_mp_(real *x, real *xc, real *phi, real *psi, real *rhol, real *rhog, int *sz, int *g);

	void mg_restriction_1_ls_(real *bc, real *r, int *sz, int *g); 
	void mg_restriction2_(real *Ac, real *bc, real *A, real *r, int *sz, int *g);
}

#endif


