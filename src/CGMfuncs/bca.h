/// @file  bca.h
/// @brief Functions to compute the elements of the coefficient matrix at the outer boundaries of the simulation box for Cartesian grid data structure (for symmetric/asymmetric 7-band matrices)

#ifndef BCA_H
#define BCA_H

#include "real.h"

extern "C" {
	void bc_aw_d_(real* Ap, real* Aw, real* b, real* xc, int* sz, int* g);
	void bc_ae_d_(real* Ap, real* Ae, real* b, real* xc, int* sz, int* g);
	void bc_as_d_(real* Ap, real* As, real* b, real* xc, int* sz, int* g);
	void bc_an_d_(real* Ap, real* An, real* b, real* xc, int* sz, int* g);
	void bc_ab_d_(real* Ap, real* Ab, real* b, real* xc, int* sz, int* g);
	void bc_at_d_(real* Ap, real* At, real* b, real* xc, int* sz, int* g); 

	void bc_aw_n_(real* Ap, real* Aw, real* b, real* xc, int* sz, int* g);
	void bc_ae_n_(real* Ap, real* Ae, real* b, real* xc, int* sz, int* g);
	void bc_as_n_(real* Ap, real* As, real* b, real* xc, int* sz, int* g);
	void bc_an_n_(real* Ap, real* An, real* b, real* xc, int* sz, int* g);
	void bc_ab_n_(real* Ap, real* Ab, real* b, real* xc, int* sz, int* g);
	void bc_at_n_(real* Ap, real* At, real* b, real* xc, int* sz, int* g);
}

#endif

