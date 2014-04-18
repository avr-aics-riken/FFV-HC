#ifndef BILS_H
#define BILS_H

#include "real.h"

extern "C" {
	void jacobi_smoother_(
					real* x1, real* x0,
					real* Ap, real* Aw, real* Ae, real* As, real* An, real* Ab, real* At,
					real* b,
					real* omega,
					int* sz, int* g);
	void rbgs_smoother_(
					real* x,
					real* Ap, real* Aw, real* Ae, real* As, real* An, real* Ab, real* At,
					real* b,
					real* omega,
					int* color, int* offset,
					int* sz, int* g);
	void jacobi_smoother2_(
					real* x1, real* x0,
					real* Ap, real* Aw, real* Ae, real* As, real* An, real* Ab, real* At,
					real* b,
					real* omega,
					int* sz, int* g);
	void calc_ax_( 
					real* Ax,
					real* Ap, real* Aw, real* Ae, real* As, real* An, real* Ab, real* At,
					real* x,
					int* sz, int* g);
	void calc_r_( 
					real* r,
					real* Ap, real* Aw, real* Ae, real* As, real* An, real* Ab, real* At,
					real* x,
					real* b,
					int* sz, int* g);
	void calc_r2_( 
					real* r,
					real* Ap, real* Aw, real* Ae, real* As, real* An, real* Ab, real* At,
					real* x,
					real* b,
					int* sz, int* g);

	void setup_mask_(
					real* m,
					int* mask,
					int* sz, int* g);
	void copy_mask_(
					real *y, real *x,
					int* mask,
					int *sz, int *g);
	void dot_mask_(
					real *xy, real *y, real *x,
					int* mask,
					int *sz, int *g);
	void calc_ax_mask_( 
					real* Ax,
					real* Ap, real* Aw, real* Ae, real* As, real* An, real* Ab, real* At,
					real* x,
					int* mask,
					int* sz, int* g);
	void jacobi_smoother_mask_(
					real* x1, real* x0,
					real* Ap, real* Aw, real* Ae, real* As, real* An, real* Ab, real* At,
					real* b,
					int* mask,
					real* omega,
					int* sz, int* g);
}

#endif

