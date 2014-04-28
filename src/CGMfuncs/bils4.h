/// @file  bils4.h
/// @brief Basic subprograms for iterative linear solvers (BILS) for Cartesian grid data structure (for symmetric matrices) 
#ifndef BILS4_H
#define BILS4_H

#include "real.h"

extern "C" {
	void rbgs_smoother_4_d_(
					real *x,
					real *A, real *b,
					real *param,
					int *color, int *offset,
					int *sz, int *g,
					int *mx,
					real *sx1, real *sx3, real *sx2, real *sx4, real *sx5, real *sx6,
					real *rx1, real *rx3, real *rx2, real *rx4, real *rx5, real *rx6,
					int *node);
	void rbgs_smoother_4_n_(
					real *x,
					real *A, real *b,
					real *param,
					int *color, int *offset,
					int *sz, int *g,
					int *mx,
					real *sx1, real *sx3, real *sx2, real *sx4, real *sx5, real *sx6,
					real *rx1, real *rx3, real *rx2, real *rx4, real *rx5, real *rx6,
					int *node);
	void rbgs_smoother_4_p_(
					real *x,
					real *A, real *b,
					real *param,
					int *color, int *offset,
					int *sz, int *g,
					int *mx,
					real *sx1, real *sx3, real *sx2, real *sx4, real *sx5, real *sx6,
					real *rx1, real *rx3, real *rx2, real *rx4, real *rx5, real *rx6,
					int *node);
	void rbgs_smoother_4_b_(
					real *x,
					real *A, real *b,
					real *param,
					int* color, int* offset,
					int *sz, int *g);
	void calc_ax_4_b_(
					real *Ax,
					real *A, real *x,
					int *sz, int *g);
	void calc_r_4_b_(
					real *r,
					real *A, real *b, real *x,
					int *sz, int *g);
	void calc_rr_4_b_(
					real *rr,
					real *A, real *b, real *x,
					int *sz, int *g);
	void rbgs_smoother_4_s_(
					real *x,
					real *A0, real *A1, real *A2, real *A3, real *b,
					real *param,
					int* color, int* offset,
					int *sz, int *g);
	void calc_ax_4_s_(
					real *Ax,
					real *A0, real *A1, real *A2, real *A3, real *x,
					int *sz, int *g);
	void calc_r_4_s_(
					real *r,
					real *A0, real *A1, real *A2, real *A3, real *b, real *x,
					int *sz, int *g);
	void calc_rr_4_s_(
					real *rr,
					real *A0, real *A1, real *A2, real *A3, real *b, real *x,
					int *sz, int *g);
}

#endif

