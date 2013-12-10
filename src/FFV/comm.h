#ifndef COMM_H
#define COMM_H

#include "real.h"

extern "C" {
	void allreduce_(
					real *total, real *a);
	void allreduce_max_(
					real *total, real *a);
	void allreduce_min_(
					real *total, real *a);

	void dot_all_(
					real *xy, real *y, real *x
				, int *sz, int *g);

	void comm_band_cells_(
					real *x
				, int *sz, int *g
				, int *mx
				, real *sx1, real *sx3, real *sx2, real *sx4, real *sx5, real *sx6
				, real *rx1, real *rx3, real *rx2, real *rx4, real *rx5, real *rx6
				, int *node);
	void comm_band_cells_nb_(
					real *x
				, int *sz, int *g
				, int *mx
				, real *sx1, real *sx3, real *sx2, real *sx4, real *sx5, real *sx6
				, real *rx1, real *rx3, real *rx2, real *rx4, real *rx5, real *rx6
				, int *node);
}

#endif

