/// @file  bsf3d.h
/// @brief Basic subprograms for 3D scalar fields for Cartesian grid data structure

#ifndef BSF3D_H
#define BSF3D_H

#include "real.h"

extern "C" {
	void sf3d_copy_x2_(
					real *x,
					real *xc,
					int *sz, int *g);

	void sf3d_calc_stats_(
					real *sum,
					real *max,
					real *min,
					real *absmax,
					real *absmin,
					real *data,
					int *sz, int *g);
}

#endif


