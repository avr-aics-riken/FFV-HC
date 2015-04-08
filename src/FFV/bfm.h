#ifndef BFM_H
#define BFM_H

#include "real.h"

extern "C" {
	void bfm_fan_(
			real* fx, real* fy, real* fz,
			real* ux0, real* uy0, real* uz0,
			int* rid,
			int* rid_target,
			real* b,
			real* nx, real* ny, real* nz,
			real* dpmax, real* umax,
			real* dx, real* dt,
			int *sz, int *g);

	void bfm_hex_(
			real* fx, real* fy, real* fz,
			real* ux0, real* uy0, real* uz0,
			int* rid,
			int* rid_target,
			real* b,
			real* nx, real* ny, real* nz,
			real* dpmax, real* umax,
			real* dx, real* dt,
			int *sz, int *g);
}

#endif

