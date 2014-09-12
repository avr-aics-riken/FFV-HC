/// @file  bsdf.h
/// @brief Basic subprograms for initializing the signed distance functions for Cartesian grid data structure

#ifndef BSDF_H
#define BSDF_H

#include "real.h"

extern "C" {
	void bsdf_liquid_(
					real *phi, real *psi,
					real *xc, real *yc, real *zc,
					real *rc, real *lc, real *wc, real *hc,
					int *sz, int *g,
					int *b, int *gsz);
	void bsdf_pipe_(
					real *phi, real *psi,
					real *xc, real *yc, real *zc,
					real *rc, real *lc, real *wc, real *hc,
					int *sz, int *g,
					int *b, int *gsz);
	void bsdf_pipe2_(
					real *phi, real *psi,
					real *xc, real *yc, real *zc,
					real *rc, real *lc, real *wc, real *hc,
					int *sz, int *g,
					int *b, int *gsz);
	void bsdf_beaker_(
					real *phi, real *psi,
					real *xc, real *yc, real *zc,
					real *rc, real *lc, real *wc, real *hc,
					int *sz, int *g,
					int *b, int *gsz);
	void bsdf_ducky_(
					real *phi, real *psi,
					real *xc, real *yc, real *zc,
					real *rc, real *lc, real *wc, real *hc,
					int *sz, int *g,
					int *b, int *gsz);
	void bsdf_bubble2d_(
					real *phi, real *psi,
					real *xc, real *yc, real *zc,
					real *rc, real *lc, real *wc, real *hc,
					int *sz, int *g,
					int *b, int *gsz);
	void bsdf_dambreak2d_(
					real *phi, real *psi,
					real *xc, real *yc, real *zc,
					real *rc, real *lc, real *wc, real *hc,
					int *sz, int *g,
					int *b, int *gsz);
	void bsdf_bubble_(
					real *phi, real *psi,
					real *xc, real *yc, real *zc,
					real *rc, real *lc, real *wc, real *hc,
					int *sz, int *g,
					int *b, int *gsz);
	void bsdf_drop_(
					real *phi, real *psi,
					real *xc, real *yc, real *zc,
					real *rc, real *lc, real *wc, real *hc,
					int *sz, int *g,
					int *b, int *gsz);
}

#endif

