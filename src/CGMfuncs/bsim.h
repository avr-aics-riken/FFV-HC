/// @file  bsim.h
/// @brief Basic subprograms for a sharp interface method (BSIM) for Cartesian grid data structure

#ifndef SIM_H
#define SIM_H

#include "real.h"

extern "C" {
	void sim_calc_c_(
					real *fc,
					real *f,
					real *v,
					real *phi, real *psi,
					real *org,
					int *sz, int *g);
	void sim_calc_c_u_(
					real *fc,
					real *f,
					real *v,
					real *phi, real *psi,
					int *dir,
					real *Omega, real *zmin, real *zmax, real *rmin, real *rmax,
					real *org,
					int *sz, int *g);
	void sim_calc_c_u_u1_(
					real *fc,
					real *f,
					real *v,
					real *phi, real *psi,
					int *dir,
					real *Omega, real *zmin, real *zmax, real *rmin, real *rmax,
					real *org,
					int *sz, int *g);
	void sim_calc_c_u_quick_(
					real *fc,
					real *f,
					real *v,
					real *phi, real *psi,
					int *dir,
					real *Omega, real *zmin, real *zmax, real *rmin, real *rmax,
					real *org,
					int *sz, int *g);
	void sim_calc_c_u_test_(
					real *fc,
					real *f,
					real *v,
					real *phi, real *psi,
					int *dir,
					real *Omega, real *zmin, real *zmax, real *rmin, real *rmax,
					real *org,
					int *sz, int *g);

	void sim_calc_abd_u_(
					real *A, real *b,
					real *u1_, real *u0_,
					real *uc0_, real *ucp_,
					real *ud0_, real *udp_,
					real *ux, real *uy, real *uz,
					real *phi, real *psi,
					real *rhol, real *rhog,
					real *mul, real *mug,
					real *nl_viscosity, real *nl_viscosity_d,
					int *dir,
					real *Omega, real *zmin, real *zmax, real *rmin, real *rmax,
					real *org,
					int *sz, int *g);

	void sim_calc_abd_u2_(
					real *A, real *b,
					real *u1_, real *u0_,
					real *uc0_, real *ucp_,
					real *ud0_, real *udp_,
					real *ux, real *uy, real *uz,
					real *phi, real *psi,
					real *work0_, real *work1_, real *work2_,
					real *rhol, real *rhog,
					real *mul, real *mug,
					real *nl_viscosity, real *nl_viscosity_d,
					int *dir,
					real *Omega, real *zmin, real *zmax, real *rmin, real *rmax,
					real *org,
					int *sz, int *g);

	void sim_calc_abd_t_(
					real *A, real *b,
					real *t1_, real *t0_,
					real *tc0_, real *tcp_,
					real *td0_, real *tdp_,
					real *phi, real *psi,
					real *rhol, real *rhog, real *rhos,
					real *cpl, real *cpg, real *cps,
					real *kl, real *kg, real *ks,
					real *org,
					int *sz, int *g);

	void sim_add_g_(
					real *ux, real *uy, real *uz,
					real *phi, real *psi,
					real *gx, real *gy, real *gz,
					real *OmegaCoriolis, real *OmegaCentrifugal,
					real *org,
					int *sz, int *g);

	void sim_calc_ab_p_(
					real *A, real *b,
					real *p1_,
					real *v,
					real *p0_,
					real *ux, real *uy, real *uz,
					real *phi, real *psi,
					real *kappa,
					real *rhol, real *rhog,
					real *sigma,
					real *Omega, real *zmin, real *zmax, real *rmin, real *rmax,
					real *org,
					int *sz, int *g);

	void sim_set_refbox_p_(
					real *A, real *b,
					real *x,
					real* xmin, real* ymin, real* zmin,
					real* xmax, real* ymax, real* zmax,
					real* pr,
					real* dx,
					real* org,
					int *sz, int *g);

	void sim_corr_u_(
					real *ux1_, real *uy1_, real *uz1_,
					real *v1_,
					real *ux0_, real *uy0_, real *uz0_,
					real *v0_,
					real *p0_,
					real *phi, real *psi,
					real *kappa,
					real *rhol, real *rhog,
					real *sigma,
					real *gx, real *gy, real *gz,
					real *Omega, real *zmin, real *zmax, real *rmin, real *rmax,
					real *org,
					int *sz, int *g);

	void sim_set_ic_u_(
					real *ux, real *uy, real *uz,
					real *phi, real *psi,
					real *Omega, real *zmin, real *zmax, real *rmin, real *rmax,
					real *org,
					int *sz, int *g);
}

#endif

