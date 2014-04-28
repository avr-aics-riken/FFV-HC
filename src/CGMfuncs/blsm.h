/// @file  blsm.h
/// @brief Basic subprograms for the level-set method (BLSM) for Cartesian grid data structure

#ifndef BLSM_H
#define BLSM_H

#include "real.h"

extern "C" {
	void lsm_calc_n_(
					real *nx, real *ny, real *nz, real *divn,
					real *phi,
					int *sz, int *g);
	void lsm_calc_vin_(
					real *vin,
					real *phi, real *psi,
					real *v,
					real *nx, real *ny, real *nz,
					int *sz, int *g);
	void lsm_calc_i_(
					real *fi,
					real *f0,
					real *phi,
					real *phinx, real *phiny, real *phinz,
					real *psi,
					int *sz, int *g);
	void lsm_update_phi_v_(
					real *phi1_,
					real *phi0_,
					real *v0_,
					real *psi_,
					int *sz, int *g);
	void lsm_update_phi_v_quick_(
					real *phi1_,
					real *phi0_,
					real *v0_,
					real *psi_,
					int *sz, int *g);
	void lsm_update_phi_v_weno3_(
					real *phi1_,
					real *phi0_,
					real *v0_,
					int *sz, int *g);
	void lsm_update_phi_f_(
					real *phi1_,
					real *phi0_,
					real *phif,
					int *sz, int *g);
	void lsm_update_phi_f_2_(
					real *phi1_,
					real *phi0_,
					real *phif,
					real *psi,
					int *sz, int *g);
	void lsm_calc_signeddistance_(
					real *d,
					real *phi,
					int *sz, int *g);
	void lsm_calc_signeddistance2_(
					real *d,
					real *phi,
					real *phase,
					int *sz, int *g);
	void lsm_reinit_core_0_(
					real *phi1_,
					real *phi0_,
					real *dtau,
					real *phii_,
					real *phid_,
					int *sz, int *g);
	void lsm_reinit_core_(
					real *phi1_,
					real *phi0_,
					real *dtau,
					real *phii_,
					real *phid_,
					real *psi_,
					int *sz, int *g);
	void lsm_reinit_s_core_(
					real *phi1_,
					real *phi0_,
					real *dtau,
					real *psi_,
					real *psinx_, real *psiny_, real *psinz_,
					real *thetas,
					int *sz, int *g);
	void lsm_smooth_core_(
					real *phi1_,
					real *phi0_,
					real *dtau,
					int *sz, int *g);
	void lsm_extend_core_(
					real *data1_,
					real *data0_,
					real *dtau,
					real *phi,
					real *nx, real *ny, real *nz,
					real *psi,
					int *sz, int *g);
	void lsm_extend_s_core_(
					real *data1_,
					real *data0_,
					real *dtau,
					real *psi,
					real *psinx, real *psiny, real *psinz,
					int *sz, int *g);
}

#endif

