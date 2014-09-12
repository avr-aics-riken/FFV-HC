/// @file  bils7.h
/// @brief Basic subprograms for iterative linear solvers (BILS) for Cartesian grid data structure (for symmetric/asymmetric 7-band matrices) 

#ifndef BILS7_H
#define BILS7_H

#include "real.h"

extern "C" {
/// @brief Red-black Gauss-Seidel smoother for 7-band matries for Cartesian grid data structure
/// @param [out] x the 3D array to be smoothed
/// @param [in]  A the 3D vector array for coefficient matries
/// @param [in]  b the 3D array for the right-hand side
/// @param [in]  param the relaxation parameter
/// @param [in]  color the color id of which the cells are smoothed
/// @param [in]  offset the offset to control the first color to be updated
/// @param [in]  sz the number of cells in the 3D array (NX, NY, NZ)
/// @param [in]  g the size of the guide cells
	void rbgs_smoother_7_b_(
					real *x,
					real *A, real *b,
					real *param,
					int* color, int* offset,
					int *sz, int *g);

/// @brief Compute matrix-vector multiplication for 7-band matries for Cartesian grid data structure \n
///        Ax = A x
/// @param [out] Ax the 3D array for the result
/// @param [in]  A the 3D vector array for coefficient matries
/// @param [in]  x the 3D array to be multiplied
/// @param [in]  sz the number of cells in the 3D array (NX, NY, NZ)
/// @param [in]  g the size of the guide cells
	void calc_ax_7_b_(
					real *Ax,
					real *A, real *x,
					int *sz, int *g);

/// @brief Compute a residual for 7-band matries for Cartesian grid data structure \n
///        r = b - A x
/// @param [out] r the 3D array for the result
/// @param [in]  A the 3D vector array for coefficient matries
/// @param [in]  b the 3D array for the right-hand side
/// @param [in]  x the 3D array to be multiplied
/// @param [in]  sz the number of cells in the 3D array (NX, NY, NZ)
/// @param [in]  g the size of the guide cells
	void calc_r_7_b_(
					real *r,
					real *A, real *b, real *x,
					int *sz, int *g);

/// @brief Compute the square of a residual for 7-band matries for Cartesian grid data structure \n
///        r = b - A x \n
///        rr += r*r
/// @param [out] rr the square of a residual
/// @param [in]  A the 3D vector array for coefficient matries
/// @param [in]  b the 3D array for the right-hand side
/// @param [in]  x the 3D array to be multiplied
/// @param [in]  sz the number of cells in the 3D array (NX, NY, NZ)
/// @param [in]  g the size of the guide cells
	void calc_rr_7_b_(
					real *rr,
					real *A, real *b, real *x,
					int *sz, int *g);

/// @brief Red-black Gauss-Seidel smoother for 7-band matries for Cartesian grid data structure (for coefficinet matrices represented by 7 vectors)
/// @param [out] x the 3D array to be smoothed
/// @param [in]  A the 3D vector array for coefficient matries
/// @param [in]  b the 3D array for the right-hand side
/// @param [in]  param the relaxation parameter
/// @param [in]  color the color id of which the cells are smoothed
/// @param [in]  offset the offset to control the first color to be updated
/// @param [in]  sz the number of cells in the 3D array (NX, NY, NZ)
/// @param [in]  g the size of the guide cells
	void rbgs_smoother_7_s_(
					real *x,
					real *A0, real *A1, real *A2, real *A3, real *A4, real *A5, real *A6, real *b,
					real *param,
					int* color, int* offset,
					int *sz, int *g);

/// @brief Compute matrix-vector multiplication for 7-band matries for Cartesian grid data structure (for coefficinet matrices represented by 7 vectors) \n
///        Ax = A x
/// @param [out] Ax the 3D array for the result
/// @param [in]  A the 3D vector array for coefficient matries
/// @param [in]  x the 3D array to be multiplied
/// @param [in]  sz the number of cells in the 3D array (NX, NY, NZ)
/// @param [in]  g the size of the guide cells
	void calc_ax_7_s_(
					real *Ax,
					real *A0, real *A1, real *A2, real *A3, real *A4, real *A5, real *A6, real *x,
					int *sz, int *g);

/// @brief Compute a residual for 7-band matries for Cartesian grid data structure (for coefficinet matrices represented by 7 vectors) \n
///        r = b - A x
/// @param [out] r the 3D array for the result
/// @param [in]  A the 3D vector array for coefficient matries
/// @param [in]  b the 3D array for the right-hand side
/// @param [in]  x the 3D array to be multiplied
/// @param [in]  sz the number of cells in the 3D array (NX, NY, NZ)
/// @param [in]  g the size of the guide cells
	void calc_r_7_s_(
					real *r,
					real *A0, real *A1, real *A2, real *A3, real *A4, real *A5, real *A6, real *b, real *x,
					int *sz, int *g);

/// @brief Red-black Gauss-Seidel smoother for 7-band matries for Cartesian grid data structure (for coefficinet matrices represented by 7 vectors) \n
///        r = b - A x \n
///        rr += r*r
/// @param [out] x the 3D array to be smoothed
/// @param [in]  A the 3D vector array for coefficient matries
/// @param [in]  b the 3D array for the right-hand side
/// @param [in]  param the relaxation parameter
/// @param [in]  color the color id of which the cells are smoothed
/// @param [in]  offset the offset to control the first color to be updated
/// @param [in]  sz the number of cells in the 3D array (NX, NY, NZ)
/// @param [in]  g the size of the guide cells
	void calc_rr_7_s_(
					real *rr,
					real *A0, real *A1, real *A2, real *A3, real *A4, real *A5, real *A6, real *b, real *x,
					int *sz, int *g);
}

#endif

