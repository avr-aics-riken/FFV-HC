/// @file  comm.h
/// @brief Wrapper functions for inter-process communications for Cartesian grid data structure

#ifndef COMM_H
#define COMM_H

#include "real.h"

extern "C" {
/// @brief Compute the sum of the values for all the processes
/// @param [out] a_sum the sum of the values for all the processes
/// @param [in]  a the value to be summed
	void comm_sum_(real *a_sum, real *a);

/// @brief Compute the maximum of the values for all the processes
/// @param [out] a_max the maximum of the values for all the processes
/// @param [in]  a the value 
	void comm_max_(real *a_max, real *a);

/// @brief Compute the minimum of the values for all the processes
/// @param [out] a_max the minimum of the values for all the processes
/// @param [in]  a the value 
	void comm_min_(real *a_min, real *a);

/// @brief Compute the sum of the inner product of two vectors of all the processes
/// @param [out] xy the sum of the inner product of x and y of all the processes
/// @param [in]  x the 3D array for the inner product
/// @param [in]  y the 3D array for the inner product
/// @param [in]  sz the number of cells in the 3D array (NX, NY, NZ)
/// @param [in]  g the size of the guide cells
	void comm_dot_(real *xy, real *y, real *x, int *sz, int *g);

/// @brief Exchange the values of the guide cells by blocking send/recv
/// @param [out] x the 3D array to be exchanged
/// @param [in]  sz the number of cells in the 3D array (NX, NY, NZ)
/// @param [in]  g the size of the guide cells
/// @param [in]  mx
/// @param [in]  sx1 the 3D array for the send-buffer of 1-face
/// @param [in]  sx3 the 3D array for the send-buffer of 3-face
/// @param [in]  sx2 the 3D array for the send-buffer of 2-face
/// @param [in]  sx4 the 3D array for the send-buffer of 4-face
/// @param [in]  sx5 the 3D array for the send-buffer of 5-face
/// @param [in]  sx6 the 3D array for the send-buffer of 6-face
/// @param [in]  rx1 the 3D array for the recv-buffer of 1-face
/// @param [in]  rx3 the 3D array for the recv-buffer of 3-face
/// @param [in]  rx2 the 3D array for the recv-buffer of 2-face
/// @param [in]  rx4 the 3D array for the recv-buffer of 4-face
/// @param [in]  rx5 the 3D array for the recv-buffer of 5-face
/// @param [in]  rx6 the 3D array for the recv-buffer of 6-face
/// @param [in]  node the array of the ranks, with which the values are exchanged
	void comm_band_cells_(
					real *x,
					int *sz, int *g,
					int *mx,
					real *sx1, real *sx3, real *sx2, real *sx4, real *sx5, real *sx6,
					real *rx1, real *rx3, real *rx2, real *rx4, real *rx5, real *rx6,
					int *node);

/// @brief Exchange the values of the guide cells by non-blocking send/recv
/// @param [out] x the 3D array to be exchanged
/// @param [in]  sz the number of cells in the 3D array (NX, NY, NZ)
/// @param [in]  g the size of the guide cells
/// @param [in]  mx
/// @param [in]  sx1 the 3D array for the send-buffer of 1-face
/// @param [in]  sx3 the 3D array for the send-buffer of 3-face
/// @param [in]  sx2 the 3D array for the send-buffer of 2-face
/// @param [in]  sx4 the 3D array for the send-buffer of 4-face
/// @param [in]  sx5 the 3D array for the send-buffer of 5-face
/// @param [in]  sx6 the 3D array for the send-buffer of 6-face
/// @param [in]  rx1 the 3D array for the recv-buffer of 1-face
/// @param [in]  rx3 the 3D array for the recv-buffer of 3-face
/// @param [in]  rx2 the 3D array for the recv-buffer of 2-face
/// @param [in]  rx4 the 3D array for the recv-buffer of 4-face
/// @param [in]  rx5 the 3D array for the recv-buffer of 5-face
/// @param [in]  rx6 the 3D array for the recv-buffer of 6-face
/// @param [in]  node the array of the ranks, with which the values are exchanged
	void comm_band_cells_nb_(
					real *x,
					int *sz, int *g,
					int *mx,
					real *sx1, real *sx3, real *sx2, real *sx4, real *sx5, real *sx6,
					real *rx1, real *rx3, real *rx2, real *rx4, real *rx5, real *rx6,
					int *node);
}

#endif

