/// @file  bplt3d.h
/// @brief Basic subprograms for I/O in PLOT3D format for Cartesian grid data structure

#ifndef BPLT3D_H
#define BPLT3D_H

extern "C" {
	void bplt3d_open_file_(char* filename, int* filenamelength, int* unit);
	void bplt3d_close_file_(int* unit);
	void bplt3d_write_xyz_header_(int* ix, int* jx, int* kx, int* ngrid, int* unit);
	void bplt3d_write_xyz_block_(float* x, float* y, float* z, int* ix, int* jx, int* kx, int* unit);
	void bplt3d_write_func_header_(int* ix, int* jx, int* kx, int* nvar, int* ngrid, int* unit);
	void bplt3d_write_func_block_(float* p, float* ux, float* uy, float* uz, int* ix, int* jx, int* kx, int* unit);
}

#endif

