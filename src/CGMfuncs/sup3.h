#ifndef SUP3_H
#define SUP3_H

#include "real.h"

extern "C" {
	void sup3_copy_from_neighbor_(
				real* data_dst,
				int* i1_dst,
				real* data_src,
				int* i1_src,
				int* sz_c,
				int* sz, int* g);

	void sup3_copy_from_neighbor_c2f_x_m_(
				real* data_dst,
				int* i1_dst,
				real* data_src,
				int* i1_src,
				int* sz_c,
				int* sz, int* g);
	void sup3_copy_from_neighbor_c2f_x_p_(
				real* data_dst,
				int* i1_dst,
				real* data_src,
				int* i1_src,
				int* sz_c,
				int* sz, int* g);
	void sup3_copy_from_neighbor_c2f_y_m_(
				real* data_dst,
				int* i1_dst,
				real* data_src,
				int* i1_src,
				int* sz_c,
				int* sz, int* g);
	void sup3_copy_from_neighbor_c2f_y_p_(
				real* data_dst,
				int* i1_dst,
				real* data_src,
				int* i1_src,
				int* sz_c,
				int* sz, int* g);
	void sup3_copy_from_neighbor_c2f_z_m_(
				real* data_dst,
				int* i1_dst,
				real* data_src,
				int* i1_src,
				int* sz_c,
				int* sz, int* g);
	void sup3_copy_from_neighbor_c2f_z_p_(
				real* data_dst,
				int* i1_dst,
				real* data_src,
				int* i1_src,
				int* sz_c,
				int* sz, int* g);

	void sup3_copy_from_neighbor_f2c_(
				real* data_dst,
				int* i1_dst,
				real* data_src,
				int* i1_src,
				int* sz_c,
				int* sz, int* g);

	void sup3_copy_to_buffer_(
				real* buffer,
				real* data_src,
				int* i1_src,
				int* sz_c,
				int* sz, int* g);

	void sup3_copy_to_buffer_c2f_(
				real* buffer,
				real* data_src,
				int* i1_src,
				int* sz_c,
				int* sz, int* g);

	void sup3_copy_to_buffer_f2c_(
				real* buffer,
				real* data_src,
				int* i1_src,
				int* sz_c,
				int* sz, int* g);

	void sup3_copy_from_buffer_(
				real* data_dst,
				int* i1_dst,
				real* buffer,
				int* sz_c,
				int* sz, int* g);

	void sup3_copy_from_buffer_c2f_(
				real* data_dst,
				int* i1_dst,
				real* buffer,
				int* sz_c,
				int* sz, int* g);

	void sup3_copy_from_buffer_f2c_(
				real* data_dst,
				int* i1_dst,
				real* buffer,
				int* sz_c,
				int* sz, int* g);
}

#endif

