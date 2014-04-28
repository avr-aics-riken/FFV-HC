#ifndef BSTL_H
#define BSTL_H

#include "real.h"

extern "C" {
	void bstl_read_cut_1_(
				real* ic0, real* ic1, real* ic2, real* ic3, real* ic4, real* ic5,
				int* cid,
				real* cut,
				int* bid,
				int* sz, int* g);
	void bstl_voxelize_1_(
				real* c0, real* c1, real* c2, real* c3, real* c4, real* c5,
				int* cid,
				int* sz, int* g);
	void bstl_cutoff_1_(
				real* c0, real* c1, real* c2, real* c3, real* c4, real* c5,
				int* cid,
				real* eps,
				int* sz, int* g);
	void bstl_symmetrize_1_(
				real* c0, real* c1, real* c2, real* c3, real* c4, real* c5,
				int* cid,
				int* sz, int* g);
	void bstl_detect_zerocut_1_(
				real* c0, real* c1, real* c2, real* c3, real* c4, real* c5,
				int* cid,
				int* count,
				int* sz, int* g);
	void bstl_fill_holes_1_(
				real* c0, real* c1, real* c2, real* c3, real* c4, real* c5,
				int* cid,
				int* count,
				int *sz, int *g);


	void bstl_read_cut_(
				real* ic0, real* ic1, real* ic2, real* ic3, real* ic4, real* ic5,
				int* cid0, int* cid1, int* cid2, int* cid3, int* cid4, int* cid5,
				real* cut,
				int* bid,
				int* sz, int* g);
	void bstl_voxelize_(
				real* c0, real* c1, real* c2, real* c3, real* c4, real* c5,
				int* cid0, int* cid1, int* cid2, int* cid3, int* cid4, int* cid5,
				int* sz, int* g);
	void bstl_cutoff_(
				real* c0, real* c1, real* c2, real* c3, real* c4, real* c5,
				int* cid0, int* cid1, int* cid2, int* cid3, int* cid4, int* cid5,
				real* eps,
				int* sz, int* g);
	void bstl_symmetrize_(
				real* c0, real* c1, real* c2, real* c3, real* c4, real* c5,
				int* cid0, int* cid1, int* cid2, int* cid3, int* cid4, int* cid5,
				int* sz, int* g);
	void bstl_detect_zerocut_(
				real* c0, real* c1, real* c2, real* c3, real* c4, real* c5,
				int* cid0, int* cid1, int* cid2, int* cid3, int* cid4, int* cid5,
				int* count,
				int* sz, int* g);
	void bstl_fill_holes_(
				real* c0, real* c1, real* c2, real* c3, real* c4, real* c5,
				int* cid0, int* cid1, int* cid2, int* cid3, int* cid4, int* cid5,
				int* count,
				int* bClose,
				int *sz, int *g);
	void bstl_fill_holes_v2_(
				real* c0, real* c1, real* c2, real* c3, real* c4, real* c5,
				int* cid0, int* cid1, int* cid2, int* cid3, int* cid4, int* cid5,
				int* count,
				int* bClose,
				int *sz, int *g);
}

#endif

