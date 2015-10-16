#ifndef FFVTA_H
#define FFVTA_H

#include <cmath>
#include <cfloat>

#include "BlockManager.h"
#include "LocalScalar3D.h"

#include "real.h"
#include "blas.h"

class FFVTA {
	public:
		FFVTA() {
			count = 0;
		}

		~FFVTA() {
		}

	private:
		int count;

	public:
		void Init() {
			count = 0;
		}

		void Update(
				BlockManager& blockManager,
				LocalScalar3D<real>* plsxa,
				LocalScalar3D<real>* plsx) {
			int vc = plsx->GetVC();
			real beta = 1.0/(1.0 + count);
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
#endif
			for (int n=0; n<blockManager.getNumBlock(); ++n) {
				BlockBase* block	= blockManager.getBlock(n);
				Vec3i blockSize		= block->getSize();
				Vec3r cellSize		= block->getCellSize();
				int sz[3]					= {blockSize.x, blockSize.y, blockSize.z};

				real* xa = plsxa->GetBlockData(block);
				real* x  = plsx->GetBlockData(block);

				avew_(xa, x, &beta, sz, &vc);
			}
			plsxa->ImposeBoundaryCondition(blockManager);
			count++;
		}

		void Update(
				BlockManager& blockManager,
				LocalScalar3D<real>* plsxya,
				LocalScalar3D<real>* plsx,
				LocalScalar3D<real>* plsy) {
			int vc = plsx->GetVC();
			real beta = 1.0/(1.0 + count);
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
#endif
			for (int n=0; n<blockManager.getNumBlock(); ++n) {
				BlockBase* block	= blockManager.getBlock(n);
				Vec3i blockSize		= block->getSize();
				Vec3r cellSize		= block->getCellSize();
				int sz[3]					= {blockSize.x, blockSize.y, blockSize.z};

				real* xya = plsxya->GetBlockData(block);
				real* x   = plsx  ->GetBlockData(block);
				real* y   = plsy  ->GetBlockData(block);

				avew_2_(xya, x, y, &beta, sz, &vc);
			}
			plsxya->ImposeBoundaryCondition(blockManager);
			count++;
		}
};

#endif

