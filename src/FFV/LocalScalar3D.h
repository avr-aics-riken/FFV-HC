#ifndef LOCALSCALAR3D_H
#define LOCALSCALAR3D_H

#include "BlockManager.h"
#include "BlockScalar3D.h"
#include "Scalar3DUpdater.h"

#include "Scalar3DUpdater1.h"
#include "Scalar3DUpdater2.h"
#include "Scalar3DUpdater3.h"
#include "Scalar3DUpdater4.h"
#include "Scalar3DUpdater5.h"
#include "Scalar3DUpdater6.h"
#include "real.h"
#include "FFVVTKWriter.h"
#include "FFVCommon.h"

template <typename T>
class LocalScalar3D {
public:
	LocalScalar3D(
					BlockManager& blockManager,
					int vc,
					std::string updateMethod,
					int* boundaryType,
					T* boundaryValue,
					int updaterType=1) {
		this->vc = vc;
		if( updaterType == 0 ) {
			this->id = blockManager.setDataClass<Scalar3D<T>, Scalar3DUpdater<T> >(this->vc);
		} else if( updaterType == 1 ) {
			this->id = blockManager.setDataClass<Scalar3D<T>, Scalar3DUpdater1<T> >(this->vc);
		} else if( updaterType == 10 ) {
			this->id = blockManager.setDataClass<Scalar3D<T>, Scalar3DUpdater1<T>, T >(this->vc);
		} else if( updaterType == 2 ) {
			this->id = blockManager.setDataClass<Scalar3D<T>, Scalar3DUpdater2<T> >(this->vc);
		} else if( updaterType == 3 ) {
			this->id = blockManager.setDataClass<Scalar3D<T>, Scalar3DUpdater3<T> >(this->vc);
		} else if( updaterType == 4 ) {
			this->id = blockManager.setDataClass<Scalar3D<T>, Scalar3DUpdater4<T> >(this->vc);
		} else if( updaterType == 5 ) {
			this->id = blockManager.setDataClass<Scalar3D<T>, Scalar3DUpdater5<T> >(this->vc);
		} else if( updaterType == 6 ) {
			this->id = blockManager.setDataClass<Scalar3D<T>, Scalar3DUpdater6<T> >(this->vc);
		} else {
			this->id = blockManager.setDataClass<Scalar3D<T>, Scalar3DUpdater1<T> >(this->vc);
		}

		this->updateMethod = VCUpdateMethod::AtOnce;
		if (updateMethod == "AtOnce") {
			this->updateMethod = VCUpdateMethod::AtOnce;
		} else if (updateMethod == "SeparateXYZ") {
			this->updateMethod = VCUpdateMethod::SeparateXYZ;
		} else if (updateMethod == "SeparateXYZ_SeparateLevelDiff") {
			this->updateMethod = VCUpdateMethod::SeparateXYZ_SeparateLevelDiff;
		}

		int tag = 100;
		blockManager.prepareForVCUpdate(this->id, tag, this->updateMethod);
		this->blockScalar3D = new BlockScalar3D<T> [blockManager.getNumBlock()];
#ifdef _BLOCK_IS_LARGE_
#else
#endif
		for (int n=0; n<blockManager.getNumBlock(); ++n) {
			BlockBase* block	= blockManager.getBlock(n);
			Vec3i size				= block->getSize();
			Vec3r origin			= block->getOrigin();
			Vec3r blockSize		= block->getBlockSize();
			Vec3r cellSize		= block->getCellSize();
			T*    blockData		= this->GetBlockData(block);
			const NeighborInfo* neighborInfo = block->getNeighborInfo();
			this->blockScalar3D[n].InitBoundaryCondition(
					size,
					origin,
					blockSize,
					cellSize,
					vc,
					blockData,
					boundaryType,
					boundaryValue,
					neighborInfo);
		}
	}
	~LocalScalar3D() {
		if( this->blockScalar3D ) {
			delete [] this->blockScalar3D;
			this->blockScalar3D = 0;
		}
	}

private:
	int vc;
	int id;
	VCUpdateMethod::Type updateMethod;
	BlockScalar3D<T>* blockScalar3D;

private:
	T sum_l;
	T max_l;
	T min_l;
	T absmax_l;
	T absmin_l;

	T sum_g;
	T max_g;
	T min_g;
	T absmax_g;
	T absmin_g;

	double totalcells_l;
	double totalcells_g;

public:
	friend void LSSwap(LocalScalar3D<T>& x1, LocalScalar3D<T>& x0) {
		int id_tmp = x1.id;
		x1.id = x0.id;
		x0.id = id_tmp;
	}

public:
	int GetVC() {
		return vc;
	};

	int GetID() {
		return id;
	}

	T GetSum() {
		return sum_g;
	}

	T GetMax() {
		return max_g;
	}

	T GetMin() {
		return min_g;
	}

	T GetAbsMax() {
		return absmax_g;
	}

	T GetAbsMin() {
		return absmin_g;
	}

	T GetSumL() {
		return sum_l;
	}

	T GetMaxL() {
		return max_l;
	}

	T GetMinL() {
		return min_l;
	}

	T GetAbsMaxL() {
		return absmax_l;
	}

	T GetAbsMinL() {
		return absmin_l;
	}

	double GetTotalCellsG() {
		return totalcells_g;
	}

	double GetTotalCellsL() {
		return totalcells_l;
	}

	T* GetBlockData(BlockBase* block) {
		Scalar3D<T>* sdata = dynamic_cast<Scalar3D<T>*>(block->getDataClass(this->id));
		return sdata->getData();
	}

	void UpdateVirtualCells(
										BlockManager& blockManager) {
		switch(this->updateMethod) {
			case VCUpdateMethod::AtOnce:
				blockManager.beginUpdateVC(this->id);
				blockManager.endUpdateVC(this->id);
				break;
			case VCUpdateMethod::SeparateXYZ:
				blockManager.beginUpdateVC_X(this->id);
				blockManager.endUpdateVC_X(this->id);
				blockManager.beginUpdateVC_Y(this->id);
				blockManager.endUpdateVC_Y(this->id);
				blockManager.beginUpdateVC_Z(this->id);
				blockManager.endUpdateVC_Z(this->id);
				break;
			case VCUpdateMethod::SeparateXYZ_SeparateLevelDiff:
				blockManager.updateVC_X_F2C(this->id);
				blockManager.updateVC_Y_F2C(this->id);
				blockManager.updateVC_Z_F2C(this->id);
				blockManager.updateVC_X_Flat(this->id);
				blockManager.updateVC_Y_Flat(this->id);
				blockManager.updateVC_Z_Flat(this->id);
				blockManager.updateVC_X_C2F(this->id);
				blockManager.updateVC_Y_C2F(this->id);
				blockManager.updateVC_Z_C2F(this->id);
				break;
		}
	}

	void ImposeBoundaryCondition(
										BlockManager& blockManager) {
		switch(this->updateMethod) {
			case VCUpdateMethod::AtOnce:
				blockManager.beginUpdateVC(this->id);
				blockManager.endUpdateVC(this->id);
				ImposeBlockBoundaryCondition(blockManager);
				break;
			case VCUpdateMethod::SeparateXYZ:
				blockManager.beginUpdateVC_X(this->id);
				blockManager.endUpdateVC_X(this->id);
				blockManager.beginUpdateVC_Y(this->id);
				blockManager.endUpdateVC_Y(this->id);
				blockManager.beginUpdateVC_Z(this->id);
				blockManager.endUpdateVC_Z(this->id);
				ImposeBlockBoundaryCondition(blockManager);
				break;
			case VCUpdateMethod::SeparateXYZ_SeparateLevelDiff:
				blockManager.updateVC_X_F2C(this->id);
				blockManager.updateVC_Y_F2C(this->id);
				blockManager.updateVC_Z_F2C(this->id);
				blockManager.updateVC_X_Flat(this->id);
				blockManager.updateVC_Y_Flat(this->id);
				blockManager.updateVC_Z_Flat(this->id);
				blockManager.updateVC_X_C2F(this->id);
				blockManager.updateVC_Y_C2F(this->id);
				blockManager.updateVC_Z_C2F(this->id);
				ImposeBlockBoundaryCondition(blockManager);
				break;
		}
	}

	void ImposeBoundaryCondition(
										BlockManager& blockManager,
										LocalScalar3D<T>* plsAp,
										LocalScalar3D<T>* plsAw,
										LocalScalar3D<T>* plsAe,
										LocalScalar3D<T>* plsAs,
										LocalScalar3D<T>* plsAn,
										LocalScalar3D<T>* plsAb,
										LocalScalar3D<T>* plsAt,
										LocalScalar3D<T>* plsb) {
#ifdef _BLOCK_IS_LARGE_
#else
#endif
		for (int n=0; n<blockManager.getNumBlock(); ++n) {
			BlockBase* block = blockManager.getBlock(n);
			real* Ap = plsAp->GetBlockData(block);
			real* Aw = plsAw->GetBlockData(block);
			real* Ae = plsAe->GetBlockData(block);
			real* As = plsAs->GetBlockData(block);
			real* An = plsAn->GetBlockData(block);
			real* Ab = plsAb->GetBlockData(block);
			real* At = plsAt->GetBlockData(block);
			real* b  = plsb ->GetBlockData(block);
			this->blockScalar3D[n].ImposeBoundaryCondition(Ap, Aw, Ae, As, An, Ab, At, b);
		}
	}

	void ResetBoundaryConditionValue(
									BlockManager& blockManager,
									int i_face,
									T boundaryValue) {
#ifdef _BLOCK_IS_LARGE_
#else
#endif
		for (int n=0; n<blockManager.getNumBlock(); ++n) {
			this->blockScalar3D[n].ResetBoundaryConditionValue(i_face, boundaryValue);
		}
	}

private:
	void ImposeBlockBoundaryCondition(BlockManager& blockManager) {
#ifdef _BLOCK_IS_LARGE_
#else
#endif
		for (int n=0; n<blockManager.getNumBlock(); ++n) {
			this->blockScalar3D[n].ImposeBoundaryCondition();
		}
	}

public:
	void WriteDataInVTKFormat(
						const char* path,
						const char* prefix,
						const char* dataname,
						int step,
						int difflevel,
						int maxLevel,
						int minLevel,
						RootGrid* rootGrid,
						BCMOctree* tree,
						Partition* partition,
						Vec3d rootBlockOrigin,
						double rootBlockLength) {
		VtkWriter writer;
		writer.writeScalar<T>(
						this->id,
						this->vc,
						string(path),
						string(prefix),
						string(dataname),
						step,
						difflevel,
						rootGrid,
						tree,
						partition,
						rootBlockOrigin,
						rootBlockLength);

		writer.writeScalar_OverlappingAMR<T>(
						this->id,
						this->vc,
						string(path),
						string(prefix),
						string(dataname),
						step,
						maxLevel,
						minLevel,
						rootGrid,
						tree,
						partition,
						rootBlockOrigin,
						rootBlockLength);
	}

	void Fill(BlockManager& blockManager, T value, T deviation=0.0) {
#ifdef _BLOCK_IS_LARGE_
#else
#endif
		for (int n=0; n<blockManager.getNumBlock(); ++n) {
			BlockBase* block = blockManager.getBlock(n);
			Vec3i size = block->getSize();
			T*      blockData = this->GetBlockData(block);
			int ix = size.x + 2*vc;
			int jx = size.y + 2*vc;
			int kx = size.z + 2*vc;
			for(int k=0; k<kx; k++) {
				for(int j=0; j<jx; j++) {
					for(int i=0; i<ix; i++) {
						int mp = i + ix*(j + jx*k);
						blockData[mp] = RandomNormal(value, deviation);
					}
				}
			}
		}
	}

	void Shift(BlockManager& blockManager, T value) {
		CalcStats(blockManager);
		double sum = GetSum();
		double totalcells_g = GetTotalCellsG();
		std::cout << sum << " " << totalcells_g << std::endl;

#ifdef _BLOCK_IS_LARGE_
#else
#endif
		for (int n=0; n<blockManager.getNumBlock(); ++n) {
			BlockBase* block = blockManager.getBlock(n);
			Vec3i size = block->getSize();
			T*      blockData = this->GetBlockData(block);
			int ix = size.x;
			int jx = size.y;
			int kx = size.z;
			for(int k=vc; k<kx+vc; k++) {
				for(int j=vc; j<jx+vc; j++) {
					for(int i=vc; i<ix+vc; i++) {
						int mp = i + (ix + 2*vc)*(j + (jx + 2*vc)*k);
						blockData[mp] -= sum/totalcells_g;
					}
				}
			}
		}
	}

	void CalcStats(BlockManager& BlockManager) {
	}

	void Dump2(BlockManager& blockManager, const int step, const char* path, const char* prefix, const char* label) {
	}

	void Load2(BlockManager& blockManager, const int step, const char* path, const char* prefix, const char* label) {
	}
};

template <>
void LocalScalar3D<real>::CalcStats(BlockManager& blockManager);

template <>
void LocalScalar3D<real>::Dump2(BlockManager& blockManager, const int step, const char* path, const char* prefix, const char* label);

template <>
void LocalScalar3D<real>::Load2(BlockManager& blockManager, const int step, const char* path, const char* prefix, const char* label);

#endif

