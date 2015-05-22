#include "LocalScalar3D.h"

#include "bsf3d.h"
#include "comm.h"
#include "bcut.h"
#include "blas.h"

template <>
void LocalScalar3D<real>::CalcStats(BlockManager& blockManager) {
	BlockBase* block0 = blockManager.getBlock(0);
	Vec3i size = block0->getSize();
	real* pData0 = GetBlockData(block0);

	int m0 = vc + (size.x+2*vc)*(vc + (size.y+2*vc)*vc);
	sum_l = 0.0;
	max_l = pData0[m0];
	min_l = pData0[m0];
	absmax_l = fabs(pData0[m0]);
	absmin_l = fabs(pData0[m0]);
	totalcells_l = 0.0;

	for (int id = 0; id < blockManager.getNumBlock(); ++id) {
		BlockBase* block	= blockManager.getBlock(id);
		Vec3i size				= block->getSize();
		Vec3r origin			= block->getOrigin();
		Vec3r blockSize		= block->getBlockSize();
		Vec3r cellSize		= block->getCellSize();

		int sz[3]		= {size.x, size.y, size.z};
		int g[1]		= {vc};
		real dx			= cellSize.x;

		real* pData = GetBlockData(block);
		real sum_b = 0.0;
		real max_b = 0.0;
		real min_b = 0.0;
		real absmax_b = 0.0;
		real absmin_b = 0.0;
		sf3d_calc_stats_(
				&sum_b,
				&max_b,
				&min_b,
				&absmax_b,
				&absmin_b,
				pData,
				sz, g);

		sum_l += sum_b;
		if( max_b > max_l ) {
			max_l = max_b;
		}
		if( min_b < min_l ) {
			min_l = min_b;
		}
		if( absmax_b > absmax_l ) {
			absmax_l = absmax_b;
		}
		if( absmin_b < absmin_l ) {
			absmin_l = absmin_b;
		}

		totalcells_l += size.x*size.y*size.z;
	}

	comm_sum_(&sum_g, &sum_l);
	comm_max_(&max_g, &max_l);
	comm_min_(&min_g, &min_l);
	comm_max_(&absmax_g, &absmax_l);
	comm_min_(&absmin_g, &absmin_l);
	comm_sum_(&totalcells_g, &totalcells_l);
}

template <>
void LocalScalar3D<real>::Dump(BlockManager& blockManager, const int step, const char* label) {
	ImposeBoundaryCondition(blockManager);
	MPI::Intracomm comm = blockManager.getCommunicator();

	ostringstream ossFileNameTime;
	ossFileNameTime << "./BIN/";
	mkdir(ossFileNameTime.str().c_str(), 0755);

#ifdef _BLOCK_IS_LARGE_
#else
#endif
	for (int id = 0; id < blockManager.getNumBlock(); ++id) {
		BlockBase* block = blockManager.getBlock(id);

		Vec3i size = block->getSize();
		Vec3r origin = block->getOrigin();
		Vec3r blockSize = block->getBlockSize();
		Vec3r cellSize = block->getCellSize();
		int level = block->getLevel();

		ostringstream ossFileName;
		ossFileName << "./BIN/";
		ossFileName << "dump-";
		ossFileName << label;
		ossFileName << "-";
		ossFileName.width(5);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << comm.Get_rank();
		ossFileName << "-";
		ossFileName.width(5);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << id;
		ossFileName << "-";
		ossFileName.width(10);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << step;
		ossFileName << ".bin";

		int cx = size.x + 2*vc;
		int cy = size.y + 2*vc;
		int cz = size.z + 2*vc;
		int iNE = 1;

		real* pData = GetBlockData(block);

		ofstream ofs;
		ofs.open(ossFileName.str().c_str(), ios::out | ios::binary);
		ofs.write((char*)&size.x, sizeof(int));
		ofs.write((char*)&size.y, sizeof(int));
		ofs.write((char*)&size.z, sizeof(int));
		ofs.write((char*)&vc    , sizeof(int));
		ofs.write((char*)&iNE   , sizeof(int));
		ofs.write((char*)pData  , sizeof(real)*cx*cy*cz);
		ofs.close();
	}
}

template <>
void LocalScalar3D<real>::Load(BlockManager& blockManager, const int step, const char* label) {
	MPI::Intracomm comm = blockManager.getCommunicator();

#ifdef _BLOCK_IS_LARGE_
#else
#endif
	for (int id = 0; id < blockManager.getNumBlock(); ++id) {
		BlockBase* block = blockManager.getBlock(id);

		Vec3i size = block->getSize();
		Vec3r origin = block->getOrigin();
		Vec3r blockSize = block->getBlockSize();
		Vec3r cellSize = block->getCellSize();
		int level = block->getLevel();

		ostringstream ossFileName;
		ossFileName << "./BIN/";
		ossFileName << "dump-";
		ossFileName << label;
		ossFileName << "-";
		ossFileName.width(5);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << comm.Get_rank();
		ossFileName << "-";
		ossFileName.width(5);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << id;
		ossFileName << "-";
		ossFileName.width(10);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << step;
		ossFileName << ".bin";

		int nx = 0;
		int ny = 0;
		int nz = 0;
		int nv = 0;
		int ne = 0;

		real* pData = GetBlockData(block);

		ifstream ifs;
		ifs.open(ossFileName.str().c_str(), ios::in | ios::binary);
		ifs.read((char*)&nx, sizeof(int));
		ifs.read((char*)&ny, sizeof(int));
		ifs.read((char*)&nz, sizeof(int));
		ifs.read((char*)&nv, sizeof(int));
		ifs.read((char*)&ne, sizeof(int));
		if( nx == size.x && ny == size.y && nz == size.z && nv == vc && ne == 1 ) {
			int cx = nx + 2*nv;
			int cy = ny + 2*nv;
			int cz = nz + 2*nv;
			ifs.read((char*)pData, sizeof(real)*cx*cy*cz);
		} else if( 2*nx == size.x && 2*ny == size.y && 2*nz == size.z && nv == vc && ne == 1 ) {
			int cx = nx + 2*nv;
			int cy = ny + 2*nv;
			int cz = nz + 2*nv;
			real *pDataS = new real [cx*cy*cz];
			ifs.read((char*)pDataS, sizeof(real)*cx*cy*cz);

			int sz[3] = {2*nx, 2*ny, 2*nz};
			sf3d_copy_x2_(
					(real*)pData,
					(real*)pDataS,
					(int*)sz,
					(int*)&vc);

			delete [] pDataS;
		} else {
			Exit(0);
		}
		ifs.close();
	}

	ImposeBoundaryCondition(blockManager);
}

template <>
void LocalScalar3D<real>::Dump2(BlockManager& blockManager, const int step, const char* label) {
	ImposeBoundaryCondition(blockManager);
	MPI::Intracomm comm = blockManager.getCommunicator();

	std::string path = g_pFFVConfig->RestartOutputPath;
	std::string prefix = g_pFFVConfig->RestartPrefix;

	ostringstream ossFileNameTime;
	ossFileNameTime << path;
	mkdir(ossFileNameTime.str().c_str(), 0755);

	ostringstream ossFileName;
	ossFileName << path;
	ossFileName << "/";
	ossFileName << prefix;
	ossFileName << label;
	ossFileName << "-";
	ossFileName.width(5);
	ossFileName.setf(ios::fixed);
	ossFileName.fill('0');
	ossFileName << comm.Get_rank();
	ossFileName << "-";
	ossFileName.width(10);
	ossFileName.setf(ios::fixed);
	ossFileName.fill('0');
	ossFileName << step;
	ossFileName << ".bin";

	BlockBase* block = blockManager.getBlock(0);
	Vec3i size = block->getSize();
	int cx = size.x + 2*vc;
	int cy = size.y + 2*vc;
	int cz = size.z + 2*vc;
	int iNE = 1;
	int iNB = blockManager.getNumBlock();

	ofstream ofs;
	ofs.open(ossFileName.str().c_str(), ios::out | ios::binary);
	ofs.write((char*)&size.x, sizeof(int));
	ofs.write((char*)&size.y, sizeof(int));
	ofs.write((char*)&size.z, sizeof(int));
	ofs.write((char*)&vc    , sizeof(int));
	ofs.write((char*)&iNE   , sizeof(int));
	ofs.write((char*)&iNB   , sizeof(int));

#ifdef _BLOCK_IS_LARGE_
#else
#endif
	for (int id = 0; id < blockManager.getNumBlock(); ++id) {
		block = blockManager.getBlock(id);

		real* pData = GetBlockData(block);

		ofs.write((char*)pData  , sizeof(real)*cx*cy*cz);
	}

	ofs.close();
}

template <>
void LocalScalar3D<real>::Load2(BlockManager& blockManager, const int step, const char* label) {
	MPI::Intracomm comm = blockManager.getCommunicator();

	std::string path = g_pFFVConfig->RestartInputPath;
	std::string prefix = g_pFFVConfig->RestartPrefix;

	ostringstream ossFileName;
	ossFileName << path;
	ossFileName << "/";
	ossFileName << prefix;
	ossFileName << label;
	ossFileName << "-";
	ossFileName.width(5);
	ossFileName.setf(ios::fixed);
	ossFileName.fill('0');
	ossFileName << comm.Get_rank();
	ossFileName << "-";
	ossFileName.width(10);
	ossFileName.setf(ios::fixed);
	ossFileName.fill('0');
	ossFileName << step;
	ossFileName << ".bin";

	int nx = 0;
	int ny = 0;
	int nz = 0;
	int nv = 0;
	int ne = 0;
	int nb = 0;

	ifstream ifs;
	ifs.open(ossFileName.str().c_str(), ios::in | ios::binary);
	ifs.read((char*)&nx, sizeof(int));
	ifs.read((char*)&ny, sizeof(int));
	ifs.read((char*)&nz, sizeof(int));
	ifs.read((char*)&nv, sizeof(int));
	ifs.read((char*)&ne, sizeof(int));
	ifs.read((char*)&nb, sizeof(int));

	int cx = nx + 2*nv;
	int cy = ny + 2*nv;
	int cz = nz + 2*nv;

	BlockBase* block = blockManager.getBlock(0);
	Vec3i size = block->getSize();

	if( nx == size.x && ny == size.y && nz == size.z && nv == vc && ne == 1 && nb == blockManager.getNumBlock() ) {
		for (int id = 0; id < blockManager.getNumBlock(); ++id) {
			block = blockManager.getBlock(id);

			real* pData = GetBlockData(block);

			ifs.read((char*)pData, sizeof(real)*cx*cy*cz);
		}
	} else if( 2*nx == size.x && 2*ny == size.y && 2*nz == size.z && nv == vc && ne == 1 && nb == blockManager.getNumBlock() ) {
		real *pDataS = new real [cx*cy*cz];
		for (int id = 0; id < blockManager.getNumBlock(); ++id) {
			block = blockManager.getBlock(id);

			real* pData = GetBlockData(block);

			ifs.read((char*)pDataS, sizeof(real)*cx*cy*cz);

			int sz[3] = {2*nx, 2*ny, 2*nz};
			sf3d_copy_x2_(
					(real*)pData,
					(real*)pDataS,
					(int*)sz,
					(int*)&vc);
		}
		delete [] pDataS;
	} else {
		Exit(0);
	}

	ifs.close();

	ImposeBoundaryCondition(blockManager);
}

template <>
void LocalScalar3D<real>::Dump3(BlockManager& blockManager, const int step, const char* label, Partition* partition, int myrank) {
	ImposeBoundaryCondition(blockManager);
	MPI::Intracomm comm = blockManager.getCommunicator();

	ostringstream ossFileNameTime;
	ossFileNameTime << "./BIN/";
	mkdir(ossFileNameTime.str().c_str(), 0755);

	ossFileNameTime.width(10);
	ossFileNameTime.setf(ios::fixed);
	ossFileNameTime.fill('0');
	ossFileNameTime << step;
	mkdir(ossFileNameTime.str().c_str(), 0755);

#ifdef _BLOCK_IS_LARGE_
#else
#endif
	for (int id = 0; id < blockManager.getNumBlock(); ++id) {
		BlockBase* block = blockManager.getBlock(id);

		Vec3i size = block->getSize();
		Vec3r origin = block->getOrigin();
		Vec3r blockSize = block->getBlockSize();
		Vec3r cellSize = block->getCellSize();
		int level = block->getLevel();

		ostringstream ossFileName;
		ossFileName << "./BIN/";
		ossFileName.width(10);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << step;
		ossFileName << "/";
		ossFileName << "dump-";
		ossFileName << label;
		ossFileName << "-";
		ossFileName.width(5);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << id + partition->getStart(myrank);
		ossFileName << "-";
		ossFileName.width(10);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << step;
		ossFileName << ".bin";

		int cx = size.x + 2*vc;
		int cy = size.y + 2*vc;
		int cz = size.z + 2*vc;
		int iNE = 1;
		int iNB = 1;

		real* pData = GetBlockData(block);

		ofstream ofs;
		ofs.open(ossFileName.str().c_str(), ios::out | ios::binary);
		ofs.write((char*)&size.x, sizeof(int));
		ofs.write((char*)&size.y, sizeof(int));
		ofs.write((char*)&size.z, sizeof(int));
		ofs.write((char*)&vc    , sizeof(int));
		ofs.write((char*)&iNE   , sizeof(int));
		ofs.write((char*)&iNB   , sizeof(int));
		ofs.write((char*)pData  , sizeof(real)*cx*cy*cz);
		ofs.close();
	}
}

template <>
void LocalScalar3D<real>::Load3(BlockManager& blockManager, const int step, const char* label, Partition* partition, int myrank) {
	MPI::Intracomm comm = blockManager.getCommunicator();

#ifdef _BLOCK_IS_LARGE_
#else
#endif
	for (int id = 0; id < blockManager.getNumBlock(); ++id) {
		ostringstream ossFileName;
		ossFileName << "./BIN/";
		ossFileName.width(10);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << step;
		ossFileName << "/";
		ossFileName << "dump-";
		ossFileName << label;
		ossFileName << "-";
		ossFileName.width(5);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << id + partition->getStart(myrank);
		ossFileName << "-";
		ossFileName.width(10);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << step;
		ossFileName << ".bin";

		int nx = 0;
		int ny = 0;
		int nz = 0;
		int nv = 0;
		int ne = 0;
		int nb = 0;

		ifstream ifs;
		ifs.open(ossFileName.str().c_str(), ios::in | ios::binary);
		ifs.read((char*)&nx, sizeof(int));
		ifs.read((char*)&ny, sizeof(int));
		ifs.read((char*)&nz, sizeof(int));
		ifs.read((char*)&nv, sizeof(int));
		ifs.read((char*)&ne, sizeof(int));
		ifs.read((char*)&nb, sizeof(int));

		int cx = nx + 2*nv;
		int cy = ny + 2*nv;
		int cz = nz + 2*nv;

		BlockBase* block = blockManager.getBlock(id);
		Vec3i size = block->getSize();
		real* pData = GetBlockData(block);
		real *pDataS = new real [cx*cy*cz];
		if( nx == size.x && ny == size.y && nz == size.z && nv == vc && ne == 1 && nb == 1 ) {
			ifs.read((char*)pData, sizeof(real)*cx*cy*cz);
		} else if( 2*nx == size.x && 2*ny == size.y && 2*nz == size.z && nv == vc && ne == 1 && nb == 1 ) {
			ifs.read((char*)pDataS, sizeof(real)*cx*cy*cz);
			int sz[3] = {2*nx, 2*ny, 2*nz};
			sf3d_copy_x2_(
					(real*)pData,
					(real*)pDataS,
					(int*)sz, (int*)&vc);
		} else {
			Exit(0);
		}
		ifs.close();
		delete [] pDataS;
	}
	ImposeBoundaryCondition(blockManager);
}

template <>
real LocalScalar3D<real>::GetValue(BlockManager& blockManager, real xr, real yr, real zr) {
	real vr_local = 0.0;
	int  flag_local = 0;
	for (int n=0; n<blockManager.getNumBlock(); ++n) {
		BlockBase* block = blockManager.getBlock(n);
		Vec3i size = block->getSize();
		Vec3r origin = block->getOrigin();
		Vec3r blockSize = block->getBlockSize();
		Vec3r cellSize = block->getCellSize();

		int sz[3] = {size.x, size.y, size.z};
		int g[1] = {vc};
		real dx = cellSize.x;
		real org[3] = {origin.x, origin.y, origin.z};

		real* pData = this->GetBlockData(block);

		real vr_block = 0.0;
		int  flag_block = 0;
		bcut_get_value_at_referencepoint_(&vr_block, &flag_block, pData, &xr, &yr, &zr, &dx, org, sz, g);

		if( flag_block == 1 ) {
			vr_local = vr_block;
			flag_local = 1;
		}
	}

	int flag_global = 0;
	MPI_Allreduce(&flag_local, &flag_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	if( flag_global != 1 ) {
		Exit(1);
	}

	real vr_global = 0.0;
	MPI_Allreduce(&vr_local, &vr_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);

	return vr_global;
}

template <>
void LocalScalar3D<real>::AddValue(BlockManager& blockManager, real vr) {
#ifdef _BLOCK_IS_LARGE_
#else
#endif
	for (int n=0; n<blockManager.getNumBlock(); ++n) {
		BlockBase* block = blockManager.getBlock(n);
		Vec3i size = block->getSize();
		Vec3r origin = block->getOrigin();
		Vec3r blockSize = block->getBlockSize();
		Vec3r cellSize = block->getCellSize();

		int sz[3] = {size.x, size.y, size.z};
		int g[1] = {vc};
		real dx = cellSize.x;
		real org[3] = {origin.x, origin.y, origin.z};

		real* pData = this->GetBlockData(block);
		adda_(pData, &vr, sz, g);
	}
	ImposeBoundaryCondition(blockManager);
}

