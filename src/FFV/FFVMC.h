#ifndef FFVMC3D_H
#define FFVMC3D_H

#include <vector>

#include "BCMTools.h"
#include "BlockManager.h"
#include "Scalar3D.h"
#include "BCMOctree.h"
#include "Partition.h"

typedef struct _Vertex {
	double x;
	double y;
	double z;
}Vertex;

class FFVMC {
public:
	FFVMC()
		: blockManager(BlockManager::getInstance()),
			comm(blockManager.getCommunicator()) {
		int myrank = comm.Get_rank();
		for(int id=0; id<blockManager.getNumBlock(); ++id) {
			BlockBase* block = blockManager.getBlock(id);
		}

		if( myrank == 0 ) {
		} else {
		}
	}
	virtual ~FFVMC() {
	}

private:
	BlockManager& blockManager;
	const MPI::Intracomm& comm;

public:
	void WriteContour(
				int dataClassID,
				int vc,
				const std::string path,
				const std::string prefix,
				const std::string name,
				int step,
				int maxLevel,
				int minLevel,
				RootGrid* rootGrid,
				BCMOctree* tree,
				Partition* partition,
				Vec3r rootOrigin,
				double rootLength,
				double threshold) {

	}

private:



private:
	int PNX;
	int PNY;
	int PNZ;
	double X0;
	double Y0;
	double Z0;
	double DX;
	float *pPointData;

	std::vector<Vertex> vVertexList;

private:
	void AllocMemoryForPointData(int dataClassID);
	void CalcPointData(int dataClassID);
	void DetectTriangles(int dataClassID, double threshold);

	void PrintPVTP(int step, const char* label);
	void PrintVTP(int step, const char* label);
};

#endif 

