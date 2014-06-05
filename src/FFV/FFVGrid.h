#ifndef FFVGRID_H
#define FFVGRID_H

#include "LocalScalar3D.h"

#include <BlockManager.h>
#include <BCMOctree.h>
#include <BCMPolylib.h>

class FFVGrid {
public:
	FFVGrid();
	~FFVGrid();

private:
	BlockManager& blockManager;
	int rank;
	int nProcesses;
	int nBlocks;

	PolylibNS::BCMPolylib*  pl;
	RootGrid*               rootGrid;
	Divider*                divider;
	BCMOctree::Ordering     ordering;
	BCMOctree*              tree;
	Partition*              partition;
	int                     levelMin;
	int                     levelMax;

	LocalScalar3D<real>     *plsCut0;
	LocalScalar3D<real>     *plsCut1;
	LocalScalar3D<real>     *plsCut2;
	LocalScalar3D<real>     *plsCut3;
	LocalScalar3D<real>     *plsCut4;
	LocalScalar3D<real>     *plsCut5;
	LocalScalar3D<int>      *plsCutId0;
	LocalScalar3D<int>      *plsCutId1;
	LocalScalar3D<int>      *plsCutId2;
	LocalScalar3D<int>      *plsCutId3;
	LocalScalar3D<int>      *plsCutId4;
	LocalScalar3D<int>      *plsCutId5;

	int                     *pNormalN;
	real                    **pNormalX;
	real                    **pNormalY;
	real                    **pNormalZ;
	LocalScalar3D<int>      *plsNormalIndex0;
	LocalScalar3D<int>      *plsNormalIndex1;
	LocalScalar3D<int>      *plsNormalIndex2;
	LocalScalar3D<int>      *plsNormalIndex3;
	LocalScalar3D<int>      *plsNormalIndex4;
	LocalScalar3D<int>      *plsNormalIndex5;

	LocalScalar3D<int>      *plsBoundaryCellFlag;
	LocalScalar3D<int>      *plsRegionId;
	LocalScalar3D<int>      *plsRegionIdTmp;

public:
	void Init();

private:
	void LoadPolygon();
	void DivideDomain();
	void DistributePolygon();
	void InitVarsForCutInfo();
	void CalcCutInfo();
	void DetectBoundaryCells();
	void DefineRegions();
	void FillRegion(double xs, double ys, double zs, int rid);
	void FillRegion2(double xs, double ys, double zs, int rid);
	void FillRemainingRegion(int rid);

public:
	void WriteCutRaw(const char* path, const char* prefix);

private:
	void WritePolygon(std::ofstream& ofs, float* pv);
};

#endif

