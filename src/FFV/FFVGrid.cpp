#include "FFVGrid.h"
#include "FFVCommon.h"
#include "FFVGlobalVars.h"
#include "BlockBoundingBox.h"

#include <bstl.h>
#include <SimpleDivider.h>
#include <Cutlib.h>
#include <CutInfo/CutInfo.h>
#include <CutInfo/CutNormalArray.h>
#include <GridAccessor/Cell.h>

#include <limits.h>

FFVGrid::FFVGrid() : blockManager(BlockManager::getInstance()) {
}

FFVGrid::~FFVGrid() {
}

void FFVGrid::Init() {
	this->rank        = GetMPIRank();
	this->nProcesses  = GetNumProcesses();
	this->nBlocks     = 0;

	LoadPolygon();
	DivideDomain();
	DistributePolygon();
	InitVarsForCutInfo();
	CalcCutInfo();
	ClasifyCell();
	FillRegion();

//	if(g_pFFVConfig->WritePolygon()) {
		WriteCutRaw("./STL", "data");
//	}

	plsPhaseId->WriteDataInVTKFormat(
									g_pFFVConfig->OutputDataFormatOptionVTKPath.c_str(),
									g_pFFVConfig->OutputDataFormatOptionVTKPrefix.c_str(),
									"phase",
									0,
									levelMax - levelMin,
									levelMax,
									levelMin,
									rootGrid,
									tree,
									partition,
									g_pFFVConfig->RootBlockOrigin,
									g_pFFVConfig->RootBlockLength);

	plsCellId->WriteDataInVTKFormat(
									g_pFFVConfig->OutputDataFormatOptionVTKPath.c_str(),
									g_pFFVConfig->OutputDataFormatOptionVTKPrefix.c_str(),
									"cell",
									0,
									levelMax - levelMin,
									levelMax,
									levelMin,
									rootGrid,
									tree,
									partition,
									g_pFFVConfig->RootBlockOrigin,
									g_pFFVConfig->RootBlockLength);
}

void FFVGrid::LoadPolygon() {
g_pFFVPerfMonitor->Start(ffv_tm_Init_LoadPolygon, 0, 0, true);
{
	pl = new PolylibNS::BCMPolylib;
	struct stat st;
	int ret = stat(g_pFFVConfig->PolylibConfig.c_str(), &st);
	if( ret == 0 ) {
		pl->load(g_pFFVConfig->PolylibConfig);
	} else {
		Exit(EX_FAILURE);
	}
}
g_pFFVPerfMonitor->Stop(ffv_tm_Init_LoadPolygon);
}

void FFVGrid::DivideDomain() {
g_pFFVPerfMonitor->Start(ffv_tm_Init_DivideDomain, 0, 0, true);
{
	if(rank == 0) {
		rootGrid = new RootGrid(g_pFFVConfig->RootBlockGrid);
		if( g_pFFVConfig->RootBlockPeriodicX ) {
			rootGrid->setPeriodicX();
		}
		if( g_pFFVConfig->RootBlockPeriodicY ) {
			rootGrid->setPeriodicY();
		}
		if( g_pFFVConfig->RootBlockPeriodicZ ) {
			rootGrid->setPeriodicZ();
		}

		if (!strcasecmp(g_pFFVConfig->TreeType.c_str(), "flat") || !strcasecmp(g_pFFVConfig->TreeType.c_str(), "uniform")) {
			divider = new FlatDivider(
												rootGrid,
												g_pFFVConfig->TreeMaxLevel);
		} else if (!strcasecmp(g_pFFVConfig->TreeType.c_str(), "simple")) {
			divider = new SimpleDivider(
												rootGrid,
												g_pFFVConfig->TreeMinLevel,
												g_pFFVConfig->TreeMaxLevel);
		} else if (!strcasecmp(g_pFFVConfig->TreeType.c_str(), "polygon")) {
			divider = new PolygonBBoxDivider(
												g_pFFVConfig->RootBlockOrigin,
												g_pFFVConfig->RootBlockLength,
												rootGrid,
												g_pFFVConfig->TreeMinLevel,
												pl,
												g_pFFVConfig->PolygonGroupList,
												g_pFFVConfig->BoundingBoxList,
												(double)((double)g_pFFVConfig->LeafBlockNumberOfMarginalCells/(double)g_pFFVConfig->LeafBlockNumberOfCells));
/*
		} else if(!strcasecmp(g_pFFVConfig->TreeType.c_str(), "sphere_old") || !strcasecmp(g_pFFVConfig->TreeType.c_str(), "sphere2")) {
			divider = new SphereDivider2(
												rootGrid,
												g_pFFVConfig->TreeMinLevel,
												g_pFFVConfig->TreeMaxLevel,
												g_pFFVConfig->TreeDividerCenter.x,
												g_pFFVConfig->TreeDividerCenter.y,
												g_pFFVConfig->TreeDividerCenter.z,
												g_pFFVConfig->TreeDividerRadius,
												g_pFFVConfig->TreeDividerDeltaR,
												g_pFFVConfig->TreeDividerHollow);
		} else if (!strcasecmp(g_pFFVConfig->TreeType.c_str(), "sphere")) {
			divider = new SphereDivider3(
												g_pFFVConfig->RootBlockOrigin,
												g_pFFVConfig->RootBlockLength,
												rootGrid,
												g_pFFVConfig->TreeMinLevel,
												pl,
												g_pFFVConfig->PolygonGroupList,
												g_pFFVConfig->BoundingBoxList,
												g_pFFVConfig->SphericalBoxList,
												(double)g_pFFVConfig->LeafBlockNumberOfVirtualCells/g_pFFVConfig->LeafBlockNumberOfCells);
*/
		} else {
			exit(EX_READ_CONFIG);
		}
	}
}
g_pFFVPerfMonitor->Stop(ffv_tm_Init_DivideDomain);

g_pFFVPerfMonitor->Start(ffv_tm_Init_OrderBlock, 0, 0, true);
{
	if(rank == 0) {
		if( g_pFFVConfig->TuningBlockOrdering == "Z") {
			ordering = BCMOctree::Z;
		} else if( g_pFFVConfig->TuningBlockOrdering == "Hilbert" ) {
			ordering = BCMOctree::HILBERT;
		} else if( g_pFFVConfig->TuningBlockOrdering == "random" ) {
			ordering = BCMOctree::RANDOM;
		} else if( g_pFFVConfig->TuningBlockOrdering == "PedigreeList" ) {
			ordering = BCMOctree::PEDIGREELIST;
		} else {
			Exit(EX_READ_CONFIG);
		}
	}
}
g_pFFVPerfMonitor->Stop(ffv_tm_Init_OrderBlock);

g_pFFVPerfMonitor->Start(ffv_tm_Init_CommOctree, 0, 0, true);
{
	if(rank == 0) {
		tree = new BCMOctree(rootGrid, divider, ordering);
		tree->broadcast();
	} else {
		tree = BCMOctree::ReceiveFromMaster();
	}
	nBlocks = tree->getNumLeafNode();
}
g_pFFVPerfMonitor->Stop(ffv_tm_Init_CommOctree);

g_pFFVPerfMonitor->Start(ffv_tm_Init_PartitionBlock, 0, 0, true);
{
	partition = new Partition(nProcesses, nBlocks);
	if( rank == 0 ) {
		for(int n=0; n<nProcesses; n++) {
			if( n==0 ) {
				PrintLog(2, "%-20s : #%05d [%04d:%04d] (%04d)", 
										"Partitions",
										n,
										partition->getStart(n), partition->getEnd(n)-1,
										partition->getEnd(n) - partition->getStart(n));
			} else {
				PrintLog(2, "%-20s : #%05d [%04d:%04d] (%04d)", 
										"          ",
										n,
										partition->getStart(n), partition->getEnd(n)-1,
										partition->getEnd(n) - partition->getStart(n));
			}
		}
	}
}
g_pFFVPerfMonitor->Stop(ffv_tm_Init_PartitionBlock);

g_pFFVPerfMonitor->Start(ffv_tm_Init_RegisterBlock, 0, 0, true);
{
	int nCells = g_pFFVConfig->LeafBlockNumberOfCells;
	Vec3r rootOrigin = g_pFFVConfig->RootBlockOrigin;
	double rootLength = g_pFFVConfig->RootBlockLength;
	Vec3i size(nCells, nCells, nCells);
	std::vector<Node*>& leafNodeArray = tree->getLeafNodeArray();
	int level_min = INT_MAX;
	int level_max = 0;
	for(int id = partition->getStart(rank); id<partition->getEnd(rank); id++) {
		Node* node = leafNodeArray[id];
		int level = node->getLevel();
		level_min = level_min < level ? level_min : level;
		level_max = level_max > level ? level_max : level;
		Vec3r origin = tree->getOrigin(node) * rootLength + rootOrigin;
		Vec3r blockSize = node->getBlockSize() * rootLength;
		NeighborInfo* neighborInfo = tree->makeNeighborInfo(node, partition); 
		for (int i = 0; i < NUM_FACE; ++i) {
			Face face = Face(i);
			switch(face) {
				case X_M:
				case X_P:
					if( g_pFFVConfig->RootBlockPeriodicX == true ) {
						continue;
					}
					break;
				case Y_M:
				case Y_P:
					if( g_pFFVConfig->RootBlockPeriodicY == true ) {
						continue;
					}
					break;
				case Z_M:
				case Z_P:
					if( g_pFFVConfig->RootBlockPeriodicZ == true ) {
						continue;
					}
					break;
			}
			bool bOuterBoundary = tree->checkOnOuterBoundary(node, face);
			neighborInfo[face].setOuterBoundary(bOuterBoundary);
		}
		BlockBase* block = new BlockBase(size, origin, blockSize, level, neighborInfo);
		blockManager.registerBlock(block);
	}
	blockManager.endRegisterBlock();
	blockManager.printBlockLayoutInfo();
	blockManager.printBlockLayoutInfo(g_pFFVConfig->OutputLogFilenameBlock.c_str());

	this->nBlocks = blockManager.getNumBlock();

	MPI_Allreduce(MPI_IN_PLACE, &level_min, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &level_max, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD);

	levelMin = level_min;
	levelMax = level_max;
}
g_pFFVPerfMonitor->Stop(ffv_tm_Init_RegisterBlock);
}

void FFVGrid::DistributePolygon() {
g_pFFVPerfMonitor->Start(ffv_tm_Init_DistributePolygon, 0, 0, true);
{
		if( g_pFFVConfig->PolygonGroupList.size() > 0 ) {
			if( rank == 0 ) {
				int nCells = g_pFFVConfig->LeafBlockNumberOfCells;
				Vec3r rootOrigin = g_pFFVConfig->RootBlockOrigin;
				double rootLength = g_pFFVConfig->RootBlockLength;
				Vec3i size(nCells, nCells, nCells);
				int margin = 2;
				BlockBoundingBox bbb(tree, rootOrigin, rootLength, size, margin);
				std::vector<Node*>& leafNodeArray = tree->getLeafNodeArray();
				for (int iRank = 0; iRank < nProcesses; iRank++) {
					BoundingBox box;
					for (int id=partition->getStart(iRank); id<partition->getEnd(iRank); id++) {
						Node* node = leafNodeArray[id];
						box.addBox(bbb.getBoundingBox(node));
					}
					pl->set_bounding_box(iRank, box);
				}
				pl->send_to_all();
			} else {
				pl->load_from_rank0();
			}
		}
}
g_pFFVPerfMonitor->Stop(ffv_tm_Init_DistributePolygon);

g_pFFVPerfMonitor->Start(ffv_tm_Init_RepairPolygon, 0, 0, true);
{
		cutlib::RepairPolygonData(pl);

		std::string file;
		pl->save_parallel(&file, "stl_b");

}
g_pFFVPerfMonitor->Stop(ffv_tm_Init_RepairPolygon);
}

void FFVGrid::InitVarsForCutInfo() {
	int vc                   = g_pFFVConfig->LeafBlockNumberOfVirtualCells;
	std::string updateMethod = g_pFFVConfig->TuningVCUpdate;

	int boundaryTypeNULL[NUM_FACE] = {
		1, 1, 1, 1, 1, 1,
	};
	real boundaryValueNULL[NUM_FACE] = {
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
	};
	int boundaryValueNULLINT[NUM_FACE] = {
		0, 0, 0, 0, 0, 0,
	};

	plsCut0   = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsCut1   = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsCut2   = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsCut3   = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsCut4   = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsCut5   = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsCutId0 = new LocalScalar3D<int>(blockManager,  vc, updateMethod, boundaryTypeNULL, boundaryValueNULLINT);
	plsCutId1 = new LocalScalar3D<int>(blockManager,  vc, updateMethod, boundaryTypeNULL, boundaryValueNULLINT);
	plsCutId2 = new LocalScalar3D<int>(blockManager,  vc, updateMethod, boundaryTypeNULL, boundaryValueNULLINT);
	plsCutId3 = new LocalScalar3D<int>(blockManager,  vc, updateMethod, boundaryTypeNULL, boundaryValueNULLINT);
	plsCutId4 = new LocalScalar3D<int>(blockManager,  vc, updateMethod, boundaryTypeNULL, boundaryValueNULLINT);
	plsCutId5 = new LocalScalar3D<int>(blockManager,  vc, updateMethod, boundaryTypeNULL, boundaryValueNULLINT);
	plsCut0->Fill(blockManager, 1.0);
	plsCut1->Fill(blockManager, 1.0);
	plsCut2->Fill(blockManager, 1.0);
	plsCut3->Fill(blockManager, 1.0);
	plsCut4->Fill(blockManager, 1.0);
	plsCut5->Fill(blockManager, 1.0);
	plsCutId0->Fill(blockManager, 0);
	plsCutId1->Fill(blockManager, 0);
	plsCutId2->Fill(blockManager, 0);
	plsCutId3->Fill(blockManager, 0);
	plsCutId4->Fill(blockManager, 0);
	plsCutId5->Fill(blockManager, 0);

	pNormalN        = new int   [blockManager.getNumBlock()];
	pNormalX        = new real* [blockManager.getNumBlock()];
	pNormalY        = new real* [blockManager.getNumBlock()];
	pNormalZ        = new real* [blockManager.getNumBlock()];
	plsNormalIndex0 = new LocalScalar3D<int>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULLINT);
	plsNormalIndex1 = new LocalScalar3D<int>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULLINT);
	plsNormalIndex2 = new LocalScalar3D<int>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULLINT);
	plsNormalIndex3 = new LocalScalar3D<int>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULLINT);
	plsNormalIndex4 = new LocalScalar3D<int>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULLINT);
	plsNormalIndex5 = new LocalScalar3D<int>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULLINT);

	plsPhaseId = new LocalScalar3D<int>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULLINT, 2);
	plsPhaseId->Fill(blockManager, -1);

	plsCellId = new LocalScalar3D<int>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULLINT, 2);
	plsCellId->Fill(blockManager, -1);
}

void FFVGrid::CalcCutInfo() {
	int vc                   = g_pFFVConfig->LeafBlockNumberOfVirtualCells;

#ifdef _BLOCK_IS_LARGE_
#else
#endif
	for (int n=0; n<blockManager.getNumBlock(); ++n) {
		BlockBase* block = blockManager.getBlock(n);
		Vec3i size       = block->getSize();
		Vec3r origin     = block->getOrigin();
		Vec3r blockSize  = block->getBlockSize();
		Vec3r cellSize   = block->getCellSize();

		int sz[3]              = {size.x, size.y, size.z};
		int g[1]               = {vc};
		double bpos[3]         = {origin.x, origin.y, origin.z};
		unsigned int bbsize[3] = {size.x, size.y, size.z};
		unsigned int gcsize[3] = {vc, vc, vc};
		double dx[3]           = {cellSize.x, cellSize.x, cellSize.x};
		size_t ncell[3];
		double org[3];
		for(int i=0; i<3; i++) {
			ncell[i] = bbsize[i] + 2*gcsize[i];
			org[i] = bpos[i] - gcsize[i]*dx[i];
		}

		cutlib::GridAccessor*   cutGrid   = new cutlib::Cell(org, dx);
		cutlib::CutPosArray*    cutPos    = new cutlib::CutPos32Array(ncell);
		cutlib::CutBidArray*    cutBid    = new cutlib::CutBid5Array(ncell);
		cutlib::CutNormalArray* cutNormal = new cutlib::CutNormalArray(ncell);

		cutlib::CalcCutInfo(cutGrid, pl, cutPos, cutBid, cutNormal);

		int nNormal = cutNormal->getNumNormal();
		this->pNormalN[n] = nNormal;
		this->pNormalX[n] = new real [nNormal];
		this->pNormalY[n] = new real [nNormal];
		this->pNormalZ[n] = new real [nNormal];
		cutlib::Normal* pNormal = cutNormal->getNormalDataPointer();
		for(int m=0; m<nNormal; m++) {
			this->pNormalX[n][m] = pNormal[m][0];
			this->pNormalY[n][m] = pNormal[m][1];
			this->pNormalZ[n][m] = pNormal[m][2];
		}
		cutlib::NormalIndex* nIdx = cutNormal->getNormalIndexDataPointer();
		int* pNormalIndex0 = plsNormalIndex0->GetBlockData(block);
		int* pNormalIndex1 = plsNormalIndex1->GetBlockData(block);
		int* pNormalIndex2 = plsNormalIndex2->GetBlockData(block);
		int* pNormalIndex3 = plsNormalIndex3->GetBlockData(block);
		int* pNormalIndex4 = plsNormalIndex4->GetBlockData(block);
		int* pNormalIndex5 = plsNormalIndex5->GetBlockData(block);
#pragma omp parallel for
		for(int k=vc; k<vc+size.z; k++) {
			for(int j=vc; j<vc+size.y; j++) {
				for(int i=vc; i<vc+size.x; i++) {
					int m = i + (2*vc + size.x)*(j + (2*vc + size.y)*k);
					pNormalIndex0[m] = nIdx[m][0];
					pNormalIndex1[m] = nIdx[m][1];
					pNormalIndex2[m] = nIdx[m][2];
					pNormalIndex3[m] = nIdx[m][3];
					pNormalIndex4[m] = nIdx[m][4];
					pNormalIndex5[m] = nIdx[m][5];
				}
			}
		}

		real* pCut0 = plsCut0->GetBlockData(block);
		real* pCut1 = plsCut1->GetBlockData(block);
		real* pCut2 = plsCut2->GetBlockData(block);
		real* pCut3 = plsCut3->GetBlockData(block);
		real* pCut4 = plsCut4->GetBlockData(block);
		real* pCut5 = plsCut5->GetBlockData(block);
		int* pCutId0 = plsCutId0->GetBlockData(block);
		int* pCutId1 = plsCutId1->GetBlockData(block);
		int* pCutId2 = plsCutId2->GetBlockData(block);
		int* pCutId3 = plsCutId3->GetBlockData(block);
		int* pCutId4 = plsCutId4->GetBlockData(block);
		int* pCutId5 = plsCutId5->GetBlockData(block);
#pragma omp parallel for
		for(int k=vc; k<vc+size.z; k++) {
			for(int j=vc; j<vc+size.y; j++) {
				for(int i=vc; i<vc+size.x; i++) {
					int m = i + (2*vc + size.x)*(j + (2*vc + size.y)*k);
					float cut0 = cutPos->getPos(i, j, k, 0);
					float cut1 = cutPos->getPos(i, j, k, 1);
					float cut2 = cutPos->getPos(i, j, k, 2);
					float cut3 = cutPos->getPos(i, j, k, 3);
					float cut4 = cutPos->getPos(i, j, k, 4);
					float cut5 = cutPos->getPos(i, j, k, 5);
					cutlib::BidType bid0 = cutBid->getBid(i, j, k, 0);
					cutlib::BidType bid1 = cutBid->getBid(i, j, k, 1);
					cutlib::BidType bid2 = cutBid->getBid(i, j, k, 2);
					cutlib::BidType bid3 = cutBid->getBid(i, j, k, 3);
					cutlib::BidType bid4 = cutBid->getBid(i, j, k, 4);
					cutlib::BidType bid5 = cutBid->getBid(i, j, k, 5);
					if( (int)bid0 != 0 ) {
						pCut0[m] = cut0;
					}
					if( (int)bid1 != 0 ) {
						pCut1[m] = cut1;
					}
					if( (int)bid2 != 0 ) {
						pCut2[m] = cut2;
					}
					if( (int)bid3 != 0 ) {
						pCut3[m] = cut3;
					}
					if( (int)bid4 != 0 ) {
						pCut4[m] = cut4;
					}
					if( (int)bid5 != 0 ) {
						pCut5[m] = cut5;
					}
					pCutId0[m] = bid0;
					pCutId1[m] = bid1;
					pCutId2[m] = bid2;
					pCutId3[m] = bid3;
					pCutId4[m] = bid4;
					pCutId5[m] = bid5;
				}
			}
		}
		delete cutGrid;
		delete cutPos;
		delete cutBid;
		delete cutNormal;

/*
		if( g_pFFVConfig->ShapeApproximationMethod == "cut" ) {
			real eps[1] = {g_pFFVConfig->ShapeApproximationCutoff};
			bstl_cutoff_(
							pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
							pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
							eps,
							sz, g);
		}
*/

		if( g_pFFVConfig->ShapeApproximationVoxelization ) {
			bstl_voxelize_(
							pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
							pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
							sz, g);
		}

		if( g_pFFVConfig->ShapeApproximationSymmetrization ) {
			bstl_symmetrize_(
							pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
							pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
							sz, g);
		}
	}

	plsCut0->ImposeBoundaryCondition(blockManager);
	plsCut1->ImposeBoundaryCondition(blockManager);
	plsCut2->ImposeBoundaryCondition(blockManager);
	plsCut3->ImposeBoundaryCondition(blockManager);
	plsCut4->ImposeBoundaryCondition(blockManager);
	plsCut5->ImposeBoundaryCondition(blockManager);
	plsCutId0->ImposeBoundaryCondition(blockManager);
	plsCutId1->ImposeBoundaryCondition(blockManager);
	plsCutId2->ImposeBoundaryCondition(blockManager);
	plsCutId3->ImposeBoundaryCondition(blockManager);
	plsCutId4->ImposeBoundaryCondition(blockManager);
	plsCutId5->ImposeBoundaryCondition(blockManager);
}

void FFVGrid::ClasifyCell() {
	int vc                   = g_pFFVConfig->LeafBlockNumberOfVirtualCells;

#ifdef _BLOCK_IS_LARGE_
#else
#endif
	for (int n=0; n<blockManager.getNumBlock(); ++n) {
		BlockBase* block = blockManager.getBlock(n);
		Vec3i size      = block->getSize();
		Vec3r origin    = block->getOrigin();
		Vec3r blockSize = block->getBlockSize();
		Vec3r cellSize  = block->getCellSize();

		int sz[3] = {size.x, size.y, size.z};
		int g[1]  = {vc};

		int* pCutId0 = plsCutId0->GetBlockData(block);
		int* pCutId1 = plsCutId1->GetBlockData(block);
		int* pCutId2 = plsCutId2->GetBlockData(block);
		int* pCutId3 = plsCutId3->GetBlockData(block);
		int* pCutId4 = plsCutId4->GetBlockData(block);
		int* pCutId5 = plsCutId5->GetBlockData(block);
		int* pCellId = plsCellId->GetBlockData(block);

#pragma omp parallel for
		for(int k=vc; k<vc+size.z; k++) {
			for(int j=vc; j<vc+size.y; j++) {
				for(int i=vc; i<vc+size.x; i++) {
					int m = i + (2*vc + size.x)*(j + (2*vc + size.y)*k);
					pCellId[m] = 0;
					if( (pCutId0[m] != 0) ||
							(pCutId1[m] != 0) ||
							(pCutId2[m] != 0) ||
							(pCutId3[m] != 0) ||
							(pCutId4[m] != 0) ||
							(pCutId5[m] != 0) ) {
						pCellId[m] = 1;
					}
				}
			}
		}

	}
}

void FFVGrid::FillRegion() {
}

void FFVGrid::WritePolygon(std::ofstream& ofs, float* pv) {
	ofs << "facet normal 0 0 0" << std::endl;
	ofs << "outer loop" << std::endl;
	ofs << "vertex " << pv[0] << " " << pv[1] << " " << pv[2] << std::endl;
	ofs << "vertex " << pv[3] << " " << pv[4] << " " << pv[5] << std::endl;
	ofs << "vertex " << pv[6] << " " << pv[7] << " " << pv[8] << std::endl;
	ofs << "endloop" << std::endl;
	ofs << "endfacet" << std::endl;
}

void FFVGrid::WriteCutRaw(const char* path, const char* prefix) {
	int vc                   = g_pFFVConfig->LeafBlockNumberOfVirtualCells;

	mkdir(path, 0755);

	std::ostringstream ossFileName;
	ossFileName << path;
	ossFileName << "/";
	ossFileName << prefix;
	ossFileName << "-";
	ossFileName << "cut";
	ossFileName << "-";
	ossFileName.width(5);
	ossFileName.setf(std::ios::fixed);
	ossFileName.fill('0');
	ossFileName << rank;
	ossFileName << ".stl";

	std::ofstream ofs;
	ofs.open(ossFileName.str().c_str(), std::ios::out);
	ofs.close();

#ifdef _BLOCK_IS_LARGE_
#else
#endif
	for (int n=0; n<blockManager.getNumBlock(); ++n) {
		BlockBase* block = blockManager.getBlock(n);
		Vec3i size      = block->getSize();
		Vec3r origin    = block->getOrigin();
		Vec3r blockSize = block->getBlockSize();
		Vec3r cellSize  = block->getCellSize();

		int sz[3] = {size.x, size.y, size.z};
		int g[1]  = {vc};

		double bpos[3] = {origin.x, origin.y, origin.z};
		unsigned int bbsize[3] = {size.x, size.y, size.z};
		unsigned int gcsize[3] = {vc, vc, vc};
		double dx[3] = {cellSize.x, cellSize.x, cellSize.x};
		size_t ncell[3];
		double org[3];
		for(int i=0; i<3; i++) {
			ncell[i] = bbsize[i] + 2*gcsize[i];
			org[i] = bpos[i] - gcsize[i]*dx[i];
		}

		cutlib::GridAccessor* grid   = new cutlib::Cell(org, dx);
		cutlib::CutPosArray*  cutPos = new cutlib::CutPos32Array(ncell);
		cutlib::CutBidArray*  cutBid = new cutlib::CutBid5Array(ncell);

		cutlib::CalcCutInfo(grid, pl, cutPos, cutBid);

		ofs.open(ossFileName.str().c_str(), std::ios::app);
		ofs << "solid" << std::endl;
		for(int k=vc; k<size.z+vc; k++) {
			for(int j=vc; j<size.y+vc; j++) {
				for(int i=vc; i<size.x+vc; i++) {
					int m = i + (2*vc + size.x)*(j + (2*vc + size.y)*k);
					float cut0 = cutPos->getPos(i, j, k, 0);
					float cut1 = cutPos->getPos(i, j, k, 1);
					float cut2 = cutPos->getPos(i, j, k, 2);
					float cut3 = cutPos->getPos(i, j, k, 3);
					float cut4 = cutPos->getPos(i, j, k, 4);
					float cut5 = cutPos->getPos(i, j, k, 5);
					cutlib::BidType bidp0 = cutBid->getBid(i, j, k, 0);
					cutlib::BidType bidp1 = cutBid->getBid(i, j, k, 1);
					cutlib::BidType bidp2 = cutBid->getBid(i, j, k, 2);
					cutlib::BidType bidp3 = cutBid->getBid(i, j, k, 3);
					cutlib::BidType bidp4 = cutBid->getBid(i, j, k, 4);
					cutlib::BidType bidp5 = cutBid->getBid(i, j, k, 5);

					float x0 = origin[0] + (i - vc)*dx[0];
					float y0 = origin[1] + (j - vc)*dx[1];
					float z0 = origin[2] + (k - vc)*dx[2];
					float x1 = origin[0] + (i + 1 - vc)*dx[0];
					float y1 = origin[1] + (j + 1 - vc)*dx[1];
					float z1 = origin[2] + (k + 1 - vc)*dx[2];

					float v[9];
					float x[2] = {x0, x1};
					float y[2] = {y0, y1};
					float z[2] = {z0, z1};

					if( bidp0 != 0 ) {
					}
					if( bidp1 != 0 ) {
						v[0] = x[1]; v[1] = y[0]; v[2] = z[0];
						v[3] = x[1]; v[4] = y[0]; v[5] = z[1];
						v[6] = x[1]; v[7] = y[1]; v[8] = z[0];
						WritePolygon(ofs, v);

						v[0] = x[1]; v[1] = y[1]; v[2] = z[1];
						v[6] = x[1]; v[7] = y[0]; v[8] = z[1];
						v[3] = x[1]; v[4] = y[1]; v[5] = z[0];
						WritePolygon(ofs, v);
					}
					if( bidp2 != 0 ) {
					}
					if( bidp3 != 0 ) {
						v[0] = x[0]; v[1] = y[1]; v[2] = z[0];
						v[6] = x[0]; v[7] = y[1]; v[8] = z[1];
						v[3] = x[1]; v[4] = y[1]; v[5] = z[0];
						WritePolygon(ofs, v);

						v[0] = x[1]; v[1] = y[1]; v[2] = z[1];
						v[6] = x[1]; v[7] = y[1]; v[8] = z[0];
						v[3] = x[0]; v[4] = y[1]; v[5] = z[1];
						WritePolygon(ofs, v);
					}
					if( bidp4 != 0 ) {
					}
					if( bidp5 != 0 ) {
						v[0] = x[0]; v[1] = y[0]; v[2] = z[1];
						v[6] = x[0]; v[7] = y[1]; v[8] = z[1];
						v[3] = x[1]; v[4] = y[0]; v[5] = z[1];
						WritePolygon(ofs, v);

						v[0] = x[1]; v[1] = y[1]; v[2] = z[1];
						v[3] = x[0]; v[4] = y[1]; v[5] = z[1];
						v[6] = x[1]; v[7] = y[0]; v[8] = z[1];
						WritePolygon(ofs, v);
					}
				}
			}
		}
		ofs << "endsolid" << std::endl;
		ofs.close();

		delete grid;
		delete cutPos;
		delete cutBid;
	}
}

