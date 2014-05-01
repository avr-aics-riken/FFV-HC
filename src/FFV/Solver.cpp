#include "Solver.h"

#include <iostream>
#include <mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "real.h"
#include "comm.h"

#include "SimpleDivider.h"
#include "PolygonBBoxDivider.h"
#include "SphereDivider2.h"
#include "SphereDivider3.h"

#include "BlockBoundingBox.h"

#include "Polylib.h"
#include "Cutlib.h"
#include "CutInfo/CutInfo.h"
#include "CutInfo/CutNormalArray.h"
#include "GridAccessor/Cell.h"

#include "bcut.h"
#include "bstl.h"
#include "bils.h"
#include "FFVPM.h"

#define GLOBAL_VALUE_DEFINE
#include "FFVGlobalVars.h"
#undef GLOBAL_VALUE_DEFINE

Solver::Solver()
	: blockManager(BlockManager::getInstance()) {
}

Solver::~Solver() {
	MPI::Finalize();
}

int Solver::Init(int argc, char** argv){
/* ---------------------------------------------------------- */
/* Init MPI                                                   */
/* ---------------------------------------------------------- */
/* ---------------------------------------------------------- */
	MPI::Init(argc, argv);
	MPI::Comm& comm = MPI::COMM_WORLD;
	this->myrank = comm.Get_rank();

srand( this->myrank + 1 );

	if( argc != 2 ) {
		PrintLog(0, "usage: %s configfile", argv[0]);
		comm.Abort(EX_USAGE);
	}
/* ---------------------------------------------------------- */
	PrintLog(1, FFV_SOLVERNAME);
	PrintLog(2, "%-20s : %s", "Version", FFV_VERSION);
	PrintLog(2, "%-20s : %s", "Revision", FFV_REVISION);
	PrintLog(2, "%-20s : %s", "Build date", BUILD_DATE);
	PrintLog(2, "%-20s : %s", "Configuration file", argv[1]);
	PrintLog(2, "%-20s : %d", "MPI processes", comm.Get_size());
#ifdef _OPENMP
	PrintLog(2, "%-20s : %d", "OpenMP threads", omp_get_max_threads());
#endif
/* ---------------------------------------------------------- */

/* ---------------------------------------------------------- */
/* Init Config                                                */
/* ---------------------------------------------------------- */
	PrintLog(1, "Loading configuration file");
/* ---------------------------------------------------------- */
	g_pFFVConfig = new FFVConfig();
	g_pFFVConfig->Load(argv[1]);
/* ---------------------------------------------------------- */
	PrintLog(2, "%-20s : %s", "tree type", g_pFFVConfig->TreeType.c_str());
	PrintLog(2, "%-20s : %s", "ordering",  g_pFFVConfig->TuningBlockOrdering.c_str());
	PrintLog(2, "%-20s : %d", "block size", g_pFFVConfig->LeafBlockNumberOfCells);
	PrintLog(2, "%-20s : %d", "vc width", g_pFFVConfig->LeafBlockNumberOfVirtualCells);
	PrintLog(2, "Completed");
/* ---------------------------------------------------------- */

/* ---------------------------------------------------------- */
/* Init PMlib                                                 */
/* ---------------------------------------------------------- */
	PrintLog(1, "Initializing PMlib");
/* ---------------------------------------------------------- */
	g_pPM = new pm_lib::PerfMonitor;
	g_pPM->initialize(tm_END);
	g_pPM->setRankInfo(this->myrank);
	g_pPM->setProperties(tm_Init_LoadSTL,      "LoadSTL",      pm_lib::PerfMonitor::COMM, true);
	g_pPM->setProperties(tm_Init_DivideDomain, "DivideDomain", pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_Init_CreateTree,   "CreateTree",   pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_Init_DistributeSTL,"DistributeSTL",pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_Init_RepairPolygonData, "RepairPolygonData",   pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_Init_CalcCutInfo,  "CalcCutInfo",  pm_lib::PerfMonitor::CALC, false);
	g_pPM->setProperties(tm_Init_CalcCutInfo01,"CalcCutInfo01",pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_Init_CalcCutInfo02,"CalcCutInfo02",pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_Init_CalcCutInfo03,"CalcCutInfo03",pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_Init_CalcCutInfo04,"CalcCutInfo04",pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_Init_CalcCutInfo05,"CalcCutInfo05",pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_Init_CalcCutInfo06,"CalcCutInfo06",pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_Init_Filling,      "Filling",      pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_Init_GeometricalProperties, "GeometricalProperties",  pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_Init_InitVars,     "InitVars",     pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_Update,     "Update",    pm_lib::PerfMonitor::CALC, false);
	g_pPM->setProperties(tm_UpdateT,    "UpdateT",   pm_lib::PerfMonitor::CALC, false);
	g_pPM->setProperties(tm_UpdateUX,   "UpdateUX",  pm_lib::PerfMonitor::CALC, false);
	g_pPM->setProperties(tm_UpdateUY,   "UpdateUY",  pm_lib::PerfMonitor::CALC, false);
	g_pPM->setProperties(tm_UpdateUZ,   "UpdateUZ",  pm_lib::PerfMonitor::CALC, false);
	g_pPM->setProperties(tm_UpdateP,    "UpdateP",   pm_lib::PerfMonitor::CALC, false);
	g_pPM->setProperties(tm_UpdateU,    "UpdateU",   pm_lib::PerfMonitor::CALC, false);
	g_pPM->setProperties(tm_UpdateT01,  "UpdateT01", pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_UpdateT02,  "UpdateT02", pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_UpdateT03,  "UpdateT03", pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_UpdateT04,  "UpdateT04", pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_UpdateT05,  "UpdateT05", pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_UpdateUX01, "UpdateUX01",pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_UpdateUX02, "UpdateUX02",pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_UpdateUX03, "UpdateUX03",pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_UpdateUX04, "UpdateUX04",pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_UpdateUX05, "UpdateUX05",pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_UpdateUY01, "UpdateUY01",pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_UpdateUY02, "UpdateUY02",pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_UpdateUY03, "UpdateUY03",pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_UpdateUY04, "UpdateUY04",pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_UpdateUY05, "UpdateUY05",pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_UpdateUZ01, "UpdateUZ01",pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_UpdateUZ02, "UpdateUZ02",pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_UpdateUZ03, "UpdateUZ03",pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_UpdateUZ04, "UpdateUZ04",pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_UpdateUZ05, "UpdateUZ05",pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_UpdateP01,  "UpdateP01", pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_UpdateP02,  "UpdateP02", pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_UpdateP03,  "UpdateP03", pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_UpdateP04,  "UpdateP04", pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_UpdateP05,  "UpdateP05", pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_UpdateP06,  "UpdateP06", pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_UpdateP07,  "UpdateP07", pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_UpdateU01,  "UpdateU01", pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_UpdateU02,  "UpdateU02", pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_UpdateU03,  "UpdateU03", pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_Print,      "Print",     pm_lib::PerfMonitor::CALC, false);
	g_pPM->setProperties(tm_PrintTime,  "PrintTime", pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_PrintILS,   "PrintILS",  pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_PrintStats, "PrintStats",pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_PrintForce, "PrintForce",pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_PrintHeatFlux, "PrintHeatFlux",pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_PrintData,  "PrintData", pm_lib::PerfMonitor::CALC, true);
	int nThreads = 1;
#ifdef _OPENMP
	nThreads = omp_get_max_threads();
#endif
	g_pPM->setParallelMode("Hybrid", nThreads, comm.Get_size());
/* ---------------------------------------------------------- */
	PrintLog(2, "Completed");
/* ---------------------------------------------------------- */

/* ---------------------------------------------------------- */
/* Load STL                                                   */
/* ---------------------------------------------------------- */
	PrintLog(1, "Loading STL file(s)");
/* ---------------------------------------------------------- */
	PM_Start(tm_Init_LoadSTL, 0, 0, true);
	PolylibNS::BCMPolylib* pl = new PolylibNS::BCMPolylib;
	struct stat st;
	int ret = stat(g_pFFVConfig->PolylibConfig.c_str(), &st);
	if( ret == 0 ) {
		pl->load(g_pFFVConfig->PolylibConfig);
	} else {
	}
	PM_Stop(tm_Init_LoadSTL);
/* ---------------------------------------------------------- */
	PrintLog(2, "Completed");
/* ---------------------------------------------------------- */

/* ---------------------------------------------------------- */
/* Divide domain                                              */
/* ---------------------------------------------------------- */
	PrintLog(1, "Dividing domain");
/* ---------------------------------------------------------- */
	PM_Start(tm_Init_DivideDomain, 0, 0, true);
	Divider* divider = 0;
	if(myrank == 0) {
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
		} else {
			exit(EX_READ_CONFIG);
		}
	}
	PM_Stop(tm_Init_DivideDomain);
/* ---------------------------------------------------------- */
	PrintLog(2, "Completed");
/* ---------------------------------------------------------- */

/* ---------------------------------------------------------- */
/* Create tree                                                */
/* ---------------------------------------------------------- */
	PrintLog(1, "Creating tree");
/* ---------------------------------------------------------- */
	PM_Start(tm_Init_CreateTree, 0, 0, true);
	if(myrank == 0) {
		BCMOctree::Ordering ordering;
		if (g_pFFVConfig->TuningBlockOrdering == "Z") {
			ordering = BCMOctree::Z;
		} else if (g_pFFVConfig->TuningBlockOrdering == "Hilbert") {
			ordering = BCMOctree::HILBERT;
		} else if (g_pFFVConfig->TuningBlockOrdering == "random") {
			ordering = BCMOctree::RANDOM;
		} else if (g_pFFVConfig->TuningBlockOrdering == "PedigreeList") {
			ordering = BCMOctree::PEDIGREELIST;
		} else {
			exit(EX_READ_CONFIG);
		}
		tree = new BCMOctree(rootGrid, divider, ordering);
    tree->broadcast();
	} else {
    tree = BCMOctree::ReceiveFromMaster();
	}
	PM_Stop(tm_Init_CreateTree);
/* ---------------------------------------------------------- */
  int numLeafNode = tree->getNumLeafNode();
	partition = new Partition(comm.Get_size(), numLeafNode);
	for(int n=0; n<comm.Get_size(); n++) {
		if( n==0 ) {
			PrintLog(2, "%-20s : #%05d [%04d:%04d] (%04d)",  "Partitions", n, partition->getStart(n), partition->getEnd(n)-1, partition->getEnd(n) - partition->getStart(n));
		} else {
			PrintLog(2, "%-20s : #%05d [%04d:%04d] (%04d)",  "", n, partition->getStart(n), partition->getEnd(n)-1, partition->getEnd(n) - partition->getStart(n));
		}
	}
	PrintLog(2, "Completed");
/* ---------------------------------------------------------- */

/* ---------------------------------------------------------- */
/* Register blocks                                            */
/* ---------------------------------------------------------- */
	PrintLog(1, "Registering block(s)");
/* ---------------------------------------------------------- */
  // ブロック内のセル数
	Vec3i size(g_pFFVConfig->LeafBlockNumberOfCells, g_pFFVConfig->LeafBlockNumberOfCells, g_pFFVConfig->LeafBlockNumberOfCells);
	Vec3r rootOrigin = g_pFFVConfig->RootBlockOrigin;
	double rootLength = g_pFFVConfig->RootBlockLength;
  std::vector<Node*>& leafNodeArray = tree->getLeafNodeArray();
  for (int id = partition->getStart(myrank); id < partition->getEnd(myrank); id++) {
    Node* node = leafNodeArray[id];
		int level = node->getLevel();
		Vec3r origin = tree->getOrigin(node) * rootLength + rootOrigin;
		Vec3r blockSize = node->getBlockSize() * rootLength;
		NeighborInfo* neighborInfo = tree->makeNeighborInfo(node, partition);

		for (int i = 0; i < NUM_FACE; ++i) {
			Face face = Face(i);
			switch(i) {
				case 0:
				case 1:
					if( g_pFFVConfig->RootBlockPeriodicX == true ) {
						continue;
					}
					break;
				case 2:
				case 3:
					if( g_pFFVConfig->RootBlockPeriodicY == true ) {
						continue;
					}
					break;
				case 4:
				case 5:
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
/* ---------------------------------------------------------- */
//  blockManager.printBlockLayoutInfo();
  blockManager.printBlockLayoutInfo("data-block.txt");
	PrintLog(2, "Completed");
/* ---------------------------------------------------------- */

/* ---------------------------------------------------------- */
/* Distribute STL                                             */
/* ---------------------------------------------------------- */
	PrintLog(1, "Distributing polygon(s)");
/* ---------------------------------------------------------- */
	PM_Start(tm_Init_DistributeSTL, 0, 0, true);
	if( myrank == 0 ) {
		Vec3r rootOrigin = g_pFFVConfig->RootBlockOrigin;
		double rootLength = g_pFFVConfig->RootBlockLength;
		int margin = g_pFFVConfig->LeafBlockNumberOfVirtualCells;
		BlockBoundingBox bbb(tree, rootOrigin, rootLength, size, margin);

		for (int iRank = 0; iRank < comm.Get_size(); iRank++) {
			BoundingBox box;
			for (int id = partition->getStart(iRank); id < partition->getEnd(iRank); id++) {
				Node* node = leafNodeArray[id];
				box.addBox(bbb.getBoundingBox(node));
			}
			pl->set_bounding_box(iRank, box);
		}
		pl->send_to_all();
	} else {
		pl->load_from_rank0();
	}
	PM_Stop(tm_Init_DistributeSTL);
/* ---------------------------------------------------------- */
	PrintLog(2, "Completed");
/* ---------------------------------------------------------- */

/* ---------------------------------------------------------- */
/* Repair polygon                                             */
/* ---------------------------------------------------------- */
	PrintLog(1, "Repairing polygon(s)");
/* ---------------------------------------------------------- */
	PM_Start(tm_Init_RepairPolygonData, 0, 0, true);
	cutlib::RepairPolygonData(pl);
	PM_Stop(tm_Init_RepairPolygonData);
/* ---------------------------------------------------------- */
	PrintLog(2, "Completed");
/* ---------------------------------------------------------- */

/* ---------------------------------------------------------- */
/* Misc.                                                      */
/* ---------------------------------------------------------- */
	PrintLog(1, "Summarizing block layout");
/* ---------------------------------------------------------- */
	maxLevel = 0;
	minLevel = INT_MAX;
	for (int n=0; n<blockManager.getNumBlock(); ++n) {
		BlockBase* block = blockManager.getBlock(n);
		int level = block->getLevel();
		if( level >= maxLevel ) {
			maxLevel = level;
		}
		if( level <= minLevel ) {
			minLevel = level;
		}
	}
	diffLevel = maxLevel - minLevel;
	vc = g_pFFVConfig->LeafBlockNumberOfVirtualCells;
	updateMethod = g_pFFVConfig->TuningVCUpdate;
/* ---------------------------------------------------------- */
	PrintLog(2, "%-20s : %d", "num", blockManager.getNumBlock());
	PrintLog(2, "%-20s : %d", "min level", minLevel);
	PrintLog(2, "%-20s : %d", "max level", maxLevel);
	PrintLog(2, "Completed");
/* ---------------------------------------------------------- */


/* ---------------------------------------------------------- */
/* Init physical parameters                                   */
/* ---------------------------------------------------------- */
	dt		= g_pFFVConfig->TimeControlTimeStepDeltaT;

	rhof	= g_pFFVConfig->MediumTableFluid[0].rho;
	cpf		= g_pFFVConfig->MediumTableFluid[0].cp;
	kf		= g_pFFVConfig->MediumTableFluid[0].k;
	mu		= g_pFFVConfig->MediumTableFluid[0].mu;

	rhos	= g_pFFVConfig->MediumTableFluid[0].rho;
	cps		= g_pFFVConfig->MediumTableFluid[0].cp;
	ks		= g_pFFVConfig->MediumTableFluid[0].k;
/* ---------------------------------------------------------- */


	int boundaryTypeNULL[NUM_FACE] = {
		1, 1, 1, 1, 1, 1,
	};

	real boundaryValueNULL[NUM_FACE] = {
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
	};

	int boundaryValueNULLINT[NUM_FACE] = {
		0, 0, 0, 0, 0, 0,
	};


/* ---------------------------------------------------------- */
/* Calc cutinfo                                               */
/* ---------------------------------------------------------- */
	PrintLog(1, "Computing cuts");
/* ---------------------------------------------------------- */
	PM_Start(tm_Init_CalcCutInfo, 0, 0, true);

	plsCut0 = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsCut1 = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsCut2 = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsCut3 = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsCut4 = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsCut5 = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsCutId0 = new LocalScalar3D<int>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULLINT);
	plsCutId1 = new LocalScalar3D<int>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULLINT);
	plsCutId2 = new LocalScalar3D<int>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULLINT);
	plsCutId3 = new LocalScalar3D<int>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULLINT);
	plsCutId4 = new LocalScalar3D<int>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULLINT);
	plsCutId5 = new LocalScalar3D<int>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULLINT);
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

	plsPhaseId = new LocalScalar3D<int>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULLINT, 2);
	plsPhaseId->Fill(blockManager, -1);

	pNormalN = new int   [blockManager.getNumBlock()];
	pNormalX = new real* [blockManager.getNumBlock()];
	pNormalY = new real* [blockManager.getNumBlock()];
	pNormalZ = new real* [blockManager.getNumBlock()];
	plsNormalIndex0 = new LocalScalar3D<int>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULLINT);
	plsNormalIndex1 = new LocalScalar3D<int>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULLINT);
	plsNormalIndex2 = new LocalScalar3D<int>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULLINT);
	plsNormalIndex3 = new LocalScalar3D<int>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULLINT);
	plsNormalIndex4 = new LocalScalar3D<int>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULLINT);
	plsNormalIndex5 = new LocalScalar3D<int>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULLINT);

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
		int g[1] = {vc};

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

		cutlib::GridAccessor*   grid   = new cutlib::Cell(org, dx);
		cutlib::CutPosArray*    cutPos = new cutlib::CutPos32Array(ncell);
		cutlib::CutBidArray*    cutBid = new cutlib::CutBid5Array(ncell);
		cutlib::CutNormalArray* cutNormal = new cutlib::CutNormalArray(ncell);

//		CutInfoCell(org, dx, pl, cutPos, cutBid);
PM_Start(tm_Init_CalcCutInfo01, 0, 0, false);
		int ret = CalcCutInfo(grid, pl, cutPos, cutBid, cutNormal);
PM_Stop(tm_Init_CalcCutInfo01);

		pNormalN[n] = cutNormal->getNumNormal();

		std::cout << "CalcCutInfo: " << ret << " " << pNormalN[n] << std::endl;

		pNormalX[n] = new real [pNormalN[n]];
		pNormalY[n] = new real [pNormalN[n]];
		pNormalZ[n] = new real [pNormalN[n]];
		cutlib::Normal* pNormal = cutNormal->getNormalDataPointer();
		for(int m=0; m<pNormalN[n]; m++) {
			pNormalX[n][m] = pNormal[m][0];
			pNormalY[n][m] = pNormal[m][1];
			pNormalZ[n][m] = pNormal[m][2];
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

PM_Start(tm_Init_CalcCutInfo02, 0, 0, false);
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
		delete grid;
		delete cutPos;
		delete cutBid;
		delete cutNormal;
PM_Stop(tm_Init_CalcCutInfo02);

PM_Start(tm_Init_CalcCutInfo03, 0, 0, false);
		if( g_pFFVConfig->ShapeApproximationMethod == "cut" ) {
			real eps[1] = {g_pFFVConfig->ShapeApproximationCutoff};
			bstl_cutoff_(
							pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
							pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
							eps,
							sz, g);
		}
PM_Stop(tm_Init_CalcCutInfo03);

PM_Start(tm_Init_CalcCutInfo04, 0, 0, false);
		if( g_pFFVConfig->ShapeApproximationVoxelization ) {
			bstl_voxelize_(
							pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
							pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
							sz, g);
		}
PM_Stop(tm_Init_CalcCutInfo04);

PM_Start(tm_Init_CalcCutInfo05, 0, 0, false);
		if( g_pFFVConfig->ShapeApproximationSymmetrization ) {
			bstl_symmetrize_(
							pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
							pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
							sz, g);
		}
PM_Stop(tm_Init_CalcCutInfo05);
	}

PM_Start(tm_Init_CalcCutInfo06, 0, 0, true);
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
PM_Stop(tm_Init_CalcCutInfo06);

	PM_Stop(tm_Init_CalcCutInfo);
/* ---------------------------------------------------------- */
	if( g_pFFVConfig->GridGenerationOutputSTL ) {
		PrintLog(1, "Printing STL files for cut info");
		PrintCut(0);
		PrintHole(0);
		WriteGrid(
					rootGrid,
					tree,
					partition);
		MPI_Barrier(MPI_COMM_WORLD);
	}
	PrintLog(2, "Completed");
/* ---------------------------------------------------------- */


	{
/* ---------------------------------------------------------- */
/* Detect zero-cut                                            */
/* ---------------------------------------------------------- */
		PrintLog(1, "Detecting zero-cut");
/* ---------------------------------------------------------- */
		int countLocal = 0;
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for reduction(+: countLocal)
#endif
		for (int n=0; n<blockManager.getNumBlock(); ++n) {
			BlockBase* block = blockManager.getBlock(n);
			Vec3i size = block->getSize();
			Vec3r origin = block->getOrigin();
			Vec3r blockSize = block->getBlockSize();
			Vec3r cellSize = block->getCellSize();

			int sz[3] = {size.x, size.y, size.z};
			int g[1] = {vc};

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

			int count = 0;
			bstl_detect_zerocut_(
							pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
							pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
							&count,
							sz, g);

			countLocal += count;
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

		int countTmp = countLocal;
		MPI_Allreduce(&countTmp, &countLocal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

/* ---------------------------------------------------------- */
		PrintLog(2, "%-20s : %d", "Zero distance", countLocal);
		PrintLog(2, "Completed");
/* ---------------------------------------------------------- */
	}

	{
/* ---------------------------------------------------------- */
/* Fill holes                                                 */
/* ---------------------------------------------------------- */
		PrintLog(1, "Filling hole(s)");
/* ---------------------------------------------------------- */
		int countLocal = 0;
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for reduction(+: countLocal)
#endif
		for (int n=0; n<blockManager.getNumBlock(); ++n) {
			BlockBase* block = blockManager.getBlock(n);
			Vec3i size = block->getSize();
			Vec3r origin = block->getOrigin();
			Vec3r blockSize = block->getBlockSize();
			Vec3r cellSize = block->getCellSize();

			int sz[3] = {size.x, size.y, size.z};
			int g[1] = {vc};

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

			int count = 0;
			int bClose = 1;
			if( g_pFFVConfig->GridGenerationHoleFilling == false ) {
				bClose = 0;
			}
			bstl_fill_holes_(
							pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
							pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
							&count,
							&bClose,
							sz, g);

			countLocal += count;
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

		int countTmp = countLocal;
		MPI_Allreduce(&countTmp, &countLocal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

/* ---------------------------------------------------------- */
		PrintLog(2, "%-20s : %d", "Hole faces(single)", countLocal);
		PrintLog(2, "Completed");
/* ---------------------------------------------------------- */
	}

	{
/* ---------------------------------------------------------- */
/* Fill holes v2                                              */
/* ---------------------------------------------------------- */
		PrintLog(1, "Filling hole(s) 2");
/* ---------------------------------------------------------- */
		int nIterationCount = 0;
		int countLocal = 0;
		int countTotal = 0;
		do {
			countLocal = 0;

#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for reduction(+: countLocal)
#endif
			for (int n=0; n<blockManager.getNumBlock(); ++n) {
				BlockBase* block = blockManager.getBlock(n);
				Vec3i size = block->getSize();
				Vec3r origin = block->getOrigin();
				Vec3r blockSize = block->getBlockSize();
				Vec3r cellSize = block->getCellSize();

				int sz[3] = {size.x, size.y, size.z};
				int g[1] = {vc};

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

				int count = 0;
				int bClose = 1;
				if( g_pFFVConfig->GridGenerationHoleFilling2 == false ) {
					bClose = 0;
				}

				bstl_fill_holes_v2_(
								pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
								pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
								&count,
								&bClose,
								sz, g);

				countLocal += count;
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

			int countTmp = countLocal;
			MPI_Allreduce(&countTmp, &countLocal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

			nIterationCount++;
			countTotal += countLocal;
		}while( countLocal > 0 && g_pFFVConfig->GridGenerationHoleFilling2 == true);
/* ---------------------------------------------------------- */
		PrintLog(2, "%-20s : %d (#Ite.: %d)", "Hole faces(multi)", countTotal, nIterationCount);
		PrintLog(2, "Completed");
/* ---------------------------------------------------------- */
	}

/* ---------------------------------------------------------- */
/* Filling                                                    */
/* ---------------------------------------------------------- */
	PrintLog(1, "Filling fluid");
/* ---------------------------------------------------------- */
	PM_Start(tm_Init_Filling, 0, 0, true);

	real xs = g_pFFVConfig->FillingOrigin.x;
	real ys = g_pFFVConfig->FillingOrigin.y;
	real zs = g_pFFVConfig->FillingOrigin.z;
	PrintLog(2, "%-20s : %f %f %f", "Seed for FLUID", xs, ys, zs);

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

		int* pPhaseId = plsPhaseId->GetBlockData(block);

		bcut_set_fluidseed_(
			pPhaseId,
			&xs, &ys, &zs,
			&dx,
			org,
			sz, g);

	}
	plsPhaseId->ImposeBoundaryCondition(blockManager);

	{
		int nIterationCount = 0;
		long int nCellsChanged = 0;
		do {
			nCellsChanged = 0;
#ifdef _BLOCK_IS_LARGE_
#else
//#pragma omp parallel for reduction(+: nCellsChanged)
#endif
			for (int n=0; n<blockManager.getNumBlock(); ++n) {
				BlockBase* block = blockManager.getBlock(n);
				Vec3i size = block->getSize();
				Vec3r origin = block->getOrigin();
				Vec3r blockSize = block->getBlockSize();
				Vec3r cellSize = block->getCellSize();

				int sz[3] = {size.x, size.y, size.z};
				int g[1] = {vc};
				int nc[3] = {size.x + 2*vc, size.y + 2*vc, size.z + 2*vc};

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

				int* pPhaseId = plsPhaseId->GetBlockData(block);
#ifdef _BLOCK_IS_LARGE_
//#pragma omp parallel for reduction(+: nCellsChanged)
#else
#endif
				for(int k=vc; k<=size.z+vc-1; k++) {
					for(int j=vc; j<=size.y+vc-1; j++) {
						for(int i=vc; i<=size.x+vc-1; i++) {
							int mp = i + nc[0]*( j + nc[1]*k );
							int mw = i-1 + nc[0]*( j + nc[1]*k );
							int me = i+1 + nc[0]*( j + nc[1]*k );
							int ms = i + nc[0]*( j-1 + nc[1]*k );
							int mn = i + nc[0]*( j+1 + nc[1]*k );
							int mb = i + nc[0]*( j + nc[1]*(k-1) );
							int mt = i + nc[0]*( j + nc[1]*(k+1) );

							int cidp0 = pCutId0[mp];
							int cidp1 = pCutId1[mp];
							int cidp2 = pCutId2[mp];
							int cidp3 = pCutId3[mp];
							int cidp4 = pCutId4[mp];
							int cidp5 = pCutId5[mp];

							if( pPhaseId[mp] > 0 ) {
								continue;
							}

							if( (pPhaseId[mw] == 1 && cidp0 == 0) ||
									(pPhaseId[me] == 1 && cidp1 == 0) ||
									(pPhaseId[ms] == 1 && cidp2 == 0) ||
									(pPhaseId[mn] == 1 && cidp3 == 0) ||
									(pPhaseId[mb] == 1 && cidp4 == 0) ||
									(pPhaseId[mt] == 1 && cidp5 == 0) ) {
								pPhaseId[mp] = 1;
								nCellsChanged++;
							}
						}
					}
				}
			}
			plsPhaseId->ImposeBoundaryCondition(blockManager);

			long int nCellsChangedTmp = nCellsChanged;
			MPI_Allreduce(&nCellsChangedTmp, &nCellsChanged, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);

			nIterationCount++;
		}while(nCellsChanged>0);

		long int count = 0;
		long int countS = 0;
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
			int nc[3] = {size.x + 2*vc, size.y + 2*vc, size.z + 2*vc};

			int* pPhaseId = plsPhaseId->GetBlockData(block);
			for(int k=vc; k<=size.z+vc-1; k++) {
				for(int j=vc; j<=size.y+vc-1; j++) {
					for(int i=vc; i<=size.x+vc-1; i++) {
						int mp = i + nc[0]*( j + nc[1]*k );
						if( pPhaseId[mp] > 0 ) {
							count++;
						} else {
							countS++;
						}
					}
				}
			}
		}

		long int countTmp = count;
		MPI_Allreduce(&countTmp, &count, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);

		countTmp = countS;
		MPI_Allreduce(&countTmp, &countS, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);

/* ---------------------------------------------------------- */
		PrintLog(2, "%-20s : %d", "Iteration", nIterationCount);
		PrintLog(2, "%-20s : %d", "FLUID cells", count);
		PrintLog(2, "%-20s : %d", "SOLID cells", countS);
/* ---------------------------------------------------------- */
	}
	PM_Stop(tm_Init_Filling);
/* ---------------------------------------------------------- */
	if( g_pFFVConfig->GridGenerationOutputSTL ) {
		PrintLog(1, "Printing STL files for cut info");
		PrintCut(1);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	PrintLog(2, "Completed");
/* ---------------------------------------------------------- */


/* ---------------------------------------------------------- */
/* Compute geometrical properties                             */
/* ---------------------------------------------------------- */
	PrintLog(1, "Computing geometrical properties of flow");
/* ---------------------------------------------------------- */
	double VGlobal = 0.0;
	double SGlobal = 0.0;
	double SGlobal2 = 0.0;
	PM_Start(tm_Init_GeometricalProperties, 0, 0, true);
	{
			double VLocal = 0.0;
			double SLocal = 0.0;
			double SLocal2 = 0.0;
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
			int g[1] = {vc};

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

			int* pPhaseId = plsPhaseId->GetBlockData(block);

			int* pNormalIndex0 = plsNormalIndex0->GetBlockData(block);
			int* pNormalIndex1 = plsNormalIndex1->GetBlockData(block);
			int* pNormalIndex2 = plsNormalIndex2->GetBlockData(block);
			int* pNormalIndex3 = plsNormalIndex3->GetBlockData(block);
			int* pNormalIndex4 = plsNormalIndex4->GetBlockData(block);
			int* pNormalIndex5 = plsNormalIndex5->GetBlockData(block);

#pragma omp parallel for reduction(+:VLocal,SLocal,SLocal2)
			for(int k=vc; k<vc+size.z; k++) {
				for(int j=vc; j<vc+size.y; j++) {
					for(int i=vc; i<vc+size.x; i++) {
						int m = i + (2*vc + size.x)*(j + (2*vc + size.y)*k);

						if( pPhaseId[m] == 1 ) {
							VLocal += dx[0]*dx[1]*dx[2];
						} else {
							continue;
						}

						if( pCutId0[m] != 0 ) {
							int nIdx = pNormalIndex0[m];
							double nx = pNormalX[n][nIdx];
							SLocal += fabs(nx)*dx[1]*dx[2];
							SLocal2 += dx[1]*dx[2];
						}
						if( pCutId1[m] != 0 ) {
							int nIdx = pNormalIndex1[m];
							double nx = pNormalX[n][nIdx];
							SLocal += fabs(nx)*dx[1]*dx[2];
							SLocal2 += dx[1]*dx[2];
						}
						if( pCutId2[m] != 0 ) {
							int nIdx = pNormalIndex2[m];
							double ny = pNormalY[n][nIdx];
							SLocal += fabs(ny)*dx[2]*dx[0];
							SLocal2 += dx[2]*dx[0];
						}
						if( pCutId3[m] != 0 ) {
							int nIdx = pNormalIndex3[m];
							double ny = pNormalY[n][nIdx];
							SLocal += fabs(ny)*dx[2]*dx[0];
							SLocal2 += dx[2]*dx[0];
						}
						if( pCutId4[m] != 0 ) {
							int nIdx = pNormalIndex4[m];
							double nz = pNormalZ[n][nIdx];
							SLocal += fabs(nz)*dx[0]*dx[1];
							SLocal2 += dx[0]*dx[1];
						}
						if( pCutId5[m] != 0 ) {
							int nIdx = pNormalIndex5[m];
							double nz = pNormalZ[n][nIdx];
							SLocal += fabs(nz)*dx[0]*dx[1];
							SLocal2 += dx[0]*dx[1];
						}
					}
				}
			}
		}
		MPI_Allreduce(&VLocal, &VGlobal, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&SLocal, &SGlobal, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&SLocal2, &SGlobal2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
	}
/* ---------------------------------------------------------- */
	PrintLog(2, "%-20s : %f", "Volume",       VGlobal);
	PrintLog(2, "%-20s : %f", "Surface area", SGlobal);
	PrintLog(2, "%-20s : %f", "Surface area (voxel)", SGlobal2);
/* ---------------------------------------------------------- */
	PM_Stop(tm_Init_GeometricalProperties, 0, 0, true);
/* ---------------------------------------------------------- */
	MPI_Barrier(MPI_COMM_WORLD);
	PrintLog(2, "Completed");
/* ---------------------------------------------------------- */

	if( g_pFFVConfig->OperationMode == "gridgeneration" ) {
		return EX_FAILURE;
		return EX_SUCCESS;
	}

/* ---------------------------------------------------------- */
/* Init vars                                                  */
/* ---------------------------------------------------------- */
		PrintLog(1, "Initializing variables");
/* ---------------------------------------------------------- */

	PM_Start(tm_Init_InitVars, 0, 0, true);

/* ---------------------------------------------------------- */
/* Init mask                                                  */
/* ---------------------------------------------------------- */
	plsMaskId = new LocalScalar3D<int>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULLINT);
	plsMaskId->Fill(blockManager, 0);
	plsM = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsM->Fill(blockManager, 0.0);
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
#endif
	for (int n=0; n<blockManager.getNumBlock(); ++n) {
		BlockBase* block = blockManager.getBlock(n);
		Vec3i size = block->getSize();
		Vec3r origin = block->getOrigin();
		Vec3r blockSize = block->getBlockSize();
		Vec3r cellSize = block->getCellSize();

		int sz[3] = {size.x, size.y, size.z};
		int g[1] = {vc};

		real* pM = plsM->GetBlockData(block);
		int* pMaskId = plsMaskId->GetBlockData(block);

		setup_mask_(
			pM, 
			pMaskId, 
			sz, g);
	}
/* ---------------------------------------------------------- */

/* ---------------------------------------------------------- */
/* Init vars                                                  */
/* ---------------------------------------------------------- */
	int boundaryTypeUX[NUM_FACE] = {
		g_pFFVConfig->OuterBCUX[X_M].type,
		g_pFFVConfig->OuterBCUX[X_P].type,
		g_pFFVConfig->OuterBCUX[Y_M].type,
		g_pFFVConfig->OuterBCUX[Y_P].type,
		g_pFFVConfig->OuterBCUX[Z_M].type,
		g_pFFVConfig->OuterBCUX[Z_P].type,
	};
	real boundaryValueUX[NUM_FACE] = {
		g_pFFVConfig->OuterBCUX[X_M].value,
		g_pFFVConfig->OuterBCUX[X_P].value,
		g_pFFVConfig->OuterBCUX[Y_M].value,
		g_pFFVConfig->OuterBCUX[Y_P].value,
		g_pFFVConfig->OuterBCUX[Z_M].value,
		g_pFFVConfig->OuterBCUX[Z_P].value,
	};
	plsUX0 = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeUX, boundaryValueUX, 1);
	plsUX1 = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeUX, boundaryValueUX, 1);
	plsUXC = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsUXCP= new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsUXD = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsUXDP= new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	double ux0 = g_pFFVConfig->InitialValueU.x;
	double dux0 = g_pFFVConfig->InitialValueDU.x;
	plsUX0->Fill(blockManager, ux0, dux0);
	plsUX1->Fill(blockManager, ux0, dux0);
	plsUXC->Fill(blockManager, 0.0);
	plsUXCP->Fill(blockManager, 0.0);
	plsUXD->Fill(blockManager, 0.0);
	plsUXDP->Fill(blockManager, 0.0);
	plsUX0->ImposeBoundaryCondition(blockManager);
	plsUX1->ImposeBoundaryCondition(blockManager);
	plsUXC->ImposeBoundaryCondition(blockManager);
	plsUXCP->ImposeBoundaryCondition(blockManager);
	plsUXD->ImposeBoundaryCondition(blockManager);
	plsUXDP->ImposeBoundaryCondition(blockManager);

	int boundaryTypeUY[NUM_FACE] = {
		g_pFFVConfig->OuterBCUY[X_M].type,
		g_pFFVConfig->OuterBCUY[X_P].type,
		g_pFFVConfig->OuterBCUY[Y_M].type,
		g_pFFVConfig->OuterBCUY[Y_P].type,
		g_pFFVConfig->OuterBCUY[Z_M].type,
		g_pFFVConfig->OuterBCUY[Z_P].type,
	};
	real boundaryValueUY[NUM_FACE] = {
		g_pFFVConfig->OuterBCUY[X_M].value,
		g_pFFVConfig->OuterBCUY[X_P].value,
		g_pFFVConfig->OuterBCUY[Y_M].value,
		g_pFFVConfig->OuterBCUY[Y_P].value,
		g_pFFVConfig->OuterBCUY[Z_M].value,
		g_pFFVConfig->OuterBCUY[Z_P].value,
	};
	plsUY0 = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeUY, boundaryValueUY, 1);
	plsUY1 = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeUY, boundaryValueUY, 1);
	plsUYC = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsUYCP= new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsUYD = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsUYDP= new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	double uy0 = g_pFFVConfig->InitialValueU.y;
	double duy0 = g_pFFVConfig->InitialValueDU.y;
	plsUY0->Fill(blockManager, uy0, duy0);
	plsUY1->Fill(blockManager, uy0, duy0);
	plsUYC->Fill(blockManager, 0.0);
	plsUYCP->Fill(blockManager, 0.0);
	plsUYD->Fill(blockManager, 0.0);
	plsUYDP->Fill(blockManager, 0.0);
	plsUY0->ImposeBoundaryCondition(blockManager);
	plsUY1->ImposeBoundaryCondition(blockManager);
	plsUYC->ImposeBoundaryCondition(blockManager);
	plsUYCP->ImposeBoundaryCondition(blockManager);
	plsUYD->ImposeBoundaryCondition(blockManager);
	plsUYDP->ImposeBoundaryCondition(blockManager);

	int boundaryTypeUZ[NUM_FACE] = {
		g_pFFVConfig->OuterBCUZ[X_M].type,
		g_pFFVConfig->OuterBCUZ[X_P].type,
		g_pFFVConfig->OuterBCUZ[Y_M].type,
		g_pFFVConfig->OuterBCUZ[Y_P].type,
		g_pFFVConfig->OuterBCUZ[Z_M].type,
		g_pFFVConfig->OuterBCUZ[Z_P].type,
	};
	real boundaryValueUZ[NUM_FACE] = {
		g_pFFVConfig->OuterBCUZ[X_M].value,
		g_pFFVConfig->OuterBCUZ[X_P].value,
		g_pFFVConfig->OuterBCUZ[Y_M].value,
		g_pFFVConfig->OuterBCUZ[Y_P].value,
		g_pFFVConfig->OuterBCUZ[Z_M].value,
		g_pFFVConfig->OuterBCUZ[Z_P].value,
	};
	plsUZ0 = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeUZ, boundaryValueUZ, 1);
	plsUZ1 = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeUZ, boundaryValueUZ, 1);
	plsUZC = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsUZCP= new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsUZD = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsUZDP= new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	double uz0 = g_pFFVConfig->InitialValueU.z;
	double duz0 = g_pFFVConfig->InitialValueDU.z;
	plsUZ0->Fill(blockManager, uz0, duz0);
	plsUZ1->Fill(blockManager, uz0, duz0);
	plsUZC->Fill(blockManager, 0.0);
	plsUZCP->Fill(blockManager, 0.0);
	plsUZD->Fill(blockManager, 0.0);
	plsUZDP->Fill(blockManager, 0.0);
	plsUZ0->ImposeBoundaryCondition(blockManager);
	plsUZ1->ImposeBoundaryCondition(blockManager);
	plsUZC->ImposeBoundaryCondition(blockManager);
	plsUZCP->ImposeBoundaryCondition(blockManager);
	plsUZD->ImposeBoundaryCondition(blockManager);
	plsUZDP->ImposeBoundaryCondition(blockManager);

	plsVw = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsVe = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsVs = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsVn = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsVb = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsVt = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsVw->Fill(blockManager, 0.0);
	plsVe->Fill(blockManager, 0.0);
	plsVs->Fill(blockManager, 0.0);
	plsVn->Fill(blockManager, 0.0);
	plsVb->Fill(blockManager, 0.0);
	plsVt->Fill(blockManager, 0.0);
	plsVw->ImposeBoundaryCondition(blockManager);
	plsVe->ImposeBoundaryCondition(blockManager);
	plsVs->ImposeBoundaryCondition(blockManager);
	plsVn->ImposeBoundaryCondition(blockManager);
	plsVb->ImposeBoundaryCondition(blockManager);
	plsVt->ImposeBoundaryCondition(blockManager);

	int boundaryTypeT[NUM_FACE] = {
		g_pFFVConfig->OuterBCT[X_M].type,
		g_pFFVConfig->OuterBCT[X_P].type,
		g_pFFVConfig->OuterBCT[Y_M].type,
		g_pFFVConfig->OuterBCT[Y_P].type,
		g_pFFVConfig->OuterBCT[Z_M].type,
		g_pFFVConfig->OuterBCT[Z_P].type,
	};
	real boundaryValueT[NUM_FACE] = {
		g_pFFVConfig->OuterBCT[X_M].value,
		g_pFFVConfig->OuterBCT[X_P].value,
		g_pFFVConfig->OuterBCT[Y_M].value,
		g_pFFVConfig->OuterBCT[Y_P].value,
		g_pFFVConfig->OuterBCT[Z_M].value,
		g_pFFVConfig->OuterBCT[Z_P].value,
	};
	plsT0 = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeT, boundaryValueT, 1);
	plsT1 = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeT, boundaryValueT, 1);
	plsTC = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsTCP= new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsTD = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsTDP= new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	double t0 = g_pFFVConfig->InitialValueT;
	plsT0->Fill(blockManager, t0);
	plsT1->Fill(blockManager, t0);
	plsTC->Fill(blockManager, 0.0);
	plsTCP->Fill(blockManager, 0.0);
	plsTD->Fill(blockManager, 0.0);
	plsTDP->Fill(blockManager, 0.0);
	plsT0->ImposeBoundaryCondition(blockManager);
	plsT1->ImposeBoundaryCondition(blockManager);
	plsTC->ImposeBoundaryCondition(blockManager);
	plsTCP->ImposeBoundaryCondition(blockManager);
	plsTD->ImposeBoundaryCondition(blockManager);
	plsTDP->ImposeBoundaryCondition(blockManager);

	int boundaryTypeP[NUM_FACE] = {
		g_pFFVConfig->OuterBCP[X_M].type,
		g_pFFVConfig->OuterBCP[X_P].type,
		g_pFFVConfig->OuterBCP[Y_M].type,
		g_pFFVConfig->OuterBCP[Y_P].type,
		g_pFFVConfig->OuterBCP[Z_M].type,
		g_pFFVConfig->OuterBCP[Z_P].type,
	};
	real boundaryValueP[NUM_FACE] = {
		g_pFFVConfig->OuterBCP[X_M].value,
		g_pFFVConfig->OuterBCP[X_P].value,
		g_pFFVConfig->OuterBCP[Y_M].value,
		g_pFFVConfig->OuterBCP[Y_P].value,
		g_pFFVConfig->OuterBCP[Z_M].value,
		g_pFFVConfig->OuterBCP[Z_P].value,
	};
	plsP0 = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeP, boundaryValueP, 1);
	plsP1 = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeP, boundaryValueP, 1);
	plsLapP = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	double p0 = g_pFFVConfig->InitialValueP;
	plsP0->Fill(blockManager, p0);
	plsP1->Fill(blockManager, p0);
	plsLapP->Fill(blockManager, 0.0);
	plsP0->ImposeBoundaryCondition(blockManager);
	plsP1->ImposeBoundaryCondition(blockManager);
	plsLapP->ImposeBoundaryCondition(blockManager);
/* ---------------------------------------------------------- */


/* ---------------------------------------------------------- */
/* Init A and b                                               */
/* ---------------------------------------------------------- */
	plsAp = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsAw = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsAe = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsAs = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsAn = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsAb = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsAt = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsb  = new LocalScalar3D<real>(blockManager, vc, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsAp->Fill(blockManager, 0.0);
	plsAw->Fill(blockManager, 0.0);
	plsAe->Fill(blockManager, 0.0);
	plsAs->Fill(blockManager, 0.0);
	plsAn->Fill(blockManager, 0.0);
	plsAb->Fill(blockManager, 0.0);
	plsAt->Fill(blockManager, 0.0);
	plsb ->Fill(blockManager, 0.0);
	plsAp->ImposeBoundaryCondition(blockManager);
	plsAw->ImposeBoundaryCondition(blockManager);
	plsAe->ImposeBoundaryCondition(blockManager);
	plsAs->ImposeBoundaryCondition(blockManager);
	plsAn->ImposeBoundaryCondition(blockManager);
	plsAb->ImposeBoundaryCondition(blockManager);
	plsAt->ImposeBoundaryCondition(blockManager);
	plsb ->ImposeBoundaryCondition(blockManager);
/* ---------------------------------------------------------- */

/* ---------------------------------------------------------- */
/* Init ILS                                                   */
/* ---------------------------------------------------------- */
	omegaU		= g_pFFVConfig->IterationOmegaU;
	countMaxU	= g_pFFVConfig->IterationMaxCountU;
	epsilonU	= g_pFFVConfig->IterationEpsilonU;
	countPreConditionerU
						= g_pFFVConfig->IterationPreCountU;
	countUX		= 0;
	residualUX= 0.0;
	countUY		= 0;
	residualUY= 0.0;
	countUZ		= 0;
	residualUZ= 0.0;

	omegaP		= g_pFFVConfig->IterationOmegaP;
	countMaxP	= g_pFFVConfig->IterationMaxCountP;
	epsilonP	= g_pFFVConfig->IterationEpsilonP;
	countPreConditionerP
						= g_pFFVConfig->IterationPreCountP;
	countP		= 0;
	residualP	= 0.0;

	omegaT		= g_pFFVConfig->IterationOmegaT;
	countMaxT	= g_pFFVConfig->IterationMaxCountT;
	epsilonT	= g_pFFVConfig->IterationEpsilonT;
	countPreConditionerT
						= g_pFFVConfig->IterationPreCountT;
	countT		= 0;
	residualT	= 0.0;

	pils = new FFVILS();

	int vc1 = 2;
	plsr  = new LocalScalar3D<real>(blockManager, vc1, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsr0 = new LocalScalar3D<real>(blockManager, vc1, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsp  = new LocalScalar3D<real>(blockManager, vc1, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsp_ = new LocalScalar3D<real>(blockManager, vc1, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsq_ = new LocalScalar3D<real>(blockManager, vc1, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plss  = new LocalScalar3D<real>(blockManager, vc1, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plss_ = new LocalScalar3D<real>(blockManager, vc1, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plst_ = new LocalScalar3D<real>(blockManager, vc1, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsr ->Fill(blockManager, 0.0);
	plsr0->Fill(blockManager, 0.0);
	plsp ->Fill(blockManager, 0.0);
	plsp_->Fill(blockManager, 0.0);
	plsq_->Fill(blockManager, 0.0);
	plss ->Fill(blockManager, 0.0);
	plss_->Fill(blockManager, 0.0);
	plst_->Fill(blockManager, 0.0);
/* ---------------------------------------------------------- */


/* ---------------------------------------------------------- */
/* Init Force                                                 */
/* ---------------------------------------------------------- */
	int vc2 = 2;
	plsFspx = new LocalScalar3D<real>(blockManager, vc2, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsFspy = new LocalScalar3D<real>(blockManager, vc2, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsFspz = new LocalScalar3D<real>(blockManager, vc2, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsFsvx = new LocalScalar3D<real>(blockManager, vc2, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsFsvy = new LocalScalar3D<real>(blockManager, vc2, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsFsvz = new LocalScalar3D<real>(blockManager, vc2, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsFspx->Fill(blockManager, 0.0);
	plsFspy->Fill(blockManager, 0.0);
	plsFspz->Fill(blockManager, 0.0);
	plsFsvx->Fill(blockManager, 0.0);
	plsFsvy->Fill(blockManager, 0.0);
	plsFsvz->Fill(blockManager, 0.0);
/* ---------------------------------------------------------- */


/* ---------------------------------------------------------- */
/* Init Heat flux                                             */
/* ---------------------------------------------------------- */
	int vc3 = 2;
	plsQx = new LocalScalar3D<real>(blockManager, vc3, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsQy = new LocalScalar3D<real>(blockManager, vc3, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsQz = new LocalScalar3D<real>(blockManager, vc3, updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsQx->Fill(blockManager, 0.0);
	plsQy->Fill(blockManager, 0.0);
	plsQz->Fill(blockManager, 0.0);
/* ---------------------------------------------------------- */

	PM_Stop(tm_Init_InitVars);

/* ---------------------------------------------------------- */
	PrintLog(2, "Completed");
/* ---------------------------------------------------------- */


/////////////////////////////////////////////
	if( g_pFFVConfig->OutputDataBasicVariablesFormatBCM ) {
		BCMFileSaverInit(
					rootGrid,
					tree,
					partition);
	}
	if( g_pFFVConfig->OutputDataBasicVariablesFormatPLOT3D ) {
		WriteXYZInPlot3DFormat(
					"xyz",
					diffLevel,
					rootGrid,
					tree,
					partition);
	}
/////////////////////////////////////////////



	for(int i=0; i<32; i++) {
		this->times[i] = 0.0;
	}

	return EX_SUCCESS;
}

int Solver::Loop() {
	int StepStart = g_pFFVConfig->TimeControlSessionStartI;
	int StepEnd   = g_pFFVConfig->TimeControlSessionEndI;

	bRestart = false;
	if(StepStart == 0) {
		Print(0);
	} else {
		Load(StepStart);
		bRestart = true;
	}

	for(int step=StepStart+1; step<=StepEnd; step++) {
		int nResult = Update(step);

		Print(step);

		if( step%g_pFFVConfig->RestartInterval == 0 ) {
			Dump(step);
		}

		switch(nResult) {
			case EX_SUCCESS:
				break;
			default:
				break;
		}
	}

	return EX_SUCCESS;
}


int Solver::Print(int step) {
PM_Start(tm_Print, 0, 0, true);

	if( g_pFFVConfig->OutputLogLaptime ) {
PM_Start(tm_PrintTime, 0, 0, true);
		PrintTime(step);
PM_Stop(tm_PrintTime);
	}

	if( g_pFFVConfig->OutputLogIteration ) {
PM_Start(tm_PrintILS, 0, 0, true);
		PrintILS(step);
PM_Stop(tm_PrintILS);
	}

	if( g_pFFVConfig->OutputLogStatistics ) {
PM_Start(tm_PrintStats, 0, 0, true);
		PrintStats(step);
PM_Stop(tm_PrintStats);
	}

	if( g_pFFVConfig->OutputLogForce ) {
PM_Start(tm_PrintForce, 0, 0, true);
		PrintForce(step);
PM_Stop(tm_PrintForce);
	}

	if( g_pFFVConfig->OutputLogHeatFlux ) {
PM_Start(tm_PrintHeatFlux, 0, 0, true);
		PrintHeatFlux(step);
PM_Stop(tm_PrintHeatFlux);
	}

PM_Start(tm_PrintData, 0, 0, true);
	PrintData(step);
PM_Stop(tm_PrintData);

PM_Stop(tm_Print);

	return EX_SUCCESS;
}

void Solver::PrintData(int step) {
	if( g_pFFVConfig->OutputDataBasicVariablesIntervalI > 0 ) {
		if( step%g_pFFVConfig->OutputDataBasicVariablesIntervalI == 0 ) {
			if( g_pFFVConfig->OutputDataBasicVariablesFormatVTK ) {
				PrintBasicVariablesVTK(step);
			}
			if( g_pFFVConfig->OutputDataBasicVariablesFormatPLOT3D ) {
				PrintBasicVariablesPLOT3D(step);
			}
			if( g_pFFVConfig->OutputDataBasicVariablesFormatBCM ) {
				PrintBasicVariablesBCM(step);
			}
			if( g_pFFVConfig->OutputDataBasicVariablesFormatSILO ) {
				PrintBasicVariablesSILO(step);
			}
		}
	}

	if( g_pFFVConfig->OutputDataDerivedVariablesIntervalI > 0 ) {
		if( step%g_pFFVConfig->OutputDataDerivedVariablesIntervalI == 0 ) {
			if( g_pFFVConfig->OutputDataDerivedVariablesFormatVTK ) {
				PrintDerivedVariablesVTK(step);
			}
			if( g_pFFVConfig->OutputDataDerivedVariablesFormatPLOT3D ) {
				PrintDerivedVariablesPLOT3D(step);
			}
			if( g_pFFVConfig->OutputDataDerivedVariablesFormatBCM ) {
				PrintDerivedVariablesBCM(step);
			}
			if( g_pFFVConfig->OutputDataDerivedVariablesFormatSILO ) {
				PrintDerivedVariablesSILO(step);
			}
		}
	}
}

void Solver::PrintTime(int step) {
	if( step%g_pFFVConfig->OutputLogFileIntervalLaptime == 0) {
	} else {
		return;
	}

	if( myrank != 0 ) {
		return;
	}

	ostringstream ossMessage;
	ossMessage.width(10);
	ossMessage.setf(std::ios::fixed);
	ossMessage.fill('0');
	ossMessage << step << " ";
	ossMessage.setf(std::ios::scientific, std::ios::floatfield);
	ossMessage.precision(16);
	ossMessage << times[0];
	PrintLog(0, ossMessage.str().c_str());

	std::string filename = "data-times.txt";

	std::ofstream ofs;
	if( step==0 ) {
		ofs.open(filename.c_str(), std::ios::out);
		ofs.close();
	}
	ofs.open(filename.c_str(), std::ios::out | std::ios::app);

	ofs.width(10);
	ofs.setf(std::ios::fixed);
	ofs.fill('0');
	ofs << step << " ";

	ofs.setf(std::ios::scientific, std::ios::floatfield);
	ofs.precision(16);
	ofs << this->times[0] << " ";
	ofs << this->times[1] << " ";
	ofs << this->times[2] << " ";
	ofs << this->times[3] << " ";
	ofs << this->times[4] << " ";
	ofs << this->times[5] << " ";
	ofs << this->times[6] << " ";
	ofs << std::endl;
	ofs.close();
}

void Solver::PrintILS(int step) {
	if( step%g_pFFVConfig->OutputLogFileIntervalIteration == 0 ) {
	} else {
		return;
	}

	if( myrank != 0 ) {
		return;
	}

	std::string filename = "data-ils.txt";

	std::ofstream ofs;
	if( step==0 ) {
		ofs.open(filename.c_str(), std::ios::out);
		ofs.close();
	}
	ofs.open(filename.c_str(), std::ios::out | std::ios::app);

	ofs.width(10);
	ofs.setf(std::ios::fixed);
	ofs.fill('0');
	ofs << step << " ";

	ofs.setf(std::ios::scientific, std::ios::floatfield);
	ofs.width(4);
	ofs.setf(std::ios::fixed);
	ofs.fill('0');
	ofs << this->countUX << " ";
	ofs.precision(16);
	ofs.setf(std::ios::fixed);
	ofs.fill('0');
	ofs << this->residualUX << " ";
	ofs.width(4);
	ofs.setf(std::ios::fixed);
	ofs.fill('0');
	ofs << this->countUY << " ";
	ofs.precision(16);
	ofs.setf(std::ios::fixed);
	ofs.fill('0');
	ofs << this->residualUY << " ";
	ofs.width(4);
	ofs.setf(std::ios::fixed);
	ofs.fill('0');
	ofs << this->countUZ << " ";
	ofs.precision(16);
	ofs.setf(std::ios::fixed);
	ofs.fill('0');
	ofs << this->residualUZ << " ";
	ofs.width(4);
	ofs.setf(std::ios::fixed);
	ofs.fill('0');
	ofs << this->countP << " ";
	ofs.precision(16);
	ofs.setf(std::ios::fixed);
	ofs.fill('0');
	ofs << this->residualP << " ";
	ofs.width(4);
	ofs.setf(std::ios::fixed);
	ofs.fill('0');
	ofs << this->countT << " ";
	ofs.precision(16);
	ofs.setf(std::ios::fixed);
	ofs.fill('0');
	ofs << this->residualT << " ";
	ofs << std::endl;
	ofs.close();
}

void Solver::PrintStats(int step) {
	if( step%g_pFFVConfig->OutputLogFileIntervalStatistics == 0 ) {
	} else {
		return;
	}

	plsUX0->CalcStats(blockManager);
	plsUY0->CalcStats(blockManager);
	plsUZ0->CalcStats(blockManager);
	plsP0->CalcStats(blockManager);
	plsT0->CalcStats(blockManager);

	if( myrank == 0 ) {
		std::string filename = "data-stats.txt";

		std::ofstream ofs;
		if( step==0 ) {
			ofs.open(filename.c_str(), std::ios::out);
			ofs.close();
		}
		ofs.open(filename.c_str(), std::ios::out | std::ios::app);

		ofs.width(10);
		ofs.setf(std::ios::fixed);
		ofs.fill('0');
		ofs << step << " ";

		ofs.setf(std::ios::scientific, std::ios::floatfield);
		ofs.precision(16);
		ofs << plsUX0->GetSum() << " ";
		ofs << plsUY0->GetSum() << " ";
		ofs << plsUZ0->GetSum() << " ";
		ofs << plsP0->GetSum() << " ";
		ofs << plsT0->GetSum() << " ";
		ofs << plsUX0->GetMax() << " ";
		ofs << plsUY0->GetMax() << " ";
		ofs << plsUZ0->GetMax() << " ";
		ofs << plsP0->GetMax() << " ";
		ofs << plsT0->GetMax() << " ";
		ofs << plsUX0->GetMin() << " ";
		ofs << plsUY0->GetMin() << " ";
		ofs << plsUZ0->GetMin() << " ";
		ofs << plsP0->GetMin() << " ";
		ofs << plsT0->GetMin() << " ";
		ofs << plsUX0->GetAbsMax() << " ";
		ofs << plsUY0->GetAbsMax() << " ";
		ofs << plsUZ0->GetAbsMax() << " ";
		ofs << plsP0->GetAbsMax() << " ";
		ofs << plsT0->GetAbsMax() << " ";
		ofs << plsUX0->GetAbsMin() << " ";
		ofs << plsUY0->GetAbsMin() << " ";
		ofs << plsUZ0->GetAbsMin() << " ";
		ofs << plsP0->GetAbsMin() << " ";
		ofs << plsT0->GetAbsMin() << " ";
		ofs << std::endl;
		ofs.close();
	}

	double dt = g_pFFVConfig->TimeControlTimeStepDeltaT;
	double dx = g_pFFVConfig->RootBlockLength/(double)(1 << maxLevel);
//	std::cout << dt << " " << dx << std::endl;
	if( plsUX0->GetAbsMaxL()*dt/dx > 0.9 ) {
		std::cout << myrank << " " << plsUX0->GetAbsMaxL()*dt/dx << std::endl;
	}

	if( plsUX0->GetAbsMax()*dt/dx > 0.9 ) {
		if( myrank==0 ) {
			std::cout << myrank << " " << plsUX0->GetAbsMax()*dt/dx << std::endl;
		}
		PrintData(step);
		exit(EX_FAILURE);
	}
}

void Solver::PrintHeatFlux(int step) {
	if( step%g_pFFVConfig->OutputLogFileIntervalHeatFlux == 0 ) {
	} else {
		return;
	}

	int MaxCID = 32;
	for(int n=1; n<MaxCID; n++) {
		PrintHeatFluxCID(step, n);
	}
}

void Solver::PrintHeatFluxCID(int step, int cid_target) {
	int bc_n[1] = {32};
	int bc_type[32];
	real bc_value[32];
	for(int m=0; m<bc_n[0]; m++) {
		bc_type[m]  = g_pFFVConfig->BCInternalBoundaryType[m];
		bc_value[m] = g_pFFVConfig->BCInternalBoundaryValue[m];
	}

	real q_local_x = 0.0;
	real q_local_y = 0.0;
	real q_local_z = 0.0;
	real sa_local = 0.0;
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for reduction(+:q_local_x, q_local_y, q_local_z, sa_local)
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
		
		real* t0   = plsT0 ->GetBlockData(block);

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

		int* pPhaseId = plsPhaseId->GetBlockData(block);

		real* qx = plsQx->GetBlockData(block);
		real* qy = plsQy->GetBlockData(block);
		real* qz = plsQz->GetBlockData(block);

		int* pNormalIndex0 = plsNormalIndex0->GetBlockData(block);
		int* pNormalIndex1 = plsNormalIndex1->GetBlockData(block);
		int* pNormalIndex2 = plsNormalIndex2->GetBlockData(block);
		int* pNormalIndex3 = plsNormalIndex3->GetBlockData(block);
		int* pNormalIndex4 = plsNormalIndex4->GetBlockData(block);
		int* pNormalIndex5 = plsNormalIndex5->GetBlockData(block);

		real q_block[3] = {0.0, 0.0, 0.0};
		real sa_block[1] = {0.0};

		bcut_calc_q_(
				qx,
				qy,
				qz,
				q_block,
				sa_block,
				&cid_target,
				t0,
				pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
				pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
				pPhaseId,
				&pNormalN[n],
				pNormalX[n],
				pNormalY[n],
				pNormalZ[n],
				pNormalIndex0, pNormalIndex1, pNormalIndex2, pNormalIndex3, pNormalIndex4, pNormalIndex5,
				&rhof, &rhos,
				&cpf, &cps,
				&kf, &ks,
				bc_n,
				bc_type,
				bc_value,
				org,
				&dx, &dt,
				sz, g);

		q_local_x += q_block[0];
		q_local_y += q_block[1];
		q_local_z += q_block[2];
		sa_local  += sa_block[0];
	}

	real q_global[3] = {0.0, 0.0, 0.0};
	real sa_global = 0.0;

	real sum = q_local_x;
	real sum_tmp = sum;
	comm_sum_(&sum, &sum_tmp);
	q_global[0] = sum;

	sum = q_local_y;
	sum_tmp = sum;
	comm_sum_(&sum, &sum_tmp);
	q_global[1] = sum;

	sum = q_local_z;
	sum_tmp = sum;
	comm_sum_(&sum, &sum_tmp);
	q_global[2] = sum;

	sum = sa_local;
	sum_tmp = sum;
	comm_sum_(&sum, &sum_tmp);
	sa_global = sum;

	if( myrank != 0 ) {
		return;
	}

	std::ostringstream ossFileName;
	ossFileName << "data-heatflux";
	ossFileName << "-";
	ossFileName.width(3);
	ossFileName.setf(std::ios::fixed);
	ossFileName.fill('0');
	ossFileName << cid_target;
	ossFileName << ".txt";

	std::string filename = ossFileName.str();

	std::ofstream ofs;
	if( step==0 ) {
		ofs.open(filename.c_str(), std::ios::out);
		ofs.close();
	}
	ofs.open(filename.c_str(), std::ios::out | std::ios::app);

	ofs.width(10);
	ofs.setf(std::ios::fixed);
	ofs.fill('0');
	ofs << step << " ";

	ofs.setf(std::ios::scientific, std::ios::floatfield);
	ofs.precision(16);
	ofs << q_global[0] << " ";
	ofs << q_global[1] << " ";
	ofs << q_global[2] << " ";
	ofs << sa_global << " ";
	ofs << std::endl;
	ofs.close();
}

void Solver::PrintForce(int step) {
	if( step%g_pFFVConfig->OutputLogFileIntervalForce == 0 ) {
	} else {
		return;
	}

	int MaxCID = 32;
	for(int n=1; n<MaxCID; n++) {
		PrintForceCID(step, n);
	}
}

void Solver::PrintForceCID(int step, int cid_target) {
	real fsp_local[3] = {0.0, 0.0, 0.0};
	real fsv_local[3] = {0.0, 0.0, 0.0};
	real fsp_local_x = 0.0;
	real fsp_local_y = 0.0;
	real fsp_local_z = 0.0;
	real fsv_local_x = 0.0;
	real fsv_local_y = 0.0;
	real fsv_local_z = 0.0;
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for reduction(+:fsp_local_x, fsp_local_y, fsp_local_z, fsv_local_x, fsv_local_y, fsv_local_z)
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
		
		real* ux0  = plsUX0->GetBlockData(block);
		real* uy0  = plsUY0->GetBlockData(block);
		real* uz0  = plsUZ0->GetBlockData(block);
		real* p0   = plsP0 ->GetBlockData(block);

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

		int* pPhaseId = plsPhaseId->GetBlockData(block);

		real* fspx = plsFspx->GetBlockData(block);
		real* fspy = plsFspy->GetBlockData(block);
		real* fspz = plsFspz->GetBlockData(block);
		real* fsvx = plsFsvx->GetBlockData(block);
		real* fsvy = plsFsvy->GetBlockData(block);
		real* fsvz = plsFsvz->GetBlockData(block);

		real Uc = 0.0;

		real fsp_block[3] = {0.0, 0.0, 0.0};
		real fsv_block[3] = {0.0, 0.0, 0.0};

		bcut_calc_f_p_(
				fspx,
				fspy,
				fspz,
				fsp_block,
				&cid_target,
				p0,
				pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
				pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
				pPhaseId,
				&dx, &dt,
				sz, g);

		bcut_calc_f_v_(
				fsvx,
				fsvy,
				fsvz,
				fsv_block,
				&cid_target,
				ux0,
				uy0,
				uz0,
				pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
				pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
				pPhaseId,
				&rhof,
				&mu,
				&dx, &dt,
				&Uc,
				sz, g);

		fsp_local_x += fsp_block[0];
		fsp_local_y += fsp_block[1];
		fsp_local_z += fsp_block[2];
		fsv_local_x += fsv_block[0];
		fsv_local_y += fsv_block[1];
		fsv_local_z += fsv_block[2];
	}

	real fsp_global[3] = {0.0, 0.0, 0.0};
	real fsv_global[3] = {0.0, 0.0, 0.0};

	real sum = fsp_local_x;
	real sum_tmp = sum;
	comm_sum_(&sum, &sum_tmp);
	fsp_global[0] = sum;

	sum = fsp_local_y;
	sum_tmp = sum;
	comm_sum_(&sum, &sum_tmp);
	fsp_global[1] = sum;

	sum = fsp_local_z;
	sum_tmp = sum;
	comm_sum_(&sum, &sum_tmp);
	fsp_global[2] = sum;

	sum = fsv_local_x;
	sum_tmp = sum;
	comm_sum_(&sum, &sum_tmp);
	fsv_global[0] = sum;

	sum = fsv_local_y;
	sum_tmp = sum;
	comm_sum_(&sum, &sum_tmp);
	fsv_global[1] = sum;

	sum = fsv_local_z;
	sum_tmp = sum;
	comm_sum_(&sum, &sum_tmp);
	fsv_global[2] = sum;

	if( myrank != 0 ) {
		return;
	}

	std::ostringstream ossFileName;
	ossFileName << "data-force";
	ossFileName << "-";
	ossFileName.width(3);
	ossFileName.setf(std::ios::fixed);
	ossFileName.fill('0');
	ossFileName << cid_target;
	ossFileName << ".txt";

	std::string filename = ossFileName.str();

	std::ofstream ofs;
	if( step==0 ) {
		ofs.open(filename.c_str(), std::ios::out);
		ofs.close();
	}
	ofs.open(filename.c_str(), std::ios::out | std::ios::app);

	ofs.width(10);
	ofs.setf(std::ios::fixed);
	ofs.fill('0');
	ofs << step << " ";

	ofs.setf(std::ios::scientific, std::ios::floatfield);
	ofs.precision(16);
	ofs << fsp_global[0] << " ";
	ofs << fsp_global[1] << " ";
	ofs << fsp_global[2] << " ";
	ofs << fsv_global[0] << " ";
	ofs << fsv_global[1] << " ";
	ofs << fsv_global[2] << " ";
	ofs << std::endl;
	ofs.close();
}

void Solver::PrintBasicVariablesVTK(int step) {
	WriteBasicVariablesInVTKFormat(step, diffLevel, rootGrid, tree, partition);
	if( step == 0 ) {
		plsPhaseId->WriteDataInVTKFormat("phase", 0, diffLevel, rootGrid, tree, partition);
	}
}

void Solver::PrintBasicVariablesPLOT3D(int step) {
	WriteBasicVariablesInPlot3DFormat("flow", step, diffLevel, rootGrid, tree, partition);
}

void Solver::PrintBasicVariablesBCM(int step) {
	BCMFileSaverPrint(step);
}

void Solver::PrintBasicVariablesSILO(int step) {
	WriteBasicVariablesInSILOFormat(step, diffLevel, rootGrid, tree, partition);
}

void Solver::PrintDerivedVariablesVTK(int step) {
	WriteDerivedVariablesInVTKFormat(step, diffLevel, rootGrid, tree, partition);
}

void Solver::PrintDerivedVariablesPLOT3D(int step) {
}

void Solver::PrintDerivedVariablesBCM(int step) {
}

void Solver::PrintDerivedVariablesSILO(int step) {
}

void Solver::PrintLog(int level, const char* format, ...) {
	if( myrank != 0 ) {
		return;
	}
	static int first = 1;

	va_list ap;
	va_start(ap, format);

	char* buffer;
	int size = vasprintf(&buffer, format, ap);
	va_end(ap);

	ostringstream ossMessage;
	switch(level) {
		case 1:
			ossMessage << "# ";
			break;
		case 2:
			ossMessage << "#   ";
			break;
		case 3:
			ossMessage << "#     ";
			break;
		case 0:
		default:
			break;
	}
	ossMessage << string(buffer);


	free(buffer);

	std::cout << ossMessage.str() << std::endl;

	std::string filename = "data-log.txt";
	std::ofstream ofs;
	if( first == 1 ) {
		ofs.open(filename.c_str(), std::ios::out);
		ofs.close();
	}
	ofs.open(filename.c_str(), std::ios::out | std::ios::app);
	ofs << ossMessage.str() << std::endl;
	ofs.close();

	first = 0;
}

void Solver::PrintCS(int step) {
	mkdir("CS", 0755);

	std::ostringstream ossFileName;
	ossFileName << "./CS/";
	ossFileName << "data-cs2";
	ossFileName << "-";
	ossFileName.width(5);
	ossFileName.setf(std::ios::fixed);
	ossFileName.fill('0');
	ossFileName << myrank;
	ossFileName << "-";
	ossFileName.width(10);
	ossFileName.setf(std::ios::fixed);
	ossFileName.fill('0');
	ossFileName << step;
	ossFileName << ".txt";

	std::ofstream ofs;
	ofs.open(ossFileName.str().c_str(), std::ios::out);
	for (int n=0; n<blockManager.getNumBlock(); ++n) {
		BlockBase* block = blockManager.getBlock(n);
		Vec3i size = block->getSize();
		Vec3r origin = block->getOrigin();
		Vec3r blockSize = block->getBlockSize();
		Vec3r cellSize = block->getCellSize();

		int sz[3] = {size.x, size.y, size.z};
		int g[1] = {vc};
		int nc[3] = {size.x + 2*vc, size.y + 2*vc, size.z + 2*vc};

		real* pUX0 = plsUX0->GetBlockData(block);
		real* pUX1 = plsUX1->GetBlockData(block);
		real* pUY0 = plsUY0->GetBlockData(block);
		real* pUZ0 = plsUZ0->GetBlockData(block);
		real* pP0 = plsP0->GetBlockData(block);

		int j = size.y/2 + vc;
		for(int k=vc; k<vc+size.z; k++) {
		for(int i=vc; i<vc+size.x; i++) {
			double x = (i - vc + 0.5)*cellSize.x + origin.x;
			double y = (j - vc + 0.5)*cellSize.y + origin.y;
			double z = (k - vc + 0.5)*cellSize.z + origin.z;
			int m = i + nc[0]*(j + nc[1]*k);
			ofs.setf(std::ios::scientific, std::ios::floatfield);
			ofs.precision(16);
			ofs << i - vc;
			ofs << " ";
			ofs << j - vc;
			ofs << " ";
			ofs << k - vc;
			ofs << " ";
			ofs << x;
			ofs << " ";
			ofs << y;
			ofs << " ";
			ofs << z;
			ofs << " ";
			ofs << pUX0[m];
			ofs << " ";
			ofs << pUX1[m];
			ofs << " ";
			ofs << pUY0[m];
			ofs << " ";
			ofs << pUZ0[m];
			ofs << " ";
			ofs << pP0[m];
			ofs << std::endl;
		}
		ofs << std::endl;
		}
	}
	ofs.close();
}

int Solver::Post() {
	g_pPM->gather();

	if( this->myrank == 0 ) {
		char hostname[MPI_MAX_PROCESSOR_NAME] = {0};
		int nameLen;
		MPI_Get_processor_name(hostname, &nameLen);
		nameLen++;

		FILE* fp = fopen("data-pm.txt", "w");
		g_pPM->print(fp, hostname, g_pFFVConfig->OperatorName);
		g_pPM->printDetail(fp);
		fclose(fp);
	}

	return EX_SUCCESS;
}

int Solver::Update(int step) {
PM_Start(tm_Update, 0, 0, true);
double t0 = GetTime();

	if( g_pFFVConfig->TimeControlAccelerationAcceleratingTimeI > 0 && step <= g_pFFVConfig->TimeControlAccelerationAcceleratingTimeI ) {
		real vb = g_pFFVConfig->OuterBCUX[X_M].value*(real)step/(real)g_pFFVConfig->TimeControlAccelerationAcceleratingTimeI;
		plsUX0->ResetBoundaryConditionValue(blockManager, 0, vb);
		plsUX1->ResetBoundaryConditionValue(blockManager, 0, vb);
	}

PM_Start(tm_UpdateT, 0, 0, true);
	if( !strcasecmp(g_pFFVConfig->TimeIntegrationMethodForFlow.c_str(), "explicit") ) {
		UpdateTe(step);
	} else {
		UpdateT(step);
	}
PM_Stop(tm_UpdateT);
double t1 = GetTime();

PM_Start(tm_UpdateUX, 0, 0, true);
	if( !strcasecmp(g_pFFVConfig->TimeIntegrationMethodForFlow.c_str(), "explicit") ) {
		UpdateUXe(step);
	} else {
		UpdateUX(step);
	}
PM_Stop(tm_UpdateUX);
double t2 = GetTime();

PM_Start(tm_UpdateUY, 0, 0, true);
	if( !strcasecmp(g_pFFVConfig->TimeIntegrationMethodForFlow.c_str(), "explicit") ) {
		UpdateUYe(step);
	} else {
		UpdateUY(step);
	}
PM_Stop(tm_UpdateUY);
double t3 = GetTime();

PM_Start(tm_UpdateUZ, 0, 0, true);
	if( !strcasecmp(g_pFFVConfig->TimeIntegrationMethodForFlow.c_str(), "explicit") ) {
		UpdateUZe(step);
	} else {
		UpdateUZ(step);
	}
PM_Stop(tm_UpdateUZ);
double t4 = GetTime();

PM_Start(tm_UpdateP, 0, 0, true);
	UpdateP(step);
PM_Stop(tm_UpdateP);
double t5 = GetTime();

PM_Start(tm_UpdateU, 0, 0, true);
	UpdateU(step);
PM_Stop(tm_UpdateU);
double t6 = GetTime();

PM_Stop(tm_Update);

	this->times[0] = t6 - t0;
	this->times[1] = t1 - t0;
	this->times[2] = t2 - t1;
	this->times[3] = t3 - t2;
	this->times[4] = t4 - t3;
	this->times[5] = t5 - t4;
	this->times[6] = t6 - t5;

	return EX_SUCCESS;
}

void Solver::UpdateUXe(int step) {
/////////////////////////////////////////////
// Calc A & b
/////////////////////////////////////////////
PM_Start(tm_UpdateUX01, 0, 0, true);
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
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
	
		real* vw = plsVw->GetBlockData(block);
		real* ve = plsVe->GetBlockData(block);
		real* vs = plsVs->GetBlockData(block);
		real* vn = plsVn->GetBlockData(block);
		real* vb = plsVb->GetBlockData(block);
		real* vt = plsVt->GetBlockData(block);

		real* ux0  = plsUX0->GetBlockData(block);
		real* uxc0 = plsUXC->GetBlockData(block);
		real* uxcp = plsUXCP->GetBlockData(block);
		real* uxd0 = plsUXD->GetBlockData(block);
		real* uxdp = plsUXDP->GetBlockData(block);
		real* p0   = plsP0->GetBlockData(block);

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

		int* pPhaseId = plsPhaseId->GetBlockData(block);

		real Uc = 0.0;

		real* uy0  = plsUY0->GetBlockData(block);
		real* uz0  = plsUZ0->GetBlockData(block);

		if( g_pFFVConfig->ConvectionTermScheme == "W3" ) {
			bcut_calc_c_f_w3_(
					uxc0,
					ux0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Uc,
					sz, g);
		} else if( g_pFFVConfig->ConvectionTermScheme == "E3" ) {
			bcut_calc_c_f_e3_(
					uxc0,
					ux0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Uc,
					sz, g);
		} else if( g_pFFVConfig->ConvectionTermScheme == "QUICK" ) {
			bcut_calc_c_f_quick_(
					uxc0,
					ux0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Uc,
					sz, g);
		} else if( g_pFFVConfig->ConvectionTermScheme == "QUICK2" ) {
			bcut_calc_c_u_quick_(
					uxc0,
					ux0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Uc,
					sz, g);
		} else if( g_pFFVConfig->ConvectionTermScheme == "Blend" ) {
			real alpha = 0.95;
			bcut_calc_c_f_blend_(
					uxc0,
					ux0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Uc,
					&alpha,
					sz, g);
		} else {
			bcut_calc_c_f_c2_(
					uxc0,
					ux0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Uc,
					sz, g);
		}

		bcut_calc_d_u_(
				uxd0,
				ux0,
				pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
				pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
				pPhaseId,
				&rhof,
				&mu,
				&dx, &dt,
				&Uc,
				sz, g);

		if( step==0 ) {
			copy_(
					uxcp,
					uxc0,
					sz, g);
			copy_(
					uxdp,
					uxd0,
					sz, g);
		}

		int axis=0;
		real gx = g_pFFVConfig->GravityX;
		real gy = g_pFFVConfig->GravityY;
		real gz = g_pFFVConfig->GravityZ;

		bcut_update_u_(
				ux0,
				uxc0, uxcp,
				uxd0, uxdp,
				p0,
				pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
				pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
				pPhaseId,
				&axis,
				&rhof,
				&mu,
				&dx, &dt,
				&Uc,
				&gx, &gy, &gz,
				sz, g);

		copy_(
				uxcp,
				uxc0,
				sz, g);
		copy_(
				uxdp,
				uxd0,
				sz, g);
	}
PM_Stop(tm_UpdateUX01);
/////////////////////////////////////////////

/////////////////////////////////////////////
// Impose B.C. (x)
/////////////////////////////////////////////
PM_Start(tm_UpdateUX04, 0, 0, true);
	plsUX0->ImposeBoundaryCondition(blockManager);
PM_Stop(tm_UpdateUX04);
/////////////////////////////////////////////
}

void Solver::UpdateUYe(int step) {
/////////////////////////////////////////////
// Calc A & b
/////////////////////////////////////////////
PM_Start(tm_UpdateUY01, 0, 0, true);
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
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
	
		real* vw = plsVw->GetBlockData(block);
		real* ve = plsVe->GetBlockData(block);
		real* vs = plsVs->GetBlockData(block);
		real* vn = plsVn->GetBlockData(block);
		real* vb = plsVb->GetBlockData(block);
		real* vt = plsVt->GetBlockData(block);

		real* uy0  = plsUY0->GetBlockData(block);
		real* uyc0 = plsUYC->GetBlockData(block);
		real* uycp = plsUYCP->GetBlockData(block);
		real* uyd0 = plsUYD->GetBlockData(block);
		real* uydp = plsUYDP->GetBlockData(block);
		real* p0   = plsP0->GetBlockData(block);

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

		int* pPhaseId = plsPhaseId->GetBlockData(block);

		real Uc = 0.0;

		real* ux0  = plsUX0->GetBlockData(block);
		real* uz0  = plsUZ0->GetBlockData(block);

		if( g_pFFVConfig->ConvectionTermScheme == "W3" ) {
			bcut_calc_c_f_w3_(
					uyc0,
					uy0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Uc,
					sz, g);
		} else if( g_pFFVConfig->ConvectionTermScheme == "E3" ) {
			bcut_calc_c_f_e3_(
					uyc0,
					uy0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Uc,
					sz, g);
		} else if( g_pFFVConfig->ConvectionTermScheme == "QUICK" ) {
			bcut_calc_c_f_quick_(
					uyc0,
					uy0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Uc,
					sz, g);
		} else if( g_pFFVConfig->ConvectionTermScheme == "QUICK2" ) {
			bcut_calc_c_u_quick_(
					uyc0,
					uy0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Uc,
					sz, g);
		} else if( g_pFFVConfig->ConvectionTermScheme == "Blend" ) {
			real alpha = 0.95;
			bcut_calc_c_f_blend_(
					uyc0,
					uy0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Uc,
					&alpha,
					sz, g);
		} else {
			bcut_calc_c_f_c2_(
					uyc0,
					uy0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Uc,
					sz, g);
		}

		bcut_calc_d_u_(
				uyd0,
				uy0,
				pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
				pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
				pPhaseId,
				&rhof,
				&mu,
				&dx, &dt,
				&Uc,
				sz, g);

		if( step==0 ) {
			copy_(
					uycp,
					uyc0,
					sz, g);
			copy_(
					uydp,
					uyd0,
					sz, g);
		}

		int axis=1;
		real gx = g_pFFVConfig->GravityX;
		real gy = g_pFFVConfig->GravityY;
		real gz = g_pFFVConfig->GravityZ;
		bcut_update_u_(
				uy0,
				uyc0, uycp,
				uyd0, uydp,
				p0,
				pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
				pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
				pPhaseId,
				&axis,
				&rhof,
				&mu,
				&dx, &dt,
				&Uc,
				&gx, &gy, &gz,
				sz, g);

		copy_(
				uycp,
				uyc0,
				sz, g);
		copy_(
				uydp,
				uyd0,
				sz, g);
	}
PM_Stop(tm_UpdateUY01);
/////////////////////////////////////////////

/////////////////////////////////////////////
// Impose B.C. (x)
/////////////////////////////////////////////
PM_Start(tm_UpdateUY04, 0, 0, true);
	plsUY0->ImposeBoundaryCondition(blockManager);
PM_Stop(tm_UpdateUY04);
/////////////////////////////////////////////
}

void Solver::UpdateUZe(int step) {
/////////////////////////////////////////////
// Calc A & b
/////////////////////////////////////////////
PM_Start(tm_UpdateUZ01, 0, 0, true);
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
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
	
		real* vw = plsVw->GetBlockData(block);
		real* ve = plsVe->GetBlockData(block);
		real* vs = plsVs->GetBlockData(block);
		real* vn = plsVn->GetBlockData(block);
		real* vb = plsVb->GetBlockData(block);
		real* vt = plsVt->GetBlockData(block);

		real* uz0  = plsUZ0->GetBlockData(block);
		real* uzc0 = plsUZC->GetBlockData(block);
		real* uzcp = plsUZCP->GetBlockData(block);
		real* uzd0 = plsUZD->GetBlockData(block);
		real* uzdp = plsUZDP->GetBlockData(block);
		real* p0   = plsP0->GetBlockData(block);

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

		int* pPhaseId = plsPhaseId->GetBlockData(block);

		real Uc = 0.0;

		real* ux0  = plsUX0->GetBlockData(block);
		real* uy0  = plsUY0->GetBlockData(block);

		if( g_pFFVConfig->ConvectionTermScheme == "W3" ) {
			bcut_calc_c_f_w3_(
					uzc0,
					uz0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Uc,
					sz, g);
		} else if( g_pFFVConfig->ConvectionTermScheme == "E3" ) {
			bcut_calc_c_f_e3_(
					uzc0,
					uz0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Uc,
					sz, g);
		} else if( g_pFFVConfig->ConvectionTermScheme == "QUICK" ) {
			bcut_calc_c_f_quick_(
					uzc0,
					uz0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Uc,
					sz, g);
		} else if( g_pFFVConfig->ConvectionTermScheme == "QUICK2" ) {
			bcut_calc_c_u_quick_(
					uzc0,
					uz0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Uc,
					sz, g);
		} else if( g_pFFVConfig->ConvectionTermScheme == "Blend" ) {
			real alpha = 0.95;
			bcut_calc_c_f_blend_(
					uzc0,
					uz0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Uc,
					&alpha,
					sz, g);
		} else {
			bcut_calc_c_f_c2_(
					uzc0,
					uz0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Uc,
					sz, g);
		}

		bcut_calc_d_u_(
				uzd0,
				uz0,
				pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
				pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
				pPhaseId,
				&rhof,
				&mu,
				&dx, &dt,
				&Uc,
				sz, g);

		if( step==0 ) {
			copy_(
					uzcp,
					uzc0,
					sz, g);
			copy_(
					uzdp,
					uzd0,
					sz, g);
		}

		int axis=2;
		real gx = g_pFFVConfig->GravityX;
		real gy = g_pFFVConfig->GravityY;
		real gz = g_pFFVConfig->GravityZ;
		bcut_update_u_(
				uz0,
				uzc0, uzcp,
				uzd0, uzdp,
				p0,
				pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
				pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
				pPhaseId,
				&axis,
				&rhof,
				&mu,
				&dx, &dt,
				&Uc,
				&gx, &gy, &gz,
				sz, g);

		copy_(
				uzcp,
				uzc0,
				sz, g);
		copy_(
				uzdp,
				uzd0,
				sz, g);
	}
PM_Stop(tm_UpdateUZ01);
/////////////////////////////////////////////

/////////////////////////////////////////////
// Impose B.C. (x)
/////////////////////////////////////////////
PM_Start(tm_UpdateUZ04, 0, 0, true);
	plsUZ0->ImposeBoundaryCondition(blockManager);
PM_Stop(tm_UpdateUZ04);
/////////////////////////////////////////////
}

void Solver::UpdateTe(int step) {
/////////////////////////////////////////////
// Calc A & b
/////////////////////////////////////////////
PM_Start(tm_UpdateT01, 0, 0, true);
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
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
		
		real* t0  = plsT0->GetBlockData(block);
		real* tc0 = plsTC->GetBlockData(block);
		real* tcp = plsTCP->GetBlockData(block);
		real* td0 = plsTD->GetBlockData(block);
		real* tdp = plsTDP->GetBlockData(block);

		real* vw = plsVw->GetBlockData(block);
		real* ve = plsVe->GetBlockData(block);
		real* vs = plsVs->GetBlockData(block);
		real* vn = plsVn->GetBlockData(block);
		real* vb = plsVb->GetBlockData(block);
		real* vt = plsVt->GetBlockData(block);

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

		int* pPhaseId = plsPhaseId->GetBlockData(block);

		int* pNormalIndex0 = plsNormalIndex0->GetBlockData(block);
		int* pNormalIndex1 = plsNormalIndex1->GetBlockData(block);
		int* pNormalIndex2 = plsNormalIndex2->GetBlockData(block);
		int* pNormalIndex3 = plsNormalIndex3->GetBlockData(block);
		int* pNormalIndex4 = plsNormalIndex4->GetBlockData(block);
		int* pNormalIndex5 = plsNormalIndex5->GetBlockData(block);

		int bc_n[1] = {32};
		int bc_type[32];
		real bc_value[32];
		for(int m=0; m<bc_n[0]; m++) {
			bc_type[m]  = g_pFFVConfig->BCInternalBoundaryType[m];
			bc_value[m] = g_pFFVConfig->BCInternalBoundaryValue[m];
//			std::cout << bc_type[n] << " " << bc_value[n] << std::endl;
		}

		real Tc = 1.0;

		if( g_pFFVConfig->ConvectionTermScheme == "W3" ) {
			bcut_calc_c_f_w3_(
					tc0,
					t0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Tc,
					sz, g);
		} else if( g_pFFVConfig->ConvectionTermScheme == "E3" ) {
			bcut_calc_c_f_e3_(
					tc0,
					t0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Tc,
					sz, g);
		} else if( g_pFFVConfig->ConvectionTermScheme == "QUICK" ) {
			bcut_calc_c_f_quick_(
					tc0,
					t0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Tc,
					sz, g);
		} else if( g_pFFVConfig->ConvectionTermScheme == "QUICK2" ) {
			bcut_calc_c_f_quick_(
					tc0,
					t0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Tc,
					sz, g);
		} else if( g_pFFVConfig->ConvectionTermScheme == "Blend" ) {
			real alpha = 0.95;
			bcut_calc_c_f_blend_(
					tc0,
					t0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Tc,
					&alpha,
					sz, g);
		} else {
			bcut_calc_c_f_c2_(
					tc0,
					t0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Tc,
					sz, g);
		}

		bcut_calc_d_t_(
				td0,
				t0,
				pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
				pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
				pPhaseId,
				&pNormalN[n],
				pNormalX[n],
				pNormalY[n],
				pNormalZ[n],
				pNormalIndex0, pNormalIndex1, pNormalIndex2, pNormalIndex3, pNormalIndex4, pNormalIndex5,
				&rhof, &rhos,
				&cpf, &cps,
				&kf, &ks,
				bc_n,
				bc_type,
				bc_value,
				org,
				&dx, &dt,
				sz, g);

		if( step==0 ) {
			copy_(
					tcp,
					tc0,
					sz, g);
			copy_(
					tdp,
					td0,
					sz, g);
		}

		bcut_update_t_(
				t0,
				tc0, tcp,
				td0, tdp,
				pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
				pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
				pPhaseId,
				&rhof, &rhos,
				&cpf, &cps,
				&kf, &ks,
				&dx, &dt,
				sz, g);

		copy_(
				tcp,
				tc0,
				sz, g);
		copy_(
				tdp,
				td0,
				sz, g);
	}
PM_Stop(tm_UpdateT01);
/////////////////////////////////////////////

/////////////////////////////////////////////
// Impose B.C. (x)
/////////////////////////////////////////////
PM_Start(tm_UpdateT04, 0, 0, true);
	plsT0->ImposeBoundaryCondition(blockManager);
PM_Stop(tm_UpdateT04);
/////////////////////////////////////////////
}

void Solver::UpdateUX(int step) {
/////////////////////////////////////////////
// Calc A & b
/////////////////////////////////////////////
PM_Start(tm_UpdateUX01, 0, 0, true);
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
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
	
		real* Ap = plsAp->GetBlockData(block);
		real* Aw = plsAw->GetBlockData(block);
		real* Ae = plsAe->GetBlockData(block);
		real* As = plsAs->GetBlockData(block);
		real* An = plsAn->GetBlockData(block);
		real* Ab = plsAb->GetBlockData(block);
		real* At = plsAt->GetBlockData(block);
		real* b  = plsb ->GetBlockData(block);

		real* vw = plsVw->GetBlockData(block);
		real* ve = plsVe->GetBlockData(block);
		real* vs = plsVs->GetBlockData(block);
		real* vn = plsVn->GetBlockData(block);
		real* vb = plsVb->GetBlockData(block);
		real* vt = plsVt->GetBlockData(block);

		real* ux0  = plsUX0->GetBlockData(block);
		real* uxc0 = plsUXC->GetBlockData(block);
		real* uxcp = plsUXCP->GetBlockData(block);
		real* uxd0 = plsUXD->GetBlockData(block);
		real* p0   = plsP0->GetBlockData(block);

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

		int* pPhaseId = plsPhaseId->GetBlockData(block);

		real Uc = 0.0;

		real* uy0  = plsUY0->GetBlockData(block);
		real* uz0  = plsUZ0->GetBlockData(block);

		if( g_pFFVConfig->ConvectionTermScheme == "W3" ) {
			bcut_calc_c_f_w3_(
					uxc0,
					ux0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Uc,
					sz, g);
		} else if( g_pFFVConfig->ConvectionTermScheme == "E3" ) {
			bcut_calc_c_f_e3_(
					uxc0,
					ux0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Uc,
					sz, g);
		} else if( g_pFFVConfig->ConvectionTermScheme == "QUICK" ) {
			bcut_calc_c_f_quick_(
					uxc0,
					ux0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Uc,
					sz, g);
		} else if( g_pFFVConfig->ConvectionTermScheme == "QUICK2" ) {
			bcut_calc_c_u_quick_(
					uxc0,
					ux0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Uc,
					sz, g);
		} else if( g_pFFVConfig->ConvectionTermScheme == "Blend" ) {
			real alpha = 0.95;
			bcut_calc_c_f_blend_(
					uxc0,
					ux0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Uc,
					&alpha,
					sz, g);
		} else {
			bcut_calc_c_f_c2_(
					uxc0,
					ux0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Uc,
					sz, g);
		}

		bcut_calc_d_u_(
				uxd0,
				ux0,
				pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
				pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
				pPhaseId,
				&rhof,
				&mu,
				&dx, &dt,
				&Uc,
				sz, g);

		int axis=0;
		real gx = g_pFFVConfig->GravityX;
		real gy = g_pFFVConfig->GravityY;
		real gz = g_pFFVConfig->GravityZ;
		bcut_calc_ab_u_(
				Ap, Aw, Ae, As, An, Ab, At, b,
				ux0,
				uxc0, uxcp,
				uxd0,
				p0,
				pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
				pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
				pPhaseId,
				&axis,
				&rhof,
				&mu,
				&dx, &dt,
				&Uc,
				&gx, &gy, &gz,
				sz, g);
		copy_(
				uxcp,
				uxc0,
				sz, g);
	}
PM_Stop(tm_UpdateUX01);
/////////////////////////////////////////////

/////////////////////////////////////////////
// Impose B.C. (A)
/////////////////////////////////////////////
PM_Start(tm_UpdateUX02, 0, 0, true);
	plsUX0->ImposeBoundaryCondition(blockManager, plsAp, plsAw, plsAe, plsAs, plsAn, plsAb, plsAt, plsb);
PM_Stop(tm_UpdateUX02);
/////////////////////////////////////////////

/////////////////////////////////////////////
// Solve Ax = b
/////////////////////////////////////////////
PM_Start(tm_UpdateUX03, 0, 0, true);
	if( g_pFFVConfig->TuningMasking == true ) {
		pils->BiCGSTAB_Mask(
							blockManager,
							plsUX0,
							plsAp,
							plsAw,
							plsAe,
							plsAs,
							plsAn,
							plsAb,
							plsAt,
							plsb,
							plsr,
							plsr0,
							plsp,
							plsp_,
							plsq_,
							plss,
							plss_,
							plst_,
							plsUX1,
							plsMaskId,
							omegaU,
							countPreConditionerU,
							countMaxU,
							epsilonU,
							countUX,
							residualUX);
	} else {
		pils->BiCGSTAB(
							blockManager,
							plsUX0,
							plsAp,
							plsAw,
							plsAe,
							plsAs,
							plsAn,
							plsAb,
							plsAt,
							plsb,
							plsr,
							plsr0,
							plsp,
							plsp_,
							plsq_,
							plss,
							plss_,
							plst_,
							plsUX1,
							omegaU,
							countPreConditionerU,
							countMaxU,
							epsilonU,
							countUX,
							residualUX);
	}
PM_Stop(tm_UpdateUX03);
/////////////////////////////////////////////

/////////////////////////////////////////////
// Impose B.C. (x)
/////////////////////////////////////////////
PM_Start(tm_UpdateUX04, 0, 0, true);
	plsUX0->ImposeBoundaryCondition(blockManager);
PM_Stop(tm_UpdateUX04);
/////////////////////////////////////////////
}

void Solver::UpdateUY(int step) {
/////////////////////////////////////////////
// Calc A & b
/////////////////////////////////////////////
PM_Start(tm_UpdateUY01, 0, 0, true);
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
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
	
		real* Ap = plsAp->GetBlockData(block);
		real* Aw = plsAw->GetBlockData(block);
		real* Ae = plsAe->GetBlockData(block);
		real* As = plsAs->GetBlockData(block);
		real* An = plsAn->GetBlockData(block);
		real* Ab = plsAb->GetBlockData(block);
		real* At = plsAt->GetBlockData(block);
		real* b  = plsb ->GetBlockData(block);

		real* vw = plsVw->GetBlockData(block);
		real* ve = plsVe->GetBlockData(block);
		real* vs = plsVs->GetBlockData(block);
		real* vn = plsVn->GetBlockData(block);
		real* vb = plsVb->GetBlockData(block);
		real* vt = plsVt->GetBlockData(block);

		real* uy0  = plsUY0->GetBlockData(block);
		real* uyc0 = plsUYC->GetBlockData(block);
		real* uycp = plsUYCP->GetBlockData(block);
		real* uyd0 = plsUYD->GetBlockData(block);
		real* p0   = plsP0->GetBlockData(block);

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

		int* pPhaseId = plsPhaseId->GetBlockData(block);

		real Uc = 0.0;

		real* ux0  = plsUX0->GetBlockData(block);
		real* uz0  = plsUZ0->GetBlockData(block);

		if( g_pFFVConfig->ConvectionTermScheme == "W3" ) {
			bcut_calc_c_f_w3_(
					uyc0,
					uy0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Uc,
					sz, g);
		} else if( g_pFFVConfig->ConvectionTermScheme == "E3" ) {
			bcut_calc_c_f_e3_(
					uyc0,
					uy0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Uc,
					sz, g);
		} else if( g_pFFVConfig->ConvectionTermScheme == "QUICK" ) {
			bcut_calc_c_f_quick_(
					uyc0,
					uy0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Uc,
					sz, g);
		} else if( g_pFFVConfig->ConvectionTermScheme == "QUICK2" ) {
			bcut_calc_c_u_quick_(
					uyc0,
					uy0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Uc,
					sz, g);
		} else if( g_pFFVConfig->ConvectionTermScheme == "Blend" ) {
			real alpha = 0.95;
			bcut_calc_c_f_blend_(
					uyc0,
					uy0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Uc,
					&alpha,
					sz, g);
		} else {
			bcut_calc_c_f_c2_(
					uyc0,
					uy0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Uc,
					sz, g);
		}

		bcut_calc_d_u_(
				uyd0,
				uy0,
				pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
				pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
				pPhaseId,
				&rhof,
				&mu,
				&dx, &dt,
				&Uc,
				sz, g);

		int axis=1;
		real gx = g_pFFVConfig->GravityX;
		real gy = g_pFFVConfig->GravityY;
		real gz = g_pFFVConfig->GravityZ;
		bcut_calc_ab_u_(
				Ap, Aw, Ae, As, An, Ab, At, b,
				uy0,
				uyc0, uycp,
				uyd0,
				p0,
				pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
				pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
				pPhaseId,
				&axis,
				&rhof,
				&mu,
				&dx, &dt,
				&Uc,
				&gx, &gy, &gz,
				sz, g);
		copy_(
				uycp,
				uyc0,
				sz, g);
	}
PM_Stop(tm_UpdateUY01);
/////////////////////////////////////////////

/////////////////////////////////////////////
// Impose B.C. (A)
/////////////////////////////////////////////
PM_Start(tm_UpdateUY02, 0, 0, true);
	plsUY0->ImposeBoundaryCondition(blockManager, plsAp, plsAw, plsAe, plsAs, plsAn, plsAb, plsAt, plsb);
PM_Stop(tm_UpdateUY02);
/////////////////////////////////////////////

/////////////////////////////////////////////
// Solve Ax = b
/////////////////////////////////////////////
PM_Start(tm_UpdateUY03, 0, 0, true);
	if( g_pFFVConfig->TuningMasking == true ) {
		pils->BiCGSTAB_Mask(
							blockManager,
							plsUY0,
							plsAp,
							plsAw,
							plsAe,
							plsAs,
							plsAn,
							plsAb,
							plsAt,
							plsb,
							plsr,
							plsr0,
							plsp,
							plsp_,
							plsq_,
							plss,
							plss_,
							plst_,
							plsUY1,
							plsMaskId,
							omegaU,
							countPreConditionerU,
							countMaxU,
							epsilonU,
							countUY,
							residualUY);
	} else {
		pils->BiCGSTAB(
							blockManager,
							plsUY0,
							plsAp,
							plsAw,
							plsAe,
							plsAs,
							plsAn,
							plsAb,
							plsAt,
							plsb,
							plsr,
							plsr0,
							plsp,
							plsp_,
							plsq_,
							plss,
							plss_,
							plst_,
							plsUY1,
							omegaU,
							countPreConditionerU,
							countMaxU,
							epsilonU,
							countUY,
							residualUY);
	}
PM_Stop(tm_UpdateUY03);
/////////////////////////////////////////////

/////////////////////////////////////////////
// Impose B.C. (x)
/////////////////////////////////////////////
PM_Start(tm_UpdateUY04, 0, 0, true);
	plsUY0->ImposeBoundaryCondition(blockManager);
PM_Stop(tm_UpdateUY04);
/////////////////////////////////////////////
}

void Solver::UpdateUZ(int step) {
/////////////////////////////////////////////
// Calc A & b
/////////////////////////////////////////////
PM_Start(tm_UpdateUZ01, 0, 0, true);
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
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
	
		real* Ap = plsAp->GetBlockData(block);
		real* Aw = plsAw->GetBlockData(block);
		real* Ae = plsAe->GetBlockData(block);
		real* As = plsAs->GetBlockData(block);
		real* An = plsAn->GetBlockData(block);
		real* Ab = plsAb->GetBlockData(block);
		real* At = plsAt->GetBlockData(block);
		real* b  = plsb ->GetBlockData(block);

		real* vw = plsVw->GetBlockData(block);
		real* ve = plsVe->GetBlockData(block);
		real* vs = plsVs->GetBlockData(block);
		real* vn = plsVn->GetBlockData(block);
		real* vb = plsVb->GetBlockData(block);
		real* vt = plsVt->GetBlockData(block);

		real* uz0  = plsUZ0->GetBlockData(block);
		real* uzc0 = plsUZC->GetBlockData(block);
		real* uzcp = plsUZCP->GetBlockData(block);
		real* uzd0 = plsUZD->GetBlockData(block);
		real* p0   = plsP0->GetBlockData(block);

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

		int* pPhaseId = plsPhaseId->GetBlockData(block);

		real Uc = 0.0;

		real* ux0  = plsUX0->GetBlockData(block);
		real* uy0  = plsUY0->GetBlockData(block);

		if( g_pFFVConfig->ConvectionTermScheme == "W3" ) {
			bcut_calc_c_f_w3_(
					uzc0,
					uz0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Uc,
					sz, g);
		} else if( g_pFFVConfig->ConvectionTermScheme == "E3" ) {
			bcut_calc_c_f_e3_(
					uzc0,
					uz0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Uc,
					sz, g);
		} else if( g_pFFVConfig->ConvectionTermScheme == "QUICK" ) {
			bcut_calc_c_f_quick_(
					uzc0,
					uz0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Uc,
					sz, g);
		} else if( g_pFFVConfig->ConvectionTermScheme == "QUICK2" ) {
			bcut_calc_c_u_quick_(
					uzc0,
					uz0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Uc,
					sz, g);
		} else if( g_pFFVConfig->ConvectionTermScheme == "Blend" ) {
			real alpha = 0.95;
			bcut_calc_c_f_blend_(
					uzc0,
					uz0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Uc,
					&alpha,
					sz, g);
		} else {
			bcut_calc_c_f_c2_(
					uzc0,
					uz0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Uc,
					sz, g);
		}

		bcut_calc_d_u_(
				uzd0,
				uz0,
				pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
				pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
				pPhaseId,
				&rhof,
				&mu,
				&dx, &dt,
				&Uc,
				sz, g);

		int axis=2;
		real gx = g_pFFVConfig->GravityX;
		real gy = g_pFFVConfig->GravityY;
		real gz = g_pFFVConfig->GravityZ;
		bcut_calc_ab_u_(
				Ap, Aw, Ae, As, An, Ab, At, b,
				uz0,
				uzc0, uzcp,
				uzd0,
				p0,
				pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
				pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
				pPhaseId,
				&axis,
				&rhof,
				&mu,
				&dx, &dt,
				&Uc,
				&gx, &gy, &gz,
				sz, g);
		copy_(
				uzcp,
				uzc0,
				sz, g);
	}
PM_Stop(tm_UpdateUZ01);
/////////////////////////////////////////////

/////////////////////////////////////////////
// Impose B.C. (A)
/////////////////////////////////////////////
PM_Start(tm_UpdateUZ02, 0, 0, true);
	plsUZ0->ImposeBoundaryCondition(blockManager, plsAp, plsAw, plsAe, plsAs, plsAn, plsAb, plsAt, plsb);
PM_Stop(tm_UpdateUZ02);
/////////////////////////////////////////////

/////////////////////////////////////////////
// Solve Ax = b
/////////////////////////////////////////////
PM_Start(tm_UpdateUZ03, 0, 0, true);
	if( g_pFFVConfig->TuningMasking == true ) {
		pils->BiCGSTAB_Mask(
							blockManager,
							plsUZ0,
							plsAp,
							plsAw,
							plsAe,
							plsAs,
							plsAn,
							plsAb,
							plsAt,
							plsb,
							plsr,
							plsr0,
							plsp,
							plsp_,
							plsq_,
							plss,
							plss_,
							plst_,
							plsUZ1,
							plsMaskId,
							omegaU,
							countPreConditionerU,
							countMaxU,
							epsilonU,
							countUZ,
							residualUZ);
	} else {
		pils->BiCGSTAB(
							blockManager,
							plsUZ0,
							plsAp,
							plsAw,
							plsAe,
							plsAs,
							plsAn,
							plsAb,
							plsAt,
							plsb,
							plsr,
							plsr0,
							plsp,
							plsp_,
							plsq_,
							plss,
							plss_,
							plst_,
							plsUZ1,
							omegaU,
							countPreConditionerU,
							countMaxU,
							epsilonU,
							countUZ,
							residualUZ);
	}
PM_Stop(tm_UpdateUZ03);
/////////////////////////////////////////////

/////////////////////////////////////////////
// Impose B.C. (x)
/////////////////////////////////////////////
PM_Start(tm_UpdateUZ04, 0, 0, true);
	plsUZ0->ImposeBoundaryCondition(blockManager);
PM_Stop(tm_UpdateUZ04);
/////////////////////////////////////////////
}

void Solver::UpdateP(int step) {
/////////////////////////////////////////////
// Remove Grad. P
/////////////////////////////////////////////
PM_Start(tm_UpdateP01, 0, 0, true);
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
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
	
		real* Ap = plsAp->GetBlockData(block);
		real* Aw = plsAw->GetBlockData(block);
		real* Ae = plsAe->GetBlockData(block);
		real* As = plsAs->GetBlockData(block);
		real* An = plsAn->GetBlockData(block);
		real* Ab = plsAb->GetBlockData(block);
		real* At = plsAt->GetBlockData(block);
		real* b  = plsb ->GetBlockData(block);

		real* vw = plsVw->GetBlockData(block);
		real* ve = plsVe->GetBlockData(block);
		real* vs = plsVs->GetBlockData(block);
		real* vn = plsVn->GetBlockData(block);
		real* vb = plsVb->GetBlockData(block);
		real* vt = plsVt->GetBlockData(block);

		real* p0 = plsP0->GetBlockData(block);
		real* t0 = plsT0->GetBlockData(block);

		real* ux = plsUX0->GetBlockData(block);
		real* uy = plsUY0->GetBlockData(block);
		real* uz = plsUZ0->GetBlockData(block);

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

		int* pPhaseId = plsPhaseId->GetBlockData(block);

		real gx = g_pFFVConfig->GravityX;
		real gy = g_pFFVConfig->GravityY;
		real gz = g_pFFVConfig->GravityZ;
		real betag = g_pFFVConfig->BetaG;
		real tr = g_pFFVConfig->Tref;

		bcut_add_g_(
				ux, uy, uz,
				t0,
				&gx, &gy, &gz,
				&betag, &tr,
				&dx, &dt,
				sz, g);

		bcut_remove_p_(
				ux, uy, uz,
				p0,
				pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
				pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
				pPhaseId,
				&rhof,
				&dx, &dt,
				sz, g);
	}
PM_Stop(tm_UpdateP01);
/////////////////////////////////////////////

/////////////////////////////////////////////
// Impose B.C UX, UY, UZ
/////////////////////////////////////////////
PM_Start(tm_UpdateP02, 0, 0, true);
	plsUX0->ImposeBoundaryCondition(blockManager);
	plsUY0->ImposeBoundaryCondition(blockManager);
	plsUZ0->ImposeBoundaryCondition(blockManager);
PM_Stop(tm_UpdateP02);
/////////////////////////////////////////////

/////////////////////////////////////////////
// Calc A & b
/////////////////////////////////////////////
PM_Start(tm_UpdateP03, 0, 0, true);
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
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
	
		real* Ap = plsAp->GetBlockData(block);
		real* Aw = plsAw->GetBlockData(block);
		real* Ae = plsAe->GetBlockData(block);
		real* As = plsAs->GetBlockData(block);
		real* An = plsAn->GetBlockData(block);
		real* Ab = plsAb->GetBlockData(block);
		real* At = plsAt->GetBlockData(block);
		real* b  = plsb ->GetBlockData(block);

		real* vw = plsVw->GetBlockData(block);
		real* ve = plsVe->GetBlockData(block);
		real* vs = plsVs->GetBlockData(block);
		real* vn = plsVn->GetBlockData(block);
		real* vb = plsVb->GetBlockData(block);
		real* vt = plsVt->GetBlockData(block);

		real* p0 = plsP0->GetBlockData(block);

		real* ux = plsUX0->GetBlockData(block);
		real* uy = plsUY0->GetBlockData(block);
		real* uz = plsUZ0->GetBlockData(block);

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

		int* pPhaseId = plsPhaseId->GetBlockData(block);

		bcut_calc_ab_p_(
				Ap, Aw, Ae, As, An, Ab, At, b,
				vw, ve, vs, vn, vb, vt,
				p0,
				ux, uy, uz,
				pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
				pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
				pPhaseId,
				&rhof,
				&dx, &dt,
				sz, g);
	}
PM_Stop(tm_UpdateP03);
/////////////////////////////////////////////

/////////////////////////////////////////////
// Set ref. P
/////////////////////////////////////////////
PM_Start(tm_UpdateP04, 0, 0, true);
	if( g_pFFVConfig->IterationReferencePressureActive ) {
		real xr = g_pFFVConfig->IterationReferencePressurePoint.x;
		real yr = g_pFFVConfig->IterationReferencePressurePoint.y;
		real zr = g_pFFVConfig->IterationReferencePressurePoint.z;
		real pr = g_pFFVConfig->IterationReferencePressureValue;
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
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
		
			real* Ap = plsAp->GetBlockData(block);
			real* Aw = plsAw->GetBlockData(block);
			real* Ae = plsAe->GetBlockData(block);
			real* As = plsAs->GetBlockData(block);
			real* An = plsAn->GetBlockData(block);
			real* Ab = plsAb->GetBlockData(block);
			real* At = plsAt->GetBlockData(block);
			real* b  = plsb ->GetBlockData(block);

			real* vw = plsVw->GetBlockData(block);
			real* ve = plsVe->GetBlockData(block);
			real* vs = plsVs->GetBlockData(block);
			real* vn = plsVn->GetBlockData(block);
			real* vb = plsVb->GetBlockData(block);
			real* vt = plsVt->GetBlockData(block);

			real* p0 = plsP0->GetBlockData(block);

			real* ux = plsUX0->GetBlockData(block);
			real* uy = plsUY0->GetBlockData(block);
			real* uz = plsUZ0->GetBlockData(block);

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

			int* pPhaseId = plsPhaseId->GetBlockData(block);

			bcut_set_reference_value_(
					Ap, Aw, Ae, As, An, Ab, At, b,
					&xr, &yr, &zr,
					&pr,
					&dx,
					org,
					sz, g);
		}
	}
PM_Stop(tm_UpdateP04);
/////////////////////////////////////////////

/////////////////////////////////////////////
// Impose B.C. (A)
/////////////////////////////////////////////
PM_Start(tm_UpdateP05, 0, 0, true);
	plsP0->ImposeBoundaryCondition(blockManager, plsAp, plsAw, plsAe, plsAs, plsAn, plsAb, plsAt, plsb);
PM_Stop(tm_UpdateP05);
/////////////////////////////////////////////

/////////////////////////////////////////////
// Solve Ax = b
/////////////////////////////////////////////
PM_Start(tm_UpdateP06, 0, 0, true);
	if( g_pFFVConfig->IterationSolverP == "RBGS" ) {
		pils->RBGS(
							blockManager,
							plsP0,
							plsAp,
							plsAw,
							plsAe,
							plsAs,
							plsAn,
							plsAb,
							plsAt,
							plsb,
							omegaP,
							countMaxP,
							epsilonP,
							countP,
							residualP);
	} else {
		if( g_pFFVConfig->TuningMasking == true ) {
			pils->BiCGSTAB_Mask(
								blockManager,
								plsP0,
								plsAp,
								plsAw,
								plsAe,
								plsAs,
								plsAn,
								plsAb,
								plsAt,
								plsb,
								plsr,
								plsr0,
								plsp,
								plsp_,
								plsq_,
								plss,
								plss_,
								plst_,
								plsP1,
								plsMaskId,
								omegaP,
								countPreConditionerP,
								countMaxP,
								epsilonP,
								countP,
								residualP);
		} else {
			pils->BiCGSTAB(
								blockManager,
								plsP0,
								plsAp,
								plsAw,
								plsAe,
								plsAs,
								plsAn,
								plsAb,
								plsAt,
								plsb,
								plsr,
								plsr0,
								plsp,
								plsp_,
								plsq_,
								plss,
								plss_,
								plst_,
								plsP1,
								omegaP,
								countPreConditionerP,
								countMaxP,
								epsilonP,
								countP,
								residualP);
		}
	}
PM_Stop(tm_UpdateP06);
/////////////////////////////////////////////

/////////////////////////////////////////////
// Impose B.C. (x)
/////////////////////////////////////////////
PM_Start(tm_UpdateP07, 0, 0, true);
	plsP0->ImposeBoundaryCondition(blockManager);
PM_Stop(tm_UpdateP07);
/////////////////////////////////////////////
}

void Solver::UpdateU(int step) {
/////////////////////////////////////////////
// Correct U
/////////////////////////////////////////////
PM_Start(tm_UpdateU01, 0, 0, true);
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
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

		real* ux0 = plsUX0->GetBlockData(block);
		real* uy0 = plsUY0->GetBlockData(block);
		real* uz0 = plsUZ0->GetBlockData(block);
		real* p0  = plsP0 ->GetBlockData(block);
		real* lapp= plsLapP->GetBlockData(block);

		real* vw = plsVw->GetBlockData(block);
		real* ve = plsVe->GetBlockData(block);
		real* vs = plsVs->GetBlockData(block);
		real* vn = plsVn->GetBlockData(block);
		real* vb = plsVb->GetBlockData(block);
		real* vt = plsVt->GetBlockData(block);

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

		int* pPhaseId = plsPhaseId->GetBlockData(block);

		bcut_corr_u_(
				ux0, uy0, uz0,
				vw, ve, vs, vn, vb, vt,
				lapp,
				p0,
				pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
				pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
				pPhaseId,
				&rhof,
				&dx, &dt,
				sz, g);
	}
PM_Stop(tm_UpdateU01);
/////////////////////////////////////////////

/////////////////////////////////////////////
// ImposeBC UX, UY, UZ
/////////////////////////////////////////////
PM_Start(tm_UpdateU02, 0, 0, true);
	plsUX0->ImposeBoundaryCondition(blockManager);
	plsUY0->ImposeBoundaryCondition(blockManager);
	plsUZ0->ImposeBoundaryCondition(blockManager);
PM_Stop(tm_UpdateU02);
/////////////////////////////////////////////

/////////////////////////////////////////////
// Calc Laplace P
/////////////////////////////////////////////
PM_Start(tm_UpdateU03, 0, 0, true);
	plsP0->ImposeBoundaryCondition(blockManager);
	plsLapP->ImposeBoundaryCondition(blockManager);
PM_Stop(tm_UpdateU03);
/////////////////////////////////////////////
}

void Solver::UpdateT(int step) {
/////////////////////////////////////////////
// Calc A & b
/////////////////////////////////////////////
PM_Start(tm_UpdateT01, 0, 0, true);
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
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
		
		real* Ap = plsAp->GetBlockData(block);
		real* Aw = plsAw->GetBlockData(block);
		real* Ae = plsAe->GetBlockData(block);
		real* As = plsAs->GetBlockData(block);
		real* An = plsAn->GetBlockData(block);
		real* Ab = plsAb->GetBlockData(block);
		real* At = plsAt->GetBlockData(block);
		real* b  = plsb ->GetBlockData(block);

		real* t0  = plsT0->GetBlockData(block);
		real* tc0 = plsTC->GetBlockData(block);
		real* tcp = plsTCP->GetBlockData(block);
		real* td0 = plsTD->GetBlockData(block);

		real* vw = plsVw->GetBlockData(block);
		real* ve = plsVe->GetBlockData(block);
		real* vs = plsVs->GetBlockData(block);
		real* vn = plsVn->GetBlockData(block);
		real* vb = plsVb->GetBlockData(block);
		real* vt = plsVt->GetBlockData(block);

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

		int* pPhaseId = plsPhaseId->GetBlockData(block);

		int* pNormalIndex0 = plsNormalIndex0->GetBlockData(block);
		int* pNormalIndex1 = plsNormalIndex1->GetBlockData(block);
		int* pNormalIndex2 = plsNormalIndex2->GetBlockData(block);
		int* pNormalIndex3 = plsNormalIndex3->GetBlockData(block);
		int* pNormalIndex4 = plsNormalIndex4->GetBlockData(block);
		int* pNormalIndex5 = plsNormalIndex5->GetBlockData(block);

		int bc_n[1] = {32};
		int bc_type[32];
		real bc_value[32];
		for(int m=0; m<bc_n[0]; m++) {
			bc_type[m]  = g_pFFVConfig->BCInternalBoundaryType[m];
			bc_value[m] = g_pFFVConfig->BCInternalBoundaryValue[m];
//			std::cout << m << " " << bc_type[m] << " " << bc_value[m] << std::endl;
//			std::cout << m << " " << BCInternalBoundaryType[n] << " " << BCInternalBoundaryValue[n] << std::endl;
		}

		real Tc = 1.0;
		if( g_pFFVConfig->ConvectionTermScheme == "W3" ) {
			bcut_calc_c_f_w3_(
					tc0,
					t0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Tc,
					sz, g);
		} else if( g_pFFVConfig->ConvectionTermScheme == "E3" ) {
			bcut_calc_c_f_e3_(
					tc0,
					t0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Tc,
					sz, g);
		} else if( g_pFFVConfig->ConvectionTermScheme == "QUICK" ) {
			bcut_calc_c_f_quick_(
					tc0,
					t0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Tc,
					sz, g);
		} else if( g_pFFVConfig->ConvectionTermScheme == "QUICK2" ) {
			bcut_calc_c_f_quick_(
					tc0,
					t0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Tc,
					sz, g);
		} else if( g_pFFVConfig->ConvectionTermScheme == "Blend" ) {
			real alpha = 0.95;
			bcut_calc_c_f_blend_(
					tc0,
					t0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Tc,
					&alpha,
					sz, g);
		} else {
			bcut_calc_c_f_c2_(
					tc0,
					t0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Tc,
					sz, g);
		}

		bcut_calc_abd_t_(
				Ap, Aw, Ae, As, An, Ab, At, b,
				t0,
				tc0, tcp,
				td0,
				pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
				pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
				pPhaseId,
				&pNormalN[n],
				pNormalX[n],
				pNormalY[n],
				pNormalZ[n],
				pNormalIndex0, pNormalIndex1, pNormalIndex2, pNormalIndex3, pNormalIndex4, pNormalIndex5,
				&rhof, &rhos,
				&cpf, &cps,
				&kf, &ks,
				bc_n,
				bc_type,
				bc_value,
				org,
				&dx, &dt,
				sz, g);
		copy_(
				tcp,
				tc0,
				sz, g);
	}
PM_Stop(tm_UpdateT01);
/////////////////////////////////////////////

/////////////////////////////////////////////
// Set ref. T
/////////////////////////////////////////////
PM_Start(tm_UpdateT02, 0, 0, true);
	if( g_pFFVConfig->IterationReferenceTemperatureActive ) {
		real xr = g_pFFVConfig->IterationReferenceTemperaturePoint.x;
		real yr = g_pFFVConfig->IterationReferenceTemperaturePoint.y;
		real zr = g_pFFVConfig->IterationReferenceTemperaturePoint.z;
		real tr = g_pFFVConfig->IterationReferenceTemperatureValue;
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
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
		
			real* Ap = plsAp->GetBlockData(block);
			real* Aw = plsAw->GetBlockData(block);
			real* Ae = plsAe->GetBlockData(block);
			real* As = plsAs->GetBlockData(block);
			real* An = plsAn->GetBlockData(block);
			real* Ab = plsAb->GetBlockData(block);
			real* At = plsAt->GetBlockData(block);
			real* b  = plsb ->GetBlockData(block);

			real* vw = plsVw->GetBlockData(block);
			real* ve = plsVe->GetBlockData(block);
			real* vs = plsVs->GetBlockData(block);
			real* vn = plsVn->GetBlockData(block);
			real* vb = plsVb->GetBlockData(block);
			real* vt = plsVt->GetBlockData(block);

			real* t0 = plsT0->GetBlockData(block);

			real* ux = plsUX0->GetBlockData(block);
			real* uy = plsUY0->GetBlockData(block);
			real* uz = plsUZ0->GetBlockData(block);

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

			int* pPhaseId = plsPhaseId->GetBlockData(block);

			bcut_set_reference_value_(
					Ap, Aw, Ae, As, An, Ab, At, b,
					&xr, &yr, &zr,
					&tr,
					&dx,
					org,
					sz, g);
		}
	}
PM_Stop(tm_UpdateT02);
/////////////////////////////////////////////

/////////////////////////////////////////////
// Impose B.C. (A)
/////////////////////////////////////////////
PM_Start(tm_UpdateT03, 0, 0, true);
	plsT0->ImposeBoundaryCondition(blockManager, plsAp, plsAw, plsAe, plsAs, plsAn, plsAb, plsAt, plsb);
PM_Stop(tm_UpdateT03);
/////////////////////////////////////////////

/////////////////////////////////////////////
// Solve Ax = b
/////////////////////////////////////////////
PM_Start(tm_UpdateT04, 0, 0, true);
	if( g_pFFVConfig->TuningMasking == true ) {
		pils->BiCGSTAB_Mask(
							blockManager,
							plsT0,
							plsAp,
							plsAw,
							plsAe,
							plsAs,
							plsAn,
							plsAb,
							plsAt,
							plsb,
							plsr,
							plsr0,
							plsp,
							plsp_,
							plsq_,
							plss,
							plss_,
							plst_,
							plsT1,
							plsMaskId,
							omegaT,
							countPreConditionerT,
							countMaxT,
							epsilonT,
							countT,
							residualT);
	} else {
		pils->BiCGSTAB(
							blockManager,
							plsT0,
							plsAp,
							plsAw,
							plsAe,
							plsAs,
							plsAn,
							plsAb,
							plsAt,
							plsb,
							plsr,
							plsr0,
							plsp,
							plsp_,
							plsq_,
							plss,
							plss_,
							plst_,
							plsT1,
							omegaT,
							countPreConditionerT,
							countMaxT,
							epsilonT,
							countT,
							residualT);
	}
PM_Stop(tm_UpdateT04);
/////////////////////////////////////////////

/////////////////////////////////////////////
// Impose B.C. (x)
/////////////////////////////////////////////
PM_Start(tm_UpdateT05, 0, 0, true);
	plsT0->ImposeBoundaryCondition(blockManager);
PM_Stop(tm_UpdateT05);
/////////////////////////////////////////////
}

#include <sys/time.h>
double Solver::GetTime() {
	struct timeval tp;
	int i = gettimeofday(&tp, 0);
	return ((double)(tp.tv_sec) + (double)(tp.tv_usec)*1.0e-6);
}

void Solver::WritePolygon(std::ofstream& ofs, float* pv) {

	ofs << "facet normal 0 0 0" << std::endl;
	ofs << "outer loop" << std::endl;
	ofs << "vertex " << pv[0] << " " << pv[1] << " " << pv[2] << std::endl;
	ofs << "vertex " << pv[3] << " " << pv[4] << " " << pv[5] << std::endl;
	ofs << "vertex " << pv[6] << " " << pv[7] << " " << pv[8] << std::endl;
	ofs << "endloop" << std::endl;
	ofs << "endfacet" << std::endl;

}

void Solver::PrintHole(int id) {
	mkdir("STL", 0755);

	std::ostringstream ossFileName;
	ossFileName << "./STL/";
	ossFileName << "data-hole";
	ossFileName << "-";
	ossFileName.width(5);
	ossFileName.setf(std::ios::fixed);
	ossFileName.fill('0');
	ossFileName << myrank;
	ossFileName << "-";
	ossFileName.width(1);
	ossFileName.setf(std::ios::fixed);
	ossFileName.fill('0');
	ossFileName << id;
	ossFileName << ".stl";

	std::ofstream ofs;
	ofs.open(ossFileName.str().c_str(), std::ios::out);
	ofs.close();

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
		int nc[3] = {size.x + 2*vc, size.y + 2*vc, size.z + 2*vc};
		real dx = cellSize.x;
	
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

		ofs.open(ossFileName.str().c_str(), std::ios::app);
		ofs << "solid" << std::endl;
		for(int k=vc; k<=size.z+vc-1; k++) {
			for(int j=vc; j<=size.y+vc-1; j++) {
				for(int i=vc; i<=size.x+vc-1; i++) {
					int mp = i + nc[0]*( j + nc[1]*k );
					int mw = i-1 + nc[0]*( j + nc[1]*k );
					int me = i+1 + nc[0]*( j + nc[1]*k );
					int ms = i + nc[0]*( j-1 + nc[1]*k );
					int mn = i + nc[0]*( j+1 + nc[1]*k );
					int mb = i + nc[0]*( j + nc[1]*(k-1) );
					int mt = i + nc[0]*( j + nc[1]*(k+1) );

					float x0 = origin[0] + (i - vc)*dx;
					float y0 = origin[1] + (j - vc)*dx;
					float z0 = origin[2] + (k - vc)*dx;
					float x1 = origin[0] + (i + 1 - vc)*dx;
					float y1 = origin[1] + (j + 1 - vc)*dx;
					float z1 = origin[2] + (k + 1 - vc)*dx;

					float v[9];
					float x[2] = {x0, x1};
					float y[2] = {y0, y1};
					float z[2] = {z0, z1};

					int cidp0 = pCutId0[mp];
					int cidp1 = pCutId1[mp];
					int cidp2 = pCutId2[mp];
					int cidp3 = pCutId3[mp];
					int cidp4 = pCutId4[mp];
					int cidp5 = pCutId5[mp];

					int cidw0 = pCutId0[mw];
					int cidw1 = pCutId1[mw];
					int cidw2 = pCutId2[mw];
					int cidw3 = pCutId3[mw];
					int cidw4 = pCutId4[mw];
					int cidw5 = pCutId5[mw];

					int cide0 = pCutId0[me];
					int cide1 = pCutId1[me];
					int cide2 = pCutId2[me];
					int cide3 = pCutId3[me];
					int cide4 = pCutId4[me];
					int cide5 = pCutId5[me];

					int cids0 = pCutId0[ms];
					int cids1 = pCutId1[ms];
					int cids2 = pCutId2[ms];
					int cids3 = pCutId3[ms];
					int cids4 = pCutId4[ms];
					int cids5 = pCutId5[ms];

					int cidn0 = pCutId0[mn];
					int cidn1 = pCutId1[mn];
					int cidn2 = pCutId2[mn];
					int cidn3 = pCutId3[mn];
					int cidn4 = pCutId4[mn];
					int cidn5 = pCutId5[mn];

					int cidb0 = pCutId0[mb];
					int cidb1 = pCutId1[mb];
					int cidb2 = pCutId2[mb];
					int cidb3 = pCutId3[mb];
					int cidb4 = pCutId4[mb];
					int cidb5 = pCutId5[mb];

					int cidt0 = pCutId0[mt];
					int cidt1 = pCutId1[mt];
					int cidt2 = pCutId2[mt];
					int cidt3 = pCutId3[mt];
					int cidt4 = pCutId4[mt];
					int cidt5 = pCutId5[mt];

					if( cidp1 == 0 ) {
						if( (cids1 != 0 || cidp2 !=0 || cide2 !=0) && 
								(cidn1 != 0 || cidp3 !=0 || cide3 !=0) && 
								(cidb1 != 0 || cidp4 !=0 || cide4 !=0) && 
								(cidt1 != 0 || cidp5 !=0 || cide5 !=0) ) {
							v[0] = x[1];
							v[1] = y[0];
							v[2] = z[0];
							v[3] = x[1];
							v[4] = y[0];
							v[5] = z[1];
							v[6] = x[1];
							v[7] = y[1];
							v[8] = z[0];
							WritePolygon(ofs, v);

							v[0] = x[1];
							v[1] = y[1];
							v[2] = z[1];
							v[6] = x[1];
							v[7] = y[0];
							v[8] = z[1];
							v[3] = x[1];
							v[4] = y[1];
							v[5] = z[0];
							WritePolygon(ofs, v);
						}
					}					

					if( cidp3 == 0 ) {
						if( (cidw3 != 0 || cidp0 !=0 || cidn0 !=0) &&
								(cide3 != 0 || cidp1 !=0 || cidn1 !=0) &&
								(cidb3 != 0 || cidp4 !=0 || cidn4 !=0) &&
								(cidt3 != 0 || cidp5 !=0 || cidn5 !=0) ) {
							v[0] = x[0];
							v[1] = y[1];
							v[2] = z[0];
							v[6] = x[0];
							v[7] = y[1];
							v[8] = z[1];
							v[3] = x[1];
							v[4] = y[1];
							v[5] = z[0];
							WritePolygon(ofs, v);

							v[0] = x[1];
							v[1] = y[1];
							v[2] = z[1];
							v[6] = x[1];
							v[7] = y[1];
							v[8] = z[0];
							v[3] = x[0];
							v[4] = y[1];
							v[5] = z[1];
							WritePolygon(ofs, v);
						}
					}

					if( cidp5 == 0 ) {
						if( (cidw5 != 0 || cidp0 !=0 || cidt0 !=0) &&
								(cide5 != 0 || cidp1 !=0 || cidt1 !=0) &&
								(cids5 != 0 || cidp2 !=0 || cidt2 !=0) &&
								(cidn5 != 0 || cidp3 !=0 || cidt3 !=0) ) {
							v[0] = x[0];
							v[1] = y[0];
							v[2] = z[1];
							v[6] = x[0];
							v[7] = y[1];
							v[8] = z[1];
							v[3] = x[1];
							v[4] = y[0];
							v[5] = z[1];
							WritePolygon(ofs, v);

							v[0] = x[1];
							v[1] = y[1];
							v[2] = z[1];
							v[3] = x[0];
							v[4] = y[1];
							v[5] = z[1];
							v[6] = x[1];
							v[7] = y[0];
							v[8] = z[1];
							WritePolygon(ofs, v);
						}
					}					
				}
			}
		}
		ofs << "endsolid" << std::endl;
		ofs.close();
	}
}

void Solver::PrintCut(int id) {
	mkdir("STL", 0755);

	std::ostringstream ossFileName;
	ossFileName << "./STL/";
	ossFileName << "data-cut";
	ossFileName << "-";
	ossFileName.width(5);
	ossFileName.setf(std::ios::fixed);
	ossFileName.fill('0');
	ossFileName << myrank;
	ossFileName << "-";
	ossFileName.width(1);
	ossFileName.setf(std::ios::fixed);
	ossFileName.fill('0');
	ossFileName << id;
	ossFileName << ".stl";

	std::ofstream ofs;
	ofs.open(ossFileName.str().c_str(), std::ios::out);
	ofs.close();

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
		int nc[3] = {size.x + 2*vc, size.y + 2*vc, size.z + 2*vc};
		real dx = cellSize.x;
	
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

		ofs.open(ossFileName.str().c_str(), std::ios::app);
		ofs << "solid" << std::endl;
		for(int k=vc; k<=size.z+vc-1; k++) {
			for(int j=vc; j<=size.y+vc-1; j++) {
				for(int i=vc; i<=size.x+vc-1; i++) {
					int mp = i + nc[0]*( j + nc[1]*k );

					float x0 = origin[0] + (i - vc)*dx;
					float y0 = origin[1] + (j - vc)*dx;
					float z0 = origin[2] + (k - vc)*dx;
					float x1 = origin[0] + (i + 1 - vc)*dx;
					float y1 = origin[1] + (j + 1 - vc)*dx;
					float z1 = origin[2] + (k + 1 - vc)*dx;

					float v[9];
					float x[2] = {x0, x1};
					float y[2] = {y0, y1};
					float z[2] = {z0, z1};

					int cidp0 = pCutId0[mp];
					int cidp1 = pCutId1[mp];
					int cidp2 = pCutId2[mp];
					int cidp3 = pCutId3[mp];
					int cidp4 = pCutId4[mp];
					int cidp5 = pCutId5[mp];

					if( cidp0 != 0 ) {
					}
					if( cidp1 != 0 ) {
						v[0] = x[1];
						v[1] = y[0];
						v[2] = z[0];
						v[3] = x[1];
						v[4] = y[0];
						v[5] = z[1];
						v[6] = x[1];
						v[7] = y[1];
						v[8] = z[0];
						WritePolygon(ofs, v);

						v[0] = x[1];
						v[1] = y[1];
						v[2] = z[1];
						v[6] = x[1];
						v[7] = y[0];
						v[8] = z[1];
						v[3] = x[1];
						v[4] = y[1];
						v[5] = z[0];
						WritePolygon(ofs, v);
					}
					if( cidp2 != 0 ) {
					}
					if( cidp3 != 0 ) {
						v[0] = x[0];
						v[1] = y[1];
						v[2] = z[0];
						v[6] = x[0];
						v[7] = y[1];
						v[8] = z[1];
						v[3] = x[1];
						v[4] = y[1];
						v[5] = z[0];
						WritePolygon(ofs, v);

						v[0] = x[1];
						v[1] = y[1];
						v[2] = z[1];
						v[6] = x[1];
						v[7] = y[1];
						v[8] = z[0];
						v[3] = x[0];
						v[4] = y[1];
						v[5] = z[1];
						WritePolygon(ofs, v);
					}
					if( cidp4 != 0 ) {
					}
					if( cidp5 != 0 ) {
						v[0] = x[0];
						v[1] = y[0];
						v[2] = z[1];
						v[6] = x[0];
						v[7] = y[1];
						v[8] = z[1];
						v[3] = x[1];
						v[4] = y[0];
						v[5] = z[1];
						WritePolygon(ofs, v);

						v[0] = x[1];
						v[1] = y[1];
						v[2] = z[1];
						v[3] = x[0];
						v[4] = y[1];
						v[5] = z[1];
						v[6] = x[1];
						v[7] = y[0];
						v[8] = z[1];
						WritePolygon(ofs, v);
					}
				
				}
			}
		}
		ofs << "endsolid" << std::endl;
		ofs.close();
	}
}

void Solver::Dump(const int step) {
	plsUX0->Dump2(blockManager, step, "ux");
	plsUY0->Dump2(blockManager, step, "uy");
	plsUZ0->Dump2(blockManager, step, "uz");

	plsP0->Dump2(blockManager, step, "p");

	plsVw->Dump2(blockManager, step, "vw");
	plsVe->Dump2(blockManager, step, "ve");
	plsVs->Dump2(blockManager, step, "vs");
	plsVn->Dump2(blockManager, step, "vn");
	plsVb->Dump2(blockManager, step, "vb");
	plsVt->Dump2(blockManager, step, "vt");

	plsT0->Dump2(blockManager, step, "t");

	plsUXCP->Dump2(blockManager, step, "uxcp");
	plsUYCP->Dump2(blockManager, step, "uycp");
	plsUZCP->Dump2(blockManager, step, "uzcp");
	plsTCP->Dump2(blockManager, step, "tcp");

	plsUXCP->Dump2(blockManager, step, "uxcp");
	plsUYCP->Dump2(blockManager, step, "uycp");
	plsUZCP->Dump2(blockManager, step, "uzcp");
	plsTCP->Dump2(blockManager, step, "tcp");

	if( !strcasecmp(g_pFFVConfig->TimeIntegrationMethodForFlow.c_str(), "explicit") ) {
		plsUXDP->Dump2(blockManager, step, "uxdp");
		plsUYDP->Dump2(blockManager, step, "uydp");
		plsUZDP->Dump2(blockManager, step, "uzdp");
		plsTDP->Dump2(blockManager, step, "tdp");
	}
}

void Solver::Load(const int step) {
	plsUX0->Load2(blockManager, step, "ux");
	plsUY0->Load2(blockManager, step, "uy");
	plsUZ0->Load2(blockManager, step, "uz");

	plsP0->Load2(blockManager, step, "p");

	plsVw->Load2(blockManager, step, "vw");
	plsVe->Load2(blockManager, step, "ve");
	plsVs->Load2(blockManager, step, "vs");
	plsVn->Load2(blockManager, step, "vn");
	plsVb->Load2(blockManager, step, "vb");
	plsVt->Load2(blockManager, step, "vt");

	plsT0->Load2(blockManager, step, "t");

	plsUXCP->Load2(blockManager, step, "uxcp");
	plsUYCP->Load2(blockManager, step, "uycp");
	plsUZCP->Load2(blockManager, step, "uzcp");
	plsTCP->Load2(blockManager, step, "tcp");

	if( !strcasecmp(g_pFFVConfig->TimeIntegrationMethodForFlow.c_str(), "explicit") ) {
		plsUXDP->Load2(blockManager, step, "uxdp");
		plsUYDP->Load2(blockManager, step, "uydp");
		plsUZDP->Load2(blockManager, step, "uzdp");
		plsTDP->Load2(blockManager, step, "tdp");
	}
}

