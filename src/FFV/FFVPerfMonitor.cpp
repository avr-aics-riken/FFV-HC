#include "FFVPerfMonitor.h"
#include "FFVCommon.h"
#include "FFVGlobalVars.h"

#include <PerfMonitor.h>

#include <mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef __FUJITSU
#include <fj_tool/fapp.h>
#endif

char* FFV_label[] = {
	"Init_LoadPolygon",
	"Init_DivideDomain",
	"Init_OrderBlock",
	"Init_CommOctree",
	"Init_PartitionBlock",
	"Init_RegisterBlock",
	"Init_DistributePolygon",
	"Init_RepairPolygon",
	"Init_CalcCutInfo",
	"Init_CalcCutInfo01",
	"Init_CalcCutInfo02",
	"Init_CalcCutInfo03",
	"Init_CalcCutInfo04",
	"Init_CalcCutInfo05",
	"Init_CalcCutInfo06",
	"Init_Filling",
	"Init_GeometricalProperties",
	"Init_InitVars",
	"Update",
	"UpdateT",
	"UpdateT01",
	"UpdateT02",
	"UpdateT03",
	"UpdateT04",
	"UpdateT05",
	"UpdateUX",
	"UpdateUX01",
	"UpdateUX02",
	"UpdateUX03",
	"UpdateUX04",
	"UpdateUX05",
	"UpdateUY",
	"UpdateUY01",
	"UpdateUY02",
	"UpdateUY03",
	"UpdateUY04",
	"UpdateUY05",
	"UpdateUZ",
	"UpdateUZ01",
	"UpdateUZ02",
	"UpdateUZ03",
	"UpdateUZ04",
	"UpdateUZ05",
	"UpdateP",
	"UpdateP01",
	"UpdateP02",
	"UpdateP03",
	"UpdateP04",
	"UpdateP05",
	"UpdateP06",
	"UpdateP07",
	"UpdateU",
	"UpdateU01",
	"UpdateU02",
	"UpdateU03",
	"Print",
	"PrintTime",
	"PrintILS",
	"PrintStats",
	"PrintForce",
	"PrintHeatFlux",
	"PrintData",
};


FFVPerfMonitor::FFVPerfMonitor() {
}

FFVPerfMonitor::~FFVPerfMonitor() {
}

void FFVPerfMonitor::Init() {
	int rank        = GetMPIRank();
	int nThreads    = GetNumThreads();
	int nProcesses  = GetNumProcesses();

	pPerfMonitor = new pm_lib::PerfMonitor;
	pPerfMonitor->initialize(ffv_tm_END);
	pPerfMonitor->setRankInfo(rank);
	pPerfMonitor->setProperties(ffv_tm_Init_LoadPolygon,           "LoadPolygon",           pm_lib::PerfMonitor::COMM, true);
	pPerfMonitor->setProperties(ffv_tm_Init_DivideDomain,          "DivideDomain",          pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_Init_OrderBlock,            "OrderBlock",            pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_Init_CommOctree,            "CommOctree",            pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_Init_PartitionBlock,        "PartitionBlock",        pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_Init_RegisterBlock,         "RegisterBlock",         pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_Init_DistributePolygon,     "DistributePolygon",     pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_Init_RepairPolygon,         "RepairPolygon",         pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_Init_CalcCutInfo,           "CalcCutInfo",           pm_lib::PerfMonitor::CALC, false);
	pPerfMonitor->setProperties(ffv_tm_Init_CalcCutInfo01,         "CalcCutInfo01",         pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_Init_CalcCutInfo02,         "CalcCutInfo02",         pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_Init_CalcCutInfo03,         "CalcCutInfo03",         pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_Init_CalcCutInfo04,         "CalcCutInfo04",         pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_Init_CalcCutInfo05,         "CalcCutInfo05",         pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_Init_CalcCutInfo06,         "CalcCutInfo06",         pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_Init_Filling,               "Filling",               pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_Init_GeometricalProperties, "GeometricalProperties", pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_Init_InitVars,              "InitVars",              pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_Update,     "Update",      pm_lib::PerfMonitor::CALC, false);
	pPerfMonitor->setProperties(ffv_tm_UpdateT,    "UpdateT",     pm_lib::PerfMonitor::CALC, false);
	pPerfMonitor->setProperties(ffv_tm_UpdateUX,   "UpdateUX",    pm_lib::PerfMonitor::CALC, false);
	pPerfMonitor->setProperties(ffv_tm_UpdateUY,   "UpdateUY",    pm_lib::PerfMonitor::CALC, false);
	pPerfMonitor->setProperties(ffv_tm_UpdateUZ,   "UpdateUZ",    pm_lib::PerfMonitor::CALC, false);
	pPerfMonitor->setProperties(ffv_tm_UpdateP,    "UpdateP",     pm_lib::PerfMonitor::CALC, false);
	pPerfMonitor->setProperties(ffv_tm_UpdateU,    "UpdateU",     pm_lib::PerfMonitor::CALC, false);
	pPerfMonitor->setProperties(ffv_tm_UpdateT01,  "UpdateT01",   pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_UpdateT02,  "UpdateT02",   pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_UpdateT03,  "UpdateT03",   pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_UpdateT04,  "UpdateT04",   pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_UpdateT05,  "UpdateT05",   pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_UpdateUX01, "UpdateUX01",  pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_UpdateUX02, "UpdateUX02",  pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_UpdateUX03, "UpdateUX03",  pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_UpdateUX04, "UpdateUX04",  pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_UpdateUX05, "UpdateUX05",  pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_UpdateUY01, "UpdateUY01",  pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_UpdateUY02, "UpdateUY02",  pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_UpdateUY03, "UpdateUY03",  pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_UpdateUY04, "UpdateUY04",  pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_UpdateUY05, "UpdateUY05",  pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_UpdateUZ01, "UpdateUZ01",  pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_UpdateUZ02, "UpdateUZ02",  pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_UpdateUZ03, "UpdateUZ03",  pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_UpdateUZ04, "UpdateUZ04",  pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_UpdateUZ05, "UpdateUZ05",  pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_UpdateP01,  "UpdateP01",   pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_UpdateP02,  "UpdateP02",   pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_UpdateP03,  "UpdateP03",   pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_UpdateP04,  "UpdateP04",   pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_UpdateP05,  "UpdateP05",   pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_UpdateP06,  "UpdateP06",   pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_UpdateP07,  "UpdateP07",   pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_UpdateU01,  "UpdateU01",   pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_UpdateU02,  "UpdateU02",   pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_UpdateU03,  "UpdateU03",   pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_Print,         "Print",         pm_lib::PerfMonitor::CALC, false);
	pPerfMonitor->setProperties(ffv_tm_PrintTime,     "PrintTime",     pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_PrintILS,      "PrintILS",      pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_PrintStats,    "PrintStats",    pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_PrintForce,    "PrintForce",    pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_PrintHeatFlux, "PrintHeatFlux", pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setProperties(ffv_tm_PrintData,     "PrintData",     pm_lib::PerfMonitor::CALC, true);
	pPerfMonitor->setParallelMode("Hybrid", nThreads, nProcesses);
}

void FFVPerfMonitor::Fin() {
	pPerfMonitor->gather();

	int rank = GetMPIRank();

	if( rank == 0 ) {
		char hostname[MPI_MAX_PROCESSOR_NAME] = {0};
		int nameLen;
		MPI_Get_processor_name(hostname, &nameLen);
		nameLen++;

		FILE* fp = fopen(g_pFFVConfig->OutputLogFilenameProfiling.c_str(), "w");
		pPerfMonitor->print(fp, hostname, g_pFFVConfig->OperatorName.c_str());
		pPerfMonitor->printDetail(fp);
		fclose(fp);
	}
}

void FFVPerfMonitor::Start(unsigned key, int level, int id, bool bBarrier) {
	if( bBarrier == true && g_pFFVConfig->OperationMode == "Benchmark" ) {
		MPI_Barrier(MPI_COMM_WORLD);
	}
	pPerfMonitor->start(key);

#ifdef __FAPP
#ifdef __FUJITSU
	fapp_start(FFV_label[key], 0, 0);
#endif
#endif
}

void FFVPerfMonitor::Stop(unsigned key, int level, int id, double flopPerTask, unsigned iterationCount) {
#ifdef __FAPP
#ifdef __FUJITSU
	fapp_stop(FFV_label[key], 0, 0);
#endif
#endif

	pPerfMonitor->stop(key, flopPerTask, iterationCount);
}

