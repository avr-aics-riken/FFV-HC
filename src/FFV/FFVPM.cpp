#include "FFVPM.h"
#include <PerfMonitor.h>
#include "FFVGlobalVars.h"

#ifdef __FUJITSU
#include <fj_tool/fapp.h>
#endif

char* label[] = {
	"Init_LoadSTL",
	"Init_DivideDomain",
	"Init_CreateTree",
	"Init_DistributeSTL",
	"Init_RepairPolygonData",
	"Init_CalcCutInfo",
	"Init_CalcCutInfo01",
	"Init_CalcCutInfo02",
	"Init_CalcCutInfo03",
	"Init_CalcCutInfo04",
	"Init_CalcCutInfo05",
	"Init_CalcCutInfo06",
	"Init_Filling",
	"Init_PartitioningRegions",
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


void PM_Start(unsigned key, int level, int id, bool bBarrier) {
	if( bBarrier == true && !strcasecmp(g_pFFVConfig->OperationMode.c_str(), "benchmark") ) {
		MPI_Barrier(MPI_COMM_WORLD);
	}
	g_pPM->start(key);

#ifdef __FAPP
#ifdef __FUJITSU
	fapp_start(label[key], 0, 0);
#endif
#endif
}

void PM_Stop(unsigned key, int level, int id, double flopPerTask, unsigned iterationCount) {
#ifdef __FAPP
#ifdef __FUJITSU
	fapp_stop(label[key], 0, 0);
#endif
#endif

	g_pPM->stop(key, flopPerTask, iterationCount);
}

