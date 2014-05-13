#ifndef FFVPERFMONITOR_H
#define FFVPERFMONITOR_H

#include <PerfMonitor.h>

class FFVPerfMonitor {
public:
	FFVPerfMonitor();
	~FFVPerfMonitor();

private:
	pm_lib::PerfMonitor* pPerfMonitor;

public:
	void Init();
	void Fin();
	void Start(unsigned key, int level=0, int id=0, bool bBarrier=false);
	void Stop(unsigned key, int level=0, int id=0, double flopPerTask=0.0, unsigned iterationCount=1);
};

enum FFV_timing_key {
	ffv_tm_Init_LoadPolygon,
	ffv_tm_Init_DivideDomain,
	ffv_tm_Init_OrderBlock,
	ffv_tm_Init_CommOctree,
	ffv_tm_Init_PartitionBlock,
	ffv_tm_Init_RegisterBlock,
	ffv_tm_Init_DistributePolygon,
	ffv_tm_Init_RepairPolygon,
	ffv_tm_Init_CalcCutInfo,
	ffv_tm_Init_CalcCutInfo01,
	ffv_tm_Init_CalcCutInfo02,
	ffv_tm_Init_CalcCutInfo03,
	ffv_tm_Init_CalcCutInfo04,
	ffv_tm_Init_CalcCutInfo05,
	ffv_tm_Init_CalcCutInfo06,
	ffv_tm_Init_Filling,
	ffv_tm_Init_GeometricalProperties,
	ffv_tm_Init_InitVars,
	ffv_tm_Update,
	ffv_tm_UpdateT,
	ffv_tm_UpdateT01,
	ffv_tm_UpdateT02,
	ffv_tm_UpdateT03,
	ffv_tm_UpdateT04,
	ffv_tm_UpdateT05,
	ffv_tm_UpdateUX,
	ffv_tm_UpdateUX01,
	ffv_tm_UpdateUX02,
	ffv_tm_UpdateUX03,
	ffv_tm_UpdateUX04,
	ffv_tm_UpdateUX05,
	ffv_tm_UpdateUY,
	ffv_tm_UpdateUY01,
	ffv_tm_UpdateUY02,
	ffv_tm_UpdateUY03,
	ffv_tm_UpdateUY04,
	ffv_tm_UpdateUY05,
	ffv_tm_UpdateUZ,
	ffv_tm_UpdateUZ01,
	ffv_tm_UpdateUZ02,
	ffv_tm_UpdateUZ03,
	ffv_tm_UpdateUZ04,
	ffv_tm_UpdateUZ05,
	ffv_tm_UpdateP,
	ffv_tm_UpdateP01,
	ffv_tm_UpdateP02,
	ffv_tm_UpdateP03,
	ffv_tm_UpdateP04,
	ffv_tm_UpdateP05,
	ffv_tm_UpdateP06,
	ffv_tm_UpdateP07,
	ffv_tm_UpdateU,
	ffv_tm_UpdateU01,
	ffv_tm_UpdateU02,
	ffv_tm_UpdateU03,
	ffv_tm_Print,
	ffv_tm_PrintTime,
	ffv_tm_PrintILS,
	ffv_tm_PrintStats,
	ffv_tm_PrintForce,
	ffv_tm_PrintHeatFlux,
	ffv_tm_PrintData,
	ffv_tm_END,
};

#endif

