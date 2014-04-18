#ifndef PM_H
#define PM_H

enum timing_key {
	tm_Init_LoadSTL,
	tm_Init_DivideDomain,
	tm_Init_CreateTree,
	tm_Init_DistributeSTL,
	tm_Init_RepairPolygonData,
	tm_Init_CalcCutInfo,
	tm_Init_CalcCutInfo01,
	tm_Init_CalcCutInfo02,
	tm_Init_CalcCutInfo03,
	tm_Init_CalcCutInfo04,
	tm_Init_CalcCutInfo05,
	tm_Init_CalcCutInfo06,
	tm_Init_Filling,
	tm_Init_GeometricalProperties,
	tm_Init_InitVars,
	tm_Update,
	tm_UpdateT,
	tm_UpdateT01,
	tm_UpdateT02,
	tm_UpdateT03,
	tm_UpdateT04,
	tm_UpdateT05,
	tm_UpdateUX,
	tm_UpdateUX01,
	tm_UpdateUX02,
	tm_UpdateUX03,
	tm_UpdateUX04,
	tm_UpdateUX05,
	tm_UpdateUY,
	tm_UpdateUY01,
	tm_UpdateUY02,
	tm_UpdateUY03,
	tm_UpdateUY04,
	tm_UpdateUY05,
	tm_UpdateUZ,
	tm_UpdateUZ01,
	tm_UpdateUZ02,
	tm_UpdateUZ03,
	tm_UpdateUZ04,
	tm_UpdateUZ05,
	tm_UpdateP,
	tm_UpdateP01,
	tm_UpdateP02,
	tm_UpdateP03,
	tm_UpdateP04,
	tm_UpdateP05,
	tm_UpdateP06,
	tm_UpdateP07,
	tm_UpdateU,
	tm_UpdateU01,
	tm_UpdateU02,
	tm_UpdateU03,
	tm_Print,
	tm_PrintTime,
	tm_PrintILS,
	tm_PrintStats,
	tm_PrintForce,
	tm_PrintHeatFlux,
	tm_PrintData,
	tm_END,
};

void PM_Start(unsigned key, int level=0, int id=0, bool bBarrier=false);
void PM_Stop(unsigned key, int level=0, int id=0, double flopPerTask=0.0, unsigned iterationCount=1);

#endif

