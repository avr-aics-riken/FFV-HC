#include "FFVConfig.h"
#include <PerfMonitor.h>

#ifdef GLOBAL_VALUE_DEFINE
#define GLOBAL_VARIABLE
#define GLOBAL_VALUE(v) = (v)
#else
#define GLOBAL_VARIABLE extern
#define GLOBAL_VALUE(v)
#endif

GLOBAL_VARIABLE FFVConfig* g_pFFVConfig;
GLOBAL_VARIABLE pm_lib::PerfMonitor* g_pPM;

#undef GLOBAL_VARIABLE
#undef GLOBAL_VALUE

void PrintLog(int myrank, int level, const char* format, ...);
double RandomUniform();
double RandomNormal(double mu, double sigma);


