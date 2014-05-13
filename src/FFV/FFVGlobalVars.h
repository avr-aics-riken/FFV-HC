#include "FFVConfig.h"
#include "FFVPerfMonitor.h"

#include <PerfMonitor.h>

#ifdef GLOBAL_VALUE_DEFINE
#define GLOBAL_VARIABLE
#define GLOBAL_VALUE(v) = (v)
#else
#define GLOBAL_VARIABLE extern
#define GLOBAL_VALUE(v)
#endif

GLOBAL_VARIABLE FFVConfig*       g_pFFVConfig;
GLOBAL_VARIABLE FFVPerfMonitor*  g_pFFVPerfMonitor;

GLOBAL_VARIABLE int g_rank;
GLOBAL_VARIABLE int g_nThreads;
GLOBAL_VARIABLE int g_nProcesses;

GLOBAL_VARIABLE pm_lib::PerfMonitor* g_pPM;

#undef GLOBAL_VARIABLE
#undef GLOBAL_VALUE

