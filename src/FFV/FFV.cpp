#include "FFV.h"
#include "FFVCommon.h"

#define GLOBAL_VALUE_DEFINE
#include "FFVGlobalVars.h"
#undef GLOBAL_VALUE_DEFINE

#include <BCMTools.h>

#include <string>
#include <mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif

FFV::FFV(const char* filenameConfig) {
	this->filenameConfig = std::string(filenameConfig);
}

FFV::~FFV() {
}

int FFV::Init() {
	double t0 = GetTime();
	Init_ParallelEnv();
	double t1 = GetTime();

	double t2 = GetTime();
	Init_Config();
	double t3 = GetTime();

	double t4 = GetTime();
	Init_PerfMonitor();
	double t5 = GetTime();

/*
	time_t now = time(NULL);
	const char* exec_date = ctime(&now);
*/

	time_t now = atol(BUILD_DATE);
	char* build_date = ctime(&now);
	build_date[strlen(build_date)-1] = ' ';

	PrintLog(1,                                     FFV_SOLVERNAME);
	PrintLog(2, "%-20s : %s", "Version"           , FFV_VERSION);
	PrintLog(2, "%-20s : %s", "Build date"        , build_date);
//	PrintLog(2, "%-20s : %s", "Execution date"    , exec_date);
	PrintLog(2, "%-20s : %s", "Configuration file", this->filenameConfig.c_str());
	PrintLog(2, "%-20s : %d", "MPI processes"     , g_nProcesses);
	PrintLog(2, "%-20s : %d", "OpenMP threads"    , g_nThreads);
	PrintLog(2, "%-20s : %d", "Real in Fortran"		, sizeof_real_());
	PrintLog(2, "%-20s : %d", "Real in C/C++  "		, sizeof(real));

	double t6 = GetTime();
	Init_Grid();
	double t7 = GetTime();

	double t8 = GetTime();
	Init_Solver();
	double t9 = GetTime();

	return EX_SUCCESS;
}

int FFV::Loop() {
	int StepStart = pFFVConfig->TimeControlSessionStartD;
	int StepEnd   = pFFVConfig->TimeControlSessionEndD;

	if(StepStart == 0) {
		pFFVSolver->Print(0);
	} else {
		pFFVSolver->Load(StepStart);
	}

	for(int step=StepStart+1; step<=StepEnd; step++) {
		pFFVSolver->Update(step);
		pFFVSolver->Print(step);
		pFFVSolver->Dump(step);
	}

	return EX_SUCCESS;
}

int FFV::Post() {
	Post_ParallelEnv();
	Post_Config();
	Post_PerfMonitor();
	Post_Grid();
	Post_Solver();
	return EX_SUCCESS;
}

void FFV::Init_ParallelEnv() {
	g_rank        = GetMPIRank();
	g_nProcesses  = GetNumProcesses();
	g_nThreads    = GetNumThreads();
}

void FFV::Init_Config() {
	pFFVConfig = new FFVConfig();
	pFFVConfig->Load(this->filenameConfig);
	pFFVConfig->Check();

	g_pFFVConfig = pFFVConfig;
}

void FFV::Init_PerfMonitor() {
	pFFVPerfMonitor = new FFVPerfMonitor();
	pFFVPerfMonitor->Init();

	g_pFFVPerfMonitor = pFFVPerfMonitor;
}

void FFV::Init_Grid() {
	pFFVGrid = new FFVGrid();
	pFFVGrid->Init();
}

void FFV::Init_Solver() {
	pFFVSolver = new FFVSolver();
	pFFVSolver->Init();
}

void FFV::Post_ParallelEnv() {
}

void FFV::Post_Config() {
}

void FFV::Post_PerfMonitor() {
	g_pFFVPerfMonitor->Fin();
}

void FFV::Post_Grid() {
}

void FFV::Post_Solver() {
}

