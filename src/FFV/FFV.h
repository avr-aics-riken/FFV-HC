#ifndef FFV_H
#define FFV_H

#include "FFVConfig.h"
#include "FFVPerfMonitor.h"
#include "FFVGrid.h"
#include "FFVSolver.h"

#include <string>

class FFV {
public:
	FFV(const char* filenameConfig);
	~FFV();

private:
	std::string filenameConfig;

private:
	FFVConfig*      pFFVConfig;
	FFVPerfMonitor* pFFVPerfMonitor;
	FFVGrid*        pFFVGrid;
	FFVSolver*      pFFVSolver;

public:
	int Init();
	int Loop();
	int Post();

private:
	void Init_ParallelEnv();
	void Init_Config();
	void Init_PerfMonitor();
	void Init_Grid();
	void Init_Solver();

	void Post_ParallelEnv();
	void Post_Config();
	void Post_PerfMonitor();
	void Post_Grid();
	void Post_Solver();
};


#endif

