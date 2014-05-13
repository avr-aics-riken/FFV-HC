#include "FFVGlobalVars.h"

#include <sstream>
#include <fstream>
#include <sys/time.h>
#include <mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif

//void PrintLog(int rank, int level, const char* format, ...) {
void PrintLog(int level, const char* format, ...) {
	if( g_rank != 0 ) {
		return;
	}

	static int first = 1;

	va_list ap;
	va_start(ap, format);

	char* buffer;
	int size = vasprintf(&buffer, format, ap);
	va_end(ap);

	std::ostringstream ossMessage;
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
	ossMessage << std::string(buffer);

	free(buffer);

	std::cout << ossMessage.str() << std::endl;

	std::string filename = g_pFFVConfig->OutputLogFilenameBase;
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

double GetTime() {
	struct timeval tp;
	int i = gettimeofday(&tp, 0);
	return ((double)(tp.tv_sec) + (double)(tp.tv_usec)*1.0e-6);
}

double RandomUniform() {
	return rand()/(1.0 + RAND_MAX);
}

double RandomNormal(double mu, double sigma) {
	static int    flag = 0;
	static double save = 0.0;

	if( flag == 0 ) {
		double u_0 = RandomUniform();
		double u_1 = RandomUniform();
		double v_0 = mu + sigma*sqrt(-2.0*log(1.0 - u_0))*cos(2.0*M_PI*u_1);
		double v_1 = mu + sigma*sqrt(-2.0*log(1.0 - u_0))*sin(2.0*M_PI*u_1);
		save = v_1;
		flag = 1;
		return v_0;
	} else {
		flag = 0;
		return save;
	}

	return 0.0;
}

void PrintInt32t(int32_t i) {
	int m[32];
	for(int n=0; n<32; n++) {
		m[n] = (i >> n)%2;
	}
	for(int n=0; n<32; n++) {
		std::cout << m[31-n];
	}
}

int GetMPIRank() {
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	return rank;
}

int GetMPISize() {
	int size = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	return size;
}

int GetNumThreads() {
	int nThreads = 1;
#ifdef _OPENMP
	nThreads = omp_get_max_threads();
#endif
	return nThreads;
}

int GetNumProcesses() {
	return GetMPISize();
}

