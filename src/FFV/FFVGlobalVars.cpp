#include "FFVGlobalVars.h"

void PrintLog(int myrank, int level, const char* format, ...) {
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

