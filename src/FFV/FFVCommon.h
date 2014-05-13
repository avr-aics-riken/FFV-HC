#ifndef FFVCOMMON_H
#define FFVCOMMON_H

void PrintLog(int rank, int level, const char* format, ...);
void PrintLog(int level, const char* format, ...);
double GetTime();
double RandomUniform();
double RandomNormal(double mu, double sigma);

int GetNumThreads();
int GetNumProcesses();
int GetMPIRank();
int GetMPISize();

#endif

