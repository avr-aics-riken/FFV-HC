#ifndef FFVSOLVER_H
#define FFVSOLVER_H

#include <BlockManager.h>

class FFVSolver {
public:
	FFVSolver();
	~FFVSolver();

public:
	void Init();
	void Update(int step);
	void Load(int step);
	void Dump(int step);
	void Print(int step);

private:
	BlockManager& blockManager;
	int rank;

private:

private:

};

#endif

