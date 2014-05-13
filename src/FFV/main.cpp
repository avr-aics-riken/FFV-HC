#include "FFV.h"

#include <BCMTools.h>

#include <iostream>
#include <mpi.h>

int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);

	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if( argc != 2 ) {
		if( rank == 0 ) {
			std::cout << "usage: " << argv[0] << " configfile" << std::endl;
		}
		MPI_Abort(MPI_COMM_WORLD, EX_USAGE);
	}

	FFV *pFFV = new FFV(argv[1]);

	int nResultInit = pFFV->Init();
	switch( nResultInit ) {
		case EX_SUCCESS : {
			break;
		}
		case EX_FAILURE : {
			break;
		}
		default : {
			break;
		}
	}

	int nResultLoop = pFFV->Loop();
	switch( nResultLoop ) {
		case EX_SUCCESS : {
			break;
		}
		case EX_FAILURE : {
			break;
		}
		default : {
			break;
		}
	}

	int nResultPost = pFFV->Post();
	switch( nResultPost ) {
		case EX_SUCCESS : {
			break;
		}
		case EX_FAILURE : {
			break;
		}
		default : {
			break;
		}
	}

	delete pFFV;

	MPI_Finalize();

  return EX_SUCCESS;
}

