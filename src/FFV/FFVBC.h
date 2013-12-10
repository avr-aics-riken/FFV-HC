#ifndef BC_H
#define BC_H

#include <string>

enum BoundaryType {
	INNER					= -1,  // inner boundary
	DIRICHLET			= 0,   // dirichlet boundary condition
	NEUMANN				= 1,   // neumann boundary condition
	PERIODIC			= 9,   // periodic boundary
	POISEUILLE_U	= 100, // special-purpose boundary condition
	POISEUILLE_P	= 101, // special-purpose boundary condition
	NONE,
};

typedef struct _OBC {
	int type;
	double value;
	double xc;
	double yc;
	double zc;
	double rc;
	double lc;
	double hc;
	double wc;
}OBC;

BoundaryType GetBoundaryType(std::string bctype);

#endif

