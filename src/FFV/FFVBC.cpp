#include "FFVBC.h"

#include <stdlib.h>
#include <string.h>

BoundaryType GetBoundaryType(std::string bctype) {
	if( !strcasecmp(bctype.c_str(), "Neumann") ) {
		return NEUMANN;
	} else if( !strcasecmp(bctype.c_str(), "Dirichlet") ) {
		return DIRICHLET;
	} else if( !strcasecmp(bctype.c_str(), "Periodic") ) {
		return PERIODIC;
	} else if( !strcasecmp(bctype.c_str(), "PoiseuilleU") ) {
		return POISEUILLE_U;
	} else if( !strcasecmp(bctype.c_str(), "PoiseuilleP") ) {
		return POISEUILLE_P;
	}
	return NONE;
}

