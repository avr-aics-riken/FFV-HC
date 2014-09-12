#############################################################################
#
# FFV-HC
# FrontFlow/violet for hierarchical Cartesian grids
#
# Copyright (c) 2013,2014 Institute of Industrial Science, The University of Tokyo. 
# All rights reserved.
#
# Copyright (c) 2013,2014 Advanced Institute for Computational Science, RIKEN. 
# All rights reserved.
#
#############################################################################

.PHONY: depend clean

all:
	(cd CGMfuncs; make)
	(cd FFV; make)

clean:
	(cd CGMfuncs; make clean)
	(cd FFV; make clean)

depend:
	(cd CGMfuncs; make depend)
	(cd FFV; make depend)

