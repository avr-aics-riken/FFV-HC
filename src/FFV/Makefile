#############################################################################
#
# FFV-HC
# FrontFlow/violet based on Building-cube Method
#
# Copyright (c) 2013,2014 Institute of Industrial Science, The University of Tokyo. 
# All rights reserved.
#
# Copyright (c) 2013,2014 Advanced Institute for Computational Science, RIKEN. 
# All rights reserved.
#
#############################################################################

-include ../make_setting

FFV_INC_PATH    := -I$(PMLIB)/include
FFV_INC_PATH    += -I$(TPLIB)/include
FFV_INC_PATH    += -I$(POLYLIB)/include
FFV_INC_PATH    += -I$(CUTLIB)/include -DCUTINFO_NORMAL_DOUBLE
FFV_INC_PATH    += -I$(BCM_DIR)/include
FFV_INC_PATH    += -I$(HDMLIB)/include
FFV_INC_PATH    += -I$(MPI_DIR)/include
FFV_INC_PATH    += -I$(CGM_FUNCS)

FFV_LIB_PATH    := -L$(PMLIB)/lib
FFV_LIB_PATH    += -L$(TPLIB)/lib
FFV_LIB_PATH    += -L$(POLYLIB)/lib
FFV_LIB_PATH    += -L$(CUTLIB)/lib
FFV_LIB_PATH    += -L$(BCM_DIR)/lib
FFV_LIB_PATH    += -L$(HDMLIB)/lib
FFV_LIB_PATH    += -L$(MPI_DIR)/lib
FFV_LIB_PATH    += -L$(CGM_FUNCS)

FFV_LIBS        := -lHDM -lbcm -lCUT -lMPIPOLY -lPM -lMPITP -lcgmfuncs -lmpi -lmpi_cxx -lmpi_f90 -lmpi_f77
FFV_LIBS        := -lHDM -lBCM -lCUT -lPOLYmpi -lPM -lTPmpi -lcgmfuncs -lmpi -lmpi_cxx -lmpi_f90 -lmpi_f77

PROGRAM = ffv

FFV_SOLVERNAME	:= "FrontFlow/violet-HC"
FFV_VERSION			:= "1.3.0"
CXXFLAGS				+= -DFFV_SOLVERNAME=\"$(FFV_SOLVERNAME)\"
CFLAGS					+= -DFFV_SOLVERNAME=\"$(FFV_SOLVERNAME)\"
CXXFLAGS				+= -DFFV_VERSION=\"$(FFV_VERSION)\"
CFLAGS					+= -DFFV_VERSION=\"$(FFV_VERSION)\"
#FFV_REVISION := $(shell sh -c 'git describe --always --dirty')
#FFV_REVISION := $(shell svnversion -n .)
#CXXFLAGS += -DFFV_REVISION=\"$(FFV_REVISION)\"
#CFLAGS   += -DFFV_REVISION=\"$(FFV_REVISION)\"

BUILD_DATE   := $(shell sh -c 'date +%s')
CXXFLAGS += -DBUILD_DATE=\"$(BUILD_DATE)\"
CFLAGS   += -DBUILD_DATE=\"$(BUILD_DATE)\"

#vpath %.cpp $(BCM_FILEIO)/src
#vpath %.cpp $(BCM_UTILS)/src

CSRC            := $(wildcard *.cpp)
CHDR            := $(wildcard *.h  )
FSRC            := $(wildcard *.f90)
OBJS            := $(CSRC:.cpp=.o) $(FSRC:.f90=.o)
#BCM_UTILS_OBJS  := BCMOctree.o BCMPolylib.o
#BCM_FILEIO_OBJS := BitVoxel.o ErrorUtil.o BCMFileSaver.o LeafBlockSaver.o FileSystemUtil.o

PROGRAM: $(OBJS) $(BCM_UTILS_OBJS) $(BCM_FILEIO_OBJS)
	$(CXX) \
	$(LDFLAGS) \
	-o $(PROGRAM) $(OBJS) $(BCM_UTILS_OBJS) $(BCM_FILEIO_OBJS) \
	$(FFV_LIB_PATH) $(FFV_LIBS) $(LIBS) 
	mv $(PROGRAM) ../../bin/

.cpp.o:
	$(CXX) $(CXXFLAGS) $(FFV_INC_PATH) -o $@ -c $<

.SUFFIXES: .f90
.f90.o:
	$(F90) $(F90FLAGS) $(FFV_INC_PATH) -o $@ -c $<

clean:
	rm -rf *.o ../../bin/$(PROGRAM) depend.inc

cleandata:
	rm -rf BIN STL VTK PLOT3D BCM_OUT CS data-*

depend:
	$(CXX) -MM $(CSRC) $(CXXFLAGS) $(FFV_INC_PATH) > depend.inc

-include depend.inc


