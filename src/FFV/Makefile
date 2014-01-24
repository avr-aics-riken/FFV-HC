#############################################################################
#
# FFV-BCM
# FrontFlow/violet based on Building-cube Method
#
# Copyright (c) 2013 Institute of Industrial Science, The University of Tokyo. 
# All rights reserved.
#
# Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN. 
# All rights reserved.
#
#############################################################################

# Default environment
FFVC=gnu

### compiler detection

ifeq ($(FFVC),gnu)
include ../make_setting.gnu
endif

ifeq ($(FFVC),fx)
include ../make_setting.fx
endif

ifeq ($(FFVC),ibm)
include ../make_setting.ibm
endif

ifeq ($(FFVC),intel)
include ../make_setting.intel
endif

#####

#UTIL    = ../util

UDEF_INC_PATH=-I$(POLYLIB_DIR)/include \
              -I$(CUTLIB_DIR)/include \
              -I$(PMLIB_DIR)/include \
              -I$(TP_DIR)/include \
              -I$(MPI_DIR)/include \
 							-I$(BCM_DIR)/include \
							-I$(BCM_FILEIO)/include \
							-I$(BCM_UTILS)/include 
             
UDEF_LIB_PATH=-L$(POLYLIB_DIR)/lib \
              -L$(CUTLIB_DIR)/lib \
              -L$(PMLIB_DIR)/lib \
              -L$(TP_DIR)/lib \
							-L$(BCM_DIR)/lib \
              -L$(MPI_DIR)/lib
              
UDEF_LIBS := -lbcm -lMPIPOLY -lcutlib -lPM -lMPITP -lTP -lxml2 -lmpi 
UDEF_LIBS := -lbcm -lMPIPolylib -lcutlib -lpmlib -lTextParser_mpi -lxml2 -lmpi 
UDEF_LIBS := -lbcm -lcutlib -lMPIPOLY -lPM -lMPITP -lmpi 
UDEF_LIBS := -lbcm -lcutlib -lMPIPolylib -lPM -lMPITP -lmpi 
UDEF_LIBS := -lbcm -lCUT -lMPIPOLY -lPM -lMPITP -lmpi 

#####

PROGRAM = ffv

FFV_SOLVERNAME	:= "FrontFlow/violet-HC (alpha)"
FFV_VERSION			:= "0.9.4"
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

vpath %.cpp $(BCM_FILEIO)/src
vpath %.cpp $(BCM_UTILS)/src

CSRC := $(wildcard *.cpp)
CHDR := $(wildcard *.h  )
FSRC := $(wildcard *.f90)
OBJS = $(CSRC:.cpp=.o) $(FSRC:.f90=.o)
BCM_UTILS_OBJS = ConfigFile.o ConfigBase.o BCMOctree.o BCMPolylib.o
BCM_UTILS_OBJS = BCMOctree.o BCMPolylib.o
BCM_FILEIO_OBJS = BitVoxel.o ErrorUtil.o BCMFileSaver.o LeafBlockSaver.o FileSystemUtil.o

PROGRAM: $(OBJS) $(BCM_UTILS_OBJS) $(BCM_FILEIO_OBJS)
	$(CXX) \
	$(LDFLAGS) \
	-o $(PROGRAM) $(OBJS) $(BCM_UTILS_OBJS) $(BCM_FILEIO_OBJS) \
	$(UDEF_LIB_PATH) $(UDEF_LIBS) \
	$(LIBS) 
	mv $(PROGRAM) ../../bin/

.cpp.o:
	$(CXX) $(CXXFLAGS) $(UDEF_INC_PATH) $(UDEF_OPT) -o $@ -c $<

.SUFFIXES: .f90
.f90.o:
	$(F90) $(F90FLAGS) $(UDEF_INC_PATH) $(UDEF_OPT) -o $@ -c $<

clean:
	rm -rf *.o ../../bin/$(PROGRAM) depend.inc

cleandata:
	rm -rf BIN STL VTK PLOT3D BCM_OUT CS data-*

depend:
	$(CXX) -MM $(CSRC) $(CXXFLAGS) $(UDEF_INC_PATH) $(UDEF_OPT) > depend.inc

-include depend.inc

