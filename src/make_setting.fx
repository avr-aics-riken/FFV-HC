#############################################################################
#
# FFV-HC
# FrontFlow/violet for hierarchical Cartesian grids
#
# Copyright (c) 2013 Institute of Industrial Science, The University of Tokyo. 
# All rights reserved.
#
# Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN. 
# All rights reserved.
#
#############################################################################

########################
# FX10
########################

FFV_HOME		:= $(HOME)/FFV
TP_DIR      := $(FFV_HOME)/TextParser
PMLIB_DIR   := $(FFV_HOME)/PMlib
POLYLIB_DIR := $(FFV_HOME)/Polylib
CUTLIB_DIR  := $(FFV_HOME)/Cutlib
BCM_DIR     := $(FFV_HOME)/BCMTools
BCM_FILEIO  := $(BCM_DIR)/FileIO
BCM_UTILS   := $(BCM_DIR)/Utils

CC          := mpifccpx
CXX         := mpiFCCpx
FC          := mpifrtpx
F90         := mpifrtpx

OPT_FLAGS   := -Kfast,ocl,preex,simd=2,uxsimd,array_private,openmp,parallel,optmsg=1
OPT_FLAGS_F := -Kfast,ocl,preex,simd=2,uxsimd,array_private,openmp,parallel,auto,optmsg=1 
DEFINES     := 
XG          := -Xg

CFLAGS      := $(OPT_FLAGS)   $(DEFINES) $(XG) -V -Nsrc -x0
FCFLAGS     := $(OPT_FLAGS_F) $(DEFINES) $(XG) -Cpp -Qt -x-
CXXFLAGS    := $(CFLAGS)
F90FLAGS    := $(FCFLAGS)
LDFLAGS     := --linkfortran $(OPT_FLAGS)
LIBS        :=

## iff double
CFLAGS     += -D_REAL_IS_DOUBLE_
CXXFLAGS   += -D_REAL_IS_DOUBLE_
FCFLAGS    += -D_REAL_IS_DOUBLE_ -fdefault-real-8
F90FLAGS   += -D_REAL_IS_DOUBLE_ -fdefault-real-8

## iff large block
#CFLAGS     += -D_BLOCK_IS_LARGE_
#CXXFLAGS   += -D_BLOCK_IS_LARGE_
#FCFLAGS    += -D_BLOCK_IS_LARGE_
#F90FLAGS   += -D_BLOCK_IS_LARGE_

## iff small block
CFLAGS     += -D_BLOCK_IS_SMALL_
CXXFLAGS   += -D_BLOCK_IS_SMALL_
FCFLAGS    += -D_BLOCK_IS_SMALL_
F90FLAGS   += -D_BLOCK_IS_SMALL_

## iff different restart with staging
#CFLAGS     += -D_STAGING_
#CXXFLAGS   += -D_STAGING_
#FCFLAGS    += 
#F90FLAGS   += 

