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
# GCC
########################

MPI_DIR			:= /usr/lib/openmpi
MPI_DIR			:= /usr/local

FFV_HOME		:= $(HOME)/FFV
TP_DIR      := $(FFV_HOME)/TextParser
PMLIB_DIR   := $(FFV_HOME)/PMlib
POLYLIB_DIR := $(FFV_HOME)/Polylib
CUTLIB_DIR  := $(FFV_HOME)/Cutlib
BCM_DIR     := $(FFV_HOME)/BCMTools
BCM_FILEIO  := $(BCM_DIR)/FileIO
BCM_UTILS   := $(BCM_DIR)/Utils
BCM_FILEIO  := /Users/jonishi/Devel/BCMTools_20130709/FileIO
BCM_UTILS   := /Users/jonishi/Devel/BCMTools_20130709/Utils

CC          := gcc
CXX         := g++
FC          := gfortran
F90         := gfortran

OPT_FLAGS   := -O3 -fopenmp
DEFINES     := 

CFLAGS      := $(OPT_FLAGS) $(DEFINES)
#CFLAGS      += -g -Wall -Wno-sign-compare -Wunknown-pragmas 
FCFLAGS     := $(OPT_FLAGS) $(DEFINES) -cpp -ffree-form 
CXXFLAGS    := $(CFLAGS)
F90FLAGS    := $(FCFLAGS)
LDFLAGS     := -fopenmp
LIBS        := -lgfortran -lmpi -lmpi_cxx -lmpi_mpifh
LIBS        := -lgfortran -lmpi -lmpi_cxx -lmpi_f77 -lmpi_f90
LIBS				+= -lz

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

