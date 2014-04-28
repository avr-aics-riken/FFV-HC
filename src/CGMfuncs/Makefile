-include ../make_setting

INC_PATH := -I$(MPI_DIR)/include
LIB_PATH := -L$(MPI_DIR)/lib
LIBS     := -lmpi -lmpi_cxx -lmpi_f90 -lmpi_f77

TARGET = libcgmfuncs.a

CSRC = $(wildcard *.cpp)
CHDR = $(wildcard *.h)
FSRC = $(wildcard *.f90)
OBJS = $(CSRC:.cpp=.o) $(FSRC:.f90=.o)

all:		$(TARGET)

$(TARGET):	$(OBJS)
	$(AR) r $@ $(OBJS) 
	$(RANLIB) $@

.SUFFIXES: .cpp
.cpp.o:
	$(CXX) $(CXXFLAGS) $(INC_APTH) -c $<

.SUFFIXES: .f90
.f90.o:
	$(F90) $(F90FLAGS) $(INC_PATH) -c $<

.PHONY: clean cleandata

clean:
	rm -rf *.o *.s *.mod $(TARGET) *.dep

depend:
	$(CXX) $(XG) -MM $(CSRC) $(CXXFLAGS) > Makefile.dep

-include Makefile.dep

