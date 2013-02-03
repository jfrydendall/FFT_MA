.SUFFIXES: .o .cc .f90 .f03 .cpp .h .hpp  
FC = gfortran
LD  = gfortran

CC = g++
LCC = g++

AR = ar rv
OPT_FLAGS = -m64 -std=c++11 -DNDEBUG -fopenmp -g -Wall -ansi -O3 -fomit-frame-pointer -fstrict-aliasing -ffast-math -msse2 -mfpmath=sse -march=native 

FOPT_FLAGS = -O3 -fprefetch-loop-arrays -funroll-loops -floop-interchange -ftree-loop-distribution -ffast-math -funroll-loops -finline-functions -ftree-vectorize -march=native -fno-signed-zeros -ffinite-math-only -fopenmp

LINKFLAGS = $(OPT_FLAGS)

FLINKFLAGS = $(FOPT_FLAGS)

LIBS = -lm -lgomp -lRmath -lfftw3 -lfftw3_omp
INCLUDE	= -I../fftw++

.cpp.o:
	$(CC) $(OPT_FLAGS) $(INCLUDE) -o $*.o -c $*.cpp

.cc.o:
	$(CC) $(OPT_FLAGS) $(INCLUDE) -o $*.o -c $*.cc

.f03.o:
	$(FC) $(FOPT_FLAGS) -o $*.o -c $*.f03

.f90.o:
	$(FC) $(FOPT_FLAGS) -o $*.o -c $*.f90

include source.files
include target.mk

OBJECTS = $(CPPFILES:.cpp=.o) $(CCFILES:.cc=.o)

F90_FILES :=  $(filter %.f90, $(FORTRANFILES))
F03_FILES := $(filter %.f03, $(FORTRANFILES))

FOBJECTS = $(F90_FILES:.f90=.o) $(F03_FILES:.f03=.o) 

all:	$(TARGET) $(FORTRAN)

$(TARGET): $(OBJECTS)
	$(LCC) $(LINKFLAGS) $(INCLUDE) -o $(TARGET) $(OBJECTS) $(LIBS)

$(FORTRAN): $(FOBJECTS)
	$(LD) $(FLINKFLAGS) -o $(FORTRAN) $(FOBJECTS) $(LIBS)
clean:
	rm -f $(TARGET)
	rm -rf *.o *.mod *.so *.a
