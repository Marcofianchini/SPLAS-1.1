#
#    SPLAS v 0.1 MAKEFILE
#

########################################## COMPILER AND FLAGS

FC=gfortran
F2PY=f2py3.5

ifeq ($(FC),ifort)
F2PYFLAG = --fcompiler=intelem
F2PYFFLAG = "-free -fPIC -g"
FFLAG = -free -fPIC -g -cpp
DEBUG_FLAG= -g -O0
LIBFLAG= -qopenmp

else ifeq ($(FC),gfortran)

F2PYFLAG =
F2PYFFLAG = "-ffree-form -fPIC -Wall -g -ffree-line-length-none"
FFLAG = -ffree-form -fPIC -Wall -ffree-line-length-none -g -cpp
DEBUG_FLAG = -ffpe-trap=zero,invalid,overflow,underflow -fbacktrace -fbounds-check -g -O0
LIBFLAG=

endif

########################################## LIBRARY FLAGS

SPLAS_LIB=./
SPLAS_INC=./

#USE SYSTEM LIB
ifndef BLAS_LIB

BLAS_LIB=/usr/lib/x86_64-linux-gnu/
BLAS_INC=/usr/include/x86_64-linux-gnu/

endif

#MKL OR OTHER BLAS
ifeq ($(FC),ifort) #running on pico

LIB= -L${MKLROOT}/lib/intel64 -lmkl_rt -lmkl_core -lmkl_intel_thread -lpthread -lm -ldl
INC= -qopenmp -I${MKLROOT}/include
MKLSPLAS = 1

else ifeq ($(FC),gfortran)

LIB= -L$(BLAS_LIB) -lblas

ifdef BLAS_INC
INC= -I$(BLAS_INC)
endif

endif

########################################## LIB AND WRAPPER

OBJECTS=splas.o solver.o functions.o heuns_method.o

libsplas.so : $(OBJECTS)
		$(FC) -shared $(LIBFLAG) $^ $(LIB) -o $@

splaspy : splas_py_wrapper.f90
		$(F2PY) $(F2PYFLAG)  --f90flags=$(F2PYFFLAG) -L$(SPLAS_LIB) -I$(SPLAS_INC) -lsplas -c $^ -m $@


########################################## DEPENDENCIES

splas.o: functions.o solver.o

functions.o:

heuns_method.o:

solver.o: heuns_method.o functions.o

########################################## DEFAULTS AND PHONIES

default :
 	%.o : %.f90
	$(FC) $(FFLAG) $(DEBUG_FLAG) $(INC) -c $<  -o $@

.PHONY: all lib clean

all : libsplas.so splaspy

lib : libsplas.so

wrapper_py : splaspy

clean :
	rm *.so *.mod *.o
