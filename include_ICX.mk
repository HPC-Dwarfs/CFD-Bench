ifeq ($(ENABLE_MPI),true)
CC = mpiicx
DEFINES  = -D_MPI
else
CC = icx
endif

GCC = gcc
LINKER = $(CC)

ifeq ($(ENABLE_OPENMP),true)
OPENMP   = -qopenmp
endif

VERSION  = --version
CFLAGS   =  -O3 -xHost -qopt-zmm-usage=high -std=c99 $(OPENMP) -Wno-unused-command-line-argument
LFLAGS   = $(OPENMP)
DEFINES  += -D_GNU_SOURCE# -DDEBUG
INCLUDES =
LIBS     =
