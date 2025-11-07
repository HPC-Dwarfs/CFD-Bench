ifeq ($(ENABLE_MPI),true)
CC   = mpicc
DEFINES  = -D_MPI
else
CC   = gcc
endif

LD = $(CC)

ifeq ($(ENABLE_OPENMP),true)
OPENMP   = -fopenmp
endif

VERSION  = --version
CFLAGS   = -O3 -ffast-math -std=c99 $(OPENMP)
LFLAGS   = $(OPENMP)
DEFINES  += -D_GNU_SOURCE
INCLUDES =
LIBS     =
