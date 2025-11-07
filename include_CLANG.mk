ifeq ($(ENABLE_MPI),true)
CC   = mpicc
DEFINES  = -D_MPI
else
CC   = cc
endif

GCC  = cc
LINKER = $(CC)

ifeq ($(ENABLE_OPENMP),true)
OPENMP   = -fopenmp
#OPENMP   = -Xpreprocessor -fopenmp #required on Macos with homebrew libomp
LIBS     = # -lomp
endif

VERSION  = --version
CFLAGS   = -Ofast -std=c17
LFLAGS   = $(OPENMP) -lm
DEFINES  += -D_GNU_SOURCE# -DDEBUG
INCLUDES = -I/opt/homebrew/include
