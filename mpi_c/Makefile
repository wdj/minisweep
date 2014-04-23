#==============================================================================
# Makefile for minisweep.
#==============================================================================

#==============================================================================
# Default build options.
#==============================================================================

ifndef COMPILER_OPTION
  COMPILER_OPTION = GNU
endif

ifndef RELEASE_OPTION
  RELEASE_OPTION = DEBUG
endif

ifndef MPI_OPTION
  MPI_OPTION = MPI
endif

ifndef OPENMP_OPTION
  OPENMP_OPTION =
endif

ifndef CUDA_OPTION
  CUDA_OPTION =
endif

ifndef ALG_OPTIONS
  ALG_OPTIONS = -DSWEEPER_KBA
endif

#==============================================================================
# Settings.
#==============================================================================

INCLUDES = -I. -I../common_c
SOURCES      = $(wildcard *.c  ../common_c/*.c  )
CUDA_SOURCES = $(wildcard *.cu ../common_c/*.cu )
EXECUTABLE = sweep.x

default: make_executable

CC =
CFLAGS :=

ifeq ($(COMPILER_OPTION),GNU)
  #---"module swap PrgEnv-pgi PrgEnv-gnu"
  CC = cc
  ifeq ($(RELEASE_OPTION),DEBUG)
    CFLAGS := $(CFLAGS) -ansi -Dinline= -Wmissing-prototypes -g
  else
    CFLAGS := $(CFLAGS) -DNDEBUG -O3 -fomit-frame-pointer -funroll-loops \
                       -finline-limit=100000000
  endif
endif

ifeq ($(COMPILER_OPTION),PGI)
  CC = cc
  CFLAGS := $(CFLAGS) -DNDEBUG -fast -Mfprelaxed -alias=ansi -Msafeptr -tp=istanbul-64
endif

ifeq ($(MPI_OPTION),MPI)
  CFLAGS := $(CFLAGS) -DUSE_MPI -I$(MPICH_DIR)/include -L$(MPICH_DIR)/lib \
                      -Wl,-rpath=$(MPICH_DIR)/lib -lmpich
endif

ifneq ($(OPENMP_OPTION),)
  CFLAGS := $(CFLAGS) -DUSE_OPENMP
  ifeq ($(COMPILER_OPTION),PGI)
    CFLAGS := $(CFLAGS) -mp=bind
  endif
  ifeq ($(COMPILER_OPTION),GNU)
    CFLAGS := $(CFLAGS) -fopenmp
  endif
  COMMA := ,
  ifneq ($(findstring $(COMMA)THREADS$(COMMA),$(COMMA)$(OPENMP_OPTION)$(COMMA)),)
    CFLAGS := $(CFLAGS) -DUSE_OPENMP_THREADS
  endif
  ifneq ($(findstring $(COMMA)VO_ATOMIC$(COMMA),$(COMMA)$(OPENMP_OPTION)$(COMMA)),)
    CFLAGS := $(CFLAGS) -DUSE_OPENMP_VO_ATOMIC
  endif
endif

ifneq ($(CUDA_OPTION),)
  ifeq ($(COMPILER_OPTION),PGI)
    $(error "CUDA requires GNU compiler to be used.")
  endif
make_cu: make_c
	nvcc -c $(CUDA_SOURCES) $(INCLUDES) $(ALG_OPTIONS)
else
make_cu: make_c
endif

#==============================================================================
# Targets.
#==============================================================================

make_c:
	$(CC) -c $(SOURCES) $(INCLUDES) $(CFLAGS) $(ALG_OPTIONS)

make_executable: make_cu
	$(CC) -o $(EXECUTABLE) *.o $(INCLUDES) $(CFLAGS) $(ALG_OPTIONS)
	rm *.o

test:
	./test.bash

#==============================================================================

#run:
#	./$(EXECUTABLE) --nx 10 --ny 10 --nz 10 --ne 30 --nm 16 --na 33
#
#runmany:
#	bash -c 'for e in `seq 1 3` ; do for i in `seq 1 7` ; do for j in `seq 1 7` ; do ./$(EXECUTABLE) 10 10 10 $$e $$i $$j ; done; done; done'
#
#test:
#	./tester
#
#run_cray:
#	aprun -n 1 $(EXECUTABLE) --nx 10 --ny 10 --nz 10 --ne 30 --nm 16 --na 64

distclean:
	rm $(EXECUTABLE)

#==============================================================================
