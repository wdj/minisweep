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
COMMON_FLAGS := 

ifeq ($(COMPILER_OPTION),GNU)
  #---"module swap PrgEnv-pgi PrgEnv-gnu"
  CC = cc
  ifeq ($(RELEASE_OPTION),DEBUG)
    CFLAGS := $(CFLAGS) -ansi -Wmissing-prototypes -g -Dinline=
    COMMON_FLAGS := $(COMMON_FLAGS)
  else
    CFLAGS := $(CFLAGS) -O3 -fomit-frame-pointer -funroll-loops \
                       -finline-limit=100000000
    COMMON_FLAGS := $(COMMON_FLAGS) -DNDEBUG
  endif
endif

ifeq ($(COMPILER_OPTION),INTELMIC)
  CC = mpiicc
  CFLAGS := $(CFLAGS) -mmic -opt-report=0 -vec-report -Dinline=__inline -finline -finline-functions -finline-limit=100 -ip -fno-alias -prec-div
  ifeq ($(RELEASE_OPTION),DEBUG)
    CFLAGS := $(CFLAGS)
    COMMON_FLAGS := $(COMMON_FLAGS)
                   # -strict-ansi -Dinline=
  else
    CFLAGS := $(CFLAGS)
    COMMON_FLAGS := $(COMMON_FLAGS) -DNDEBUG
  endif
endif

ifeq ($(COMPILER_OPTION),PGI)
  CC = cc
  CFLAGS := $(CFLAGS) -fast -Mfprelaxed -alias=ansi -Msafeptr -tp=istanbul-64
  COMMON_FLAGS := $(COMMON_FLAGS) -DNDEBUG
endif

ifeq ($(MPI_OPTION),MPI)
  ifneq ($(COMPILER_OPTION),INTELMIC)
    CFLAGS := $(CFLAGS) -L$(MPICH_DIR)/lib -Wl,-rpath=$(MPICH_DIR)/lib -lmpich
    COMMON_FLAGS := $(COMMON_FLAGS) -DUSE_MPI -I$(MPICH_DIR)/include
  else
    CFLAGS := $(CFLAGS)
    COMMON_FLAGS := $(COMMON_FLAGS) -DUSE_MPI
  endif
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

ifneq ($(NM_VALUE),)
  COMMON_FLAGS := $(COMMON_FLAGS) -DNM_VALUE=$(NM_VALUE)
endif

ifeq ($(RELAXED_TESTING),1)
  COMMON_FLAGS := $(COMMON_FLAGS) -DRELAXED_TESTING
endif

ifneq ($(CUDA_OPTION),)
  CUDA_ARGS := -arch=sm_35 -O3 -use_fast_math --maxrregcount 128 -Xcompiler "-fstrict-aliasing" -Xcompiler "-fargument-noalias-global"
  ifneq ($(RELEASE_OPTION),DEBUG)
    CUDA_ARGS := $(CUDA_ARGS) -DNDEBUG
  endif
  ifeq ($(COMPILER_OPTION),PGI)
    $(error "CUDA requires GNU compiler to be used.")
  endif
make_cu: make_c
	nvcc -c $(CUDA_SOURCES) $(COMMON_FLAGS) $(INCLUDES) $(ALG_OPTIONS) -DUSE_CUDA $(CUDA_ARGS)
else
make_cu: make_c
endif

#==============================================================================
# Targets.
#==============================================================================

make_c:
	$(CC) -c $(SOURCES) $(COMMON_FLAGS) $(INCLUDES) $(CFLAGS) $(ALG_OPTIONS)

make_executable: make_cu
	$(CC) -o $(EXECUTABLE) *.o $(COMMON_FLAGS) $(INCLUDES) $(CFLAGS) $(ALG_OPTIONS)
	rm *.o

test:
	./test.bash

#==============================================================================

#make COMPILER_OPTION=INTELMIC NM_VALUE=16 RELEASE_OPTION=
#
#if [ $? = 0 ] ; then
#  \cp $HOME/minisweep/mpi_c/sweep.x $TMPDIR/mic0
#  ARGS="--nx 1 --ny 1 --nz 1 --ne 1 --na 32 --nproc_x 1 --nproc_y 1"
#  micmpiexec -n 1 -wdir $TMPDIR -host ${HOSTNAME}-mic0 $TMPDIR/sweep.x $ARGS
#fi

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
