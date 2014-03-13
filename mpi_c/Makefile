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

ifndef ALG_OPTIONS
  ALG_OPTIONS = -DSWEEPER_KBA
endif

#==============================================================================

INCLUDES = -I. -I../common_c
SOURCES = *.c ../common_c/*.c
EXECUTABLE = sweep.x

CC =
CFLAGS :=

ifeq ($(COMPILER_OPTION),GNU)
  #---"module swap PrgEnv-pgi PrgEnv-gnu"
  CC = gcc
  ifeq ($(RELEASE_OPTION),DEBUG)
    CFLAGS := $(CFLAGS) -ansi -Dinline= -Wmissing-prototypes -g
  else
    CFLAGS := $(CFLAGS) -DNDEBUG -O3 -fomit-frame-pointer -funroll-loops \
                       -finline-limit=100000000
  endif
endif

ifeq ($(COMPILER_OPTION),PGI)
  CC = pgcc
  CFLAGS := $(CFLAGS) -DNDEBUG -fast -Mfprelaxed -alias=ansi -Msafeptr -tp=istanbul-64
endif

ifeq ($(MPI_OPTION),MPI)
  CFLAGS := $(CFLAGS) -DUSE_MPI -I$(MPICH_DIR)/include -L$(MPICH_DIR)/lib \
                      -lmpich
endif

#==============================================================================

default:
	$(CC) -o $(EXECUTABLE) $(SOURCES) $(INCLUDES) $(CFLAGS) $(ALG_OPTIONS)

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
