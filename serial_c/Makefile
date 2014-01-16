#------------------------------------------------------------------------------

INCLUDES = -I. -I../common_c
SOURCES = *.c ../common_c/*.c
EXECUTABLE = sweep

debug:
	gcc -o $(EXECUTABLE) $(SOURCES) $(INCLUDES)

release: release_gnu

release_gnu:
	gcc -o $(EXECUTABLE) $(SOURCES) $(INCLUDES) -DNDEBUG -O3 -fomit-frame-pointer -funroll-loops

release_pgi:
	pgcc -o $(EXECUTABLE) $(SOURCES) $(INCLUDES) -DNDEBUG -fast -Mfprelaxed -alias=ansi -Msafeptr -tp=istanbul-64

#---
#---Build for run on Cray backend node.
#---Do "module swap PrgEnv-pgi PrgEnv-gnu" to build this way.
#---
release_cray:
	cc -o $(EXECUTABLE) $(SOURCES) $(INCLUDES) -DUSE_MPI -DNDEBUG -O3 -fomit-frame-pointer -funroll-loops

#------------------------------------------------------------------------------

run:
	./$(EXECUTABLE) 10 10 10 30 16 64

runmany:
	bash -c 'for e in `seq 1 3` ; do for i in `seq 1 7` ; do for j in `seq 1 7` ; do ./$(EXECUTABLE) 10 10 10 $$e $$i $$j ; done; done; done'

test:
	./tester

run_cray:
	aprun -n 1 $(EXECUTABLE) 10 10 10 30 16 64

#------------------------------------------------------------------------------
