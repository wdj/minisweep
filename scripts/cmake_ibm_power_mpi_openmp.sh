#!/bin/bash -l
#------------------------------------------------------------------------------

# CLEANUP
rm -rf CMakeCache.txt
rm -rf CMakeFiles

# SOURCE AND INSTALL
if [ "$SOURCE" = "" ] ; then
  SOURCE=../minisweep
fi
if [ "$INSTALL" = "" ] ; then
  INSTALL=../install
fi

if [ "$BUILD" = "" ] ; then
  BUILD=Debug
  #BUILD=Release
fi

if [ "$NM_VALUE" = "" ] ; then
  NM_VALUE=4
fi

CC=gcc
OMP_ARGS="-fopenmp"
OPT_ARGS="-O3 -fomit-frame-pointer -funroll-loops -finline-limit=10000000"

#------------------------------------------------------------------------------

# See also `which mpcc`
MPI_INCLUDE_DIR=/opt/ibmhpc/pecurrent/mpich/gnu/include64
MPI_LIB=/opt/ibmhpc/pecurrent/mpich/gnu/lib64/libmpi.so

cmake \
  -DCMAKE_BUILD_TYPE:STRING="$BUILD" \
  -DCMAKE_INSTALL_PREFIX:PATH="$INSTALL" \
 \
  -DCMAKE_C_COMPILER:STRING="gcc" \
  -DCMAKE_C_FLAGS:STRING="-DNM_VALUE=$NM_VALUE -I$MPI_INCLUDE_DIR -DUSE_OPENMP -DUSE_OPENMP_THREADS $OMP_ARGS" \
  -DCMAKE_C_FLAGS_DEBUG:STRING="-g" \
  -DCMAKE_C_FLAGS_RELEASE:STRING="$OPT_ARGS" \
 \
  -DUSE_MPI:BOOL=ON \
  -DMPI_C_INCLUDE_PATH:STRING=$MPI_INCLUDE_DIR \
  -DMPI_C_LIBRARIES:STRING=$MPI_LIB \
 \
  -DCMAKE_EXE_LINKER_FLAGS:STRING=$MPI_LIB \
 \
  $SOURCE

#------------------------------------------------------------------------------
