#!/bin/bash -l
#------------------------------------------------------------------------------

if [ "$PE_ENV" = "PGI" ] ; then
  module swap PrgEnv-pgi PrgEnv-gnu
fi
module load cudatoolkit

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

#------------------------------------------------------------------------------

#  -DCMAKE_C_FLAGS:STRING="-DNM_VALUE=$NM_VALUE -O3 -fomit-frame-pointer -funroll-loops -finline-limit=100000000" \

cmake \
  -DCMAKE_BUILD_TYPE:STRING="$BUILD" \
  -DCMAKE_INSTALL_PREFIX:PATH="$INSTALL" \
 \
  -DCMAKE_C_COMPILER:STRING="$(which cc)" \
  -DMPI_C_COMPILER="$(which cc)" \
  -DCMAKE_C_FLAGS:STRING="-DNM_VALUE=$NM_VALUE" \
 \
  -DUSE_MPI:BOOL=ON \
 \
  -DUSE_CUDA:BOOL=ON \
  -DCUDA_NVCC_FLAGS:STRING="-I$MPICH_DIR/include;-arch=sm_35;-O3;-use_fast_math;--maxrregcount;128;-DNDEBUG;-Xcompiler;-fstrict-aliasing;-Xcompiler;-fargument-noalias-global" \
  -DCUDA_HOST_COMPILER:STRING=/usr/bin/gcc \
  -DCUDA_PROPAGATE_HOST_FLAGS:BOOL=ON \
 \
  $SOURCE

#------------------------------------------------------------------------------

#  -DMPI_EXEC="aprun" \
#  -DMPI_EXEC_MAX_NUMPROCS:STRING=16 \
#  -DMPI_EXEC_NUMPROCS_FLAG:STRING=-n \
