#!/bin/bash -l
#------------------------------------------------------------------------------

module load cmake3
module load cuda

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

cmake \
  -DCMAKE_BUILD_TYPE:STRING="$BUILD" \
  -DCMAKE_INSTALL_PREFIX:PATH="$INSTALL" \
 \
  -DCMAKE_C_COMPILER:STRING="gcc" \
  -DCMAKE_C_FLAGS:STRING="-DNM_VALUE=$NM_VALUE" \
 \
  -DUSE_CUDA:BOOL=ON \
  -DCUDA_NVCC_FLAGS:STRING="-I$MPICH_DIR/include;-arch=sm_35;-O3;-use_fast_math;-DNDEBUG;--maxrregcount;128;-Xcompiler;-fstrict-aliasing;-Xcompiler;-fargument-noalias-global;-Xcompiler;-O3;-Xcompiler;-fomit-frame-pointer;-Xcompiler;-funroll-loops;-Xcompiler;-finline-limit=100000000;-Xptxas=-v" \
  -DCUDA_HOST_COMPILER:STRING=/usr/bin/gcc \
  -DCUDA_PROPAGATE_HOST_FLAGS:BOOL=ON \
 \
  $SOURCE

#------------------------------------------------------------------------------

#  -DMPI_EXEC="aprun" \
#  -DMPI_EXEC_MAX_NUMPROCS:STRING=16 \
#  -DMPI_EXEC_NUMPROCS_FLAG:STRING=-n \
