#!/bin/bash -l
#------------------------------------------------------------------------------

# CLEANUP
rm -rf CMakeCache.txt
rm -rf CMakeFiles

if [ "$COMPILER" = "xlC" ] ; then
  module load xlC
  CC="xlc_r"
  CXX="xlC_r"
  OMP_ARGS="-qsmp"
  OPT_ARGS=""
elif [ "$COMPILER" = "pgi" ] ; then
  module load pgi
  CC="pgcc"
  CXX="pgc++"
  OMP_ARGS="-mp"
  OPT_ARGS=""
elif [ "$COMPILER" = "clang" ] ; then
  module load clang
  module load cuda
  CC="clang"
  CXX="clang++"
  OMP_ARGS="-fopenmp -I$CLANG_OMP_INCLUDE -L$CLANG_OMP_LIB"
  OPT_ARGS=""
else
  CC="gcc"
  CXX="g++"
  OMP_ARGS="-fopenmp"
  OPT_ARGS="-O3 -fomit-frame-pointer -funroll-loops -finline-limit=10000000"
fi

if [ "$NO_OPENMP" != "" ] ; then
  OMP_ARGS=""
else
  OMP_ARGS="$OMP_ARGS -DUSE_OPENMP -DUSE_OPENMP_THREADS"
fi


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

# See also `which mpcc`
MPI_INCLUDE_DIR=/opt/ibmhpc/pecurrent/mpich/gnu/include64
MPI_LIB=/opt/ibmhpc/pecurrent/mpich/gnu/lib64/libmpi.so

cmake \
  -DCMAKE_BUILD_TYPE:STRING="$BUILD" \
  -DCMAKE_INSTALL_PREFIX:PATH="$INSTALL" \
 \
  -DCMAKE_C_COMPILER:STRING="$CC" \
  -DCMAKE_CXX_COMPILER:STRING="$CXX" \
  -DCMAKE_C_FLAGS:STRING="-DNM_VALUE=$NM_VALUE -I$MPI_INCLUDE_DIR $OMP_ARGS" \
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

#  -DMPI_C_COMPILER:STRING=xlc_r \

