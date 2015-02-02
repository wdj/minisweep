#!/bin/bash -l
#------------------------------------------------------------------------------

# CLEANUP
rm -rf CMakeCache.txt
rm -rf CMakeFiles

# SOURCE AND INSTALL
SOURCE=../minisweep
INSTALL=../install

BUILD=Debug
#BUILD=Release

OPT_ARGS="-mmic -vec-report1 -ip -prec-div -O3 -align -ansi-alias -fargument-noalias -fno-alias -fargument-noalias"

#------------------------------------------------------------------------------

cmake \
  -DCMAKE_BUILD_TYPE:STRING="$BUILD" \
  -DCMAKE_INSTALL_PREFIX:PATH="$INSTALL" \
  -DCMAKE_C_COMPILER:STRING=mpiicc \
  -DCMAKE_C_FLAGS:STRING="-DNM_VALUE=4 -DUSE_OPENMP -DUSE_OPENMP_THREADS -openmp -DUSE_MPI" \
  -DCMAKE_C_FLAGS_DEBUG:STRING="-mmic -g $OPT_ARGS" \
  -DCMAKE_C_FLAGS_RELEASE:STRING="-mmic $OPT_ARGS" \
 \
  -DUSE_MPI:BOOL=ON \
 \
  $SOURCE

#------------------------------------------------------------------------------
