#!/bin/bash -l
#------------------------------------------------------------------------------

if [ "$PE_ENV" = "PGI" ] ; then
  module swap PrgEnv-pgi PrgEnv-gnu
fi

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
  -DCMAKE_C_COMPILER:STRING="$(which cc)" \
  -DMPI_C_COMPILER="$(which cc)" \
  -DCMAKE_C_FLAGS:STRING="-DNM_VALUE=$NM_VALUE" \
  -DCMAKE_C_FLAGS_DEBUG:STRING="-g" \
  -DCMAKE_C_FLAGS_RELEASE:STRING="-O3 -fomit-frame-pointer -funroll-loops -finline-limit=10000000" \
 \
  -DUSE_MPI:BOOL=ON \
 \
  $SOURCE

#------------------------------------------------------------------------------

#  -DMPI_EXEC="aprun" \
#  -DMPI_EXEC_MAX_NUMPROCS:STRING=16 \
#  -DMPI_EXEC_NUMPROCS_FLAG:STRING=-n \
