#!/bin/bash -l
#------------------------------------------------------------------------------

# CLEANUP
rm -rf CMakeCache.txt
rm -rf CMakeFiles

# SOURCE AND INSTALL
if [ "$SOURCE" = "" ] ; then
  SOURCE=../../
fi
if [ "$INSTALL" = "" ] ; then
  INSTALL=./install
fi

if [ "$BUILD" = "" ] ; then
  #BUILD=Debug
  BUILD=Release
fi

if [ "$NM_VALUE" = "" ] ; then
  NM_VALUE=4
fi

#------------------------------------------------------------------------------

MPI_INCLUDE_DIR=/usr/include/openmpi
MPI_LIB=/usr/lib/openmpi/lib/libmpi.so

cmake \
  -DCMAKE_BUILD_TYPE:STRING="$BUILD" \
  -DCMAKE_INSTALL_PREFIX:PATH="$INSTALL" \
 \
  -DCMAKE_C_COMPILER:STRING=pgcc \
\
  -DCMAKE_C_FLAGS:STRING="-O3 -DNM_VALUE=$NM_VALUE $ALG_OPTIONS -DUSE_ACC -I$MPI_INCLUDE_DIR -acc -Minfo=accel -Msafeptr -ta=tesla,cc70,cc60" \
  -DCMAKE_C_FLAGS_DEBUG:STRING="-g" \
  -DCMAKE_C_FLAGS_RELEASE:STRING="" \
\
  -DUSE_MPI:BOOL=ON \
  -DMPI_C_INCLUDE_PATH:STRING=$MPI_INCLUDE_DIR \
  -DMPI_C_LIBRARIES:STRING=$MPI_LIB \
  -DCMAKE_EXE_LINKER_FLAGS:STRING=$MPI_LIB \
 \
  $SOURCE

#------------------------------------------------------------------------------
