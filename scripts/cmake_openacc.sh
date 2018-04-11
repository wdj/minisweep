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

cmake \
  -DCMAKE_BUILD_TYPE:STRING="$BUILD" \
  -DCMAKE_INSTALL_PREFIX:PATH="$INSTALL" \
 \
  -DCMAKE_C_COMPILER:STRING=pgcc \
  -DCMAKE_CXX_COMPILER:STRING=pgc++ \
  -DCMAKE_C_FLAGS:STRING="-O3 -DNM_VALUE=$NM_VALUE $ALG_OPTIONS -DUSE_ACC=ON -acc -Minfo=accel -Msafeptr -ta=tesla,cc35" \
\
  -DCMAKE_C_FLAGS_DEBUG:STRING="-g" \
  -DCMAKE_C_FLAGS_RELEASE:STRING="" \
  -DCMAKE_CXX_FLAGS:STRING="-O3 -DNM_VALUE=$NM_VALUE $ALG_OPTIONS _DUSE_ACC=ON -acc -Minfo=accel -Msafeptr -ta=tesla,cc35" \
\
  -DCMAKE_CXX_FLAGS_DEBUG:STRING="-g" \
  -DCMAKE_CXX_FLAGS_RELEASE:STRING="" \
 \
  $SOURCE

#------------------------------------------------------------------------------
