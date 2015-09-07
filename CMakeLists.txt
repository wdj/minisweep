#==============================================================================

cmake_minimum_required(VERSION 2.8)

project("minisweep")

SET(INCLUDE_DIRS
  ${CMAKE_SOURCE_DIR}/src_01_base
  ${CMAKE_SOURCE_DIR}/src_02_sweeper_base
  ${CMAKE_SOURCE_DIR}/src_03_sweeper
  ${CMAKE_SOURCE_DIR}/src_04_driver
  )

SET(SOURCES
  src_01_base/arguments.c
  src_02_sweeper_base/dimensions.c
  src_02_sweeper_base/pointer.c
  src_03_sweeper/faces_kba.c
  src_03_sweeper/quantities.c
  src_03_sweeper/stepscheduler_kba.c
  src_03_sweeper/sweeper.c
  src_03_sweeper/sweeper_kernels.c
  src_04_driver/run_tools.c
  )

SET(CUDA_SOURCES)
FOREACH(FILE IN LISTS SOURCES)
  STRING(REPLACE ".c" ".cu" FILE2 ${FILE})
  SET(CUDA_SOURCES ${CUDA_SOURCES} ${FILE2})
ENDFOREACH()

IF(USE_MPI)
  find_package(MPI REQUIRED)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_MPI")
ENDIF()

IF(USE_CUDA)
  find_package(CUDA REQUIRED)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_CUDA")
  STRING(REPLACE " " ";" CMAKE_C_FLAGS_ADJUSTED ${CMAKE_C_FLAGS})
  SET(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};${CMAKE_C_FLAGS_ADJUSTED}")
#  STRING(REGEX MATCH " -DNM_VALUE=[0-9]*" NM_VALUE_DEF_ " ${CMAKE_C_FLAGS} ")
#  STRING(REPLACE " " ";" NM_VALUE_DEF "${NM_VALUE_DEF_}")
#  SET(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS}${NM_VALUE_DEF};-DUSE_CUDA")
  CUDA_INCLUDE_DIRECTORIES(${INCLUDE_DIRS})
  CUDA_ADD_LIBRARY(sweeper STATIC ${CUDA_SOURCES})
  CUDA_ADD_EXECUTABLE(sweep src_04_driver/sweep.cu)
  TARGET_LINK_LIBRARIES(sweep sweeper)
  CUDA_ADD_EXECUTABLE(tester src_04_driver/tester.cu)
  TARGET_LINK_LIBRARIES(tester sweeper)
ELSE()
  INCLUDE_DIRECTORIES(${INCLUDE_DIRS})
  ADD_LIBRARY(sweeper STATIC ${SOURCES})
  ADD_EXECUTABLE(sweep src_04_driver/sweep.c)
  TARGET_LINK_LIBRARIES(sweep sweeper)
  ADD_EXECUTABLE(tester src_04_driver/tester.c)
  TARGET_LINK_LIBRARIES(tester sweeper)
ENDIF()

install(TARGETS sweep DESTINATION bin)
#install(TARGETS tester DESTINATION bin)

SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)

IF(CMAKE_C_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ansi -Wmissing-prototypes -Dinline=")
ENDIF()

IF(CMAKE_C_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Dinline=__forceinline -D__restrict__=restrict -restrict")
ENDIF()

MESSAGE("Compiler ID: ${CMAKE_C_COMPILER_ID}")

# enable_testing()
# add_test(test1 tester)

#==============================================================================
