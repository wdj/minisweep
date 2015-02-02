#==============================================================================

cmake_minimum_required(VERSION 2.8)

project("minisweep")

SET(INCLUDE_DIRS
  ${CMAKE_SOURCE_DIR}/driver
  ${CMAKE_SOURCE_DIR}/base
  ${CMAKE_SOURCE_DIR}/sn_base
  ${CMAKE_SOURCE_DIR}/sn
  )

SET(SOURCES
  driver/run_tools.c
  base/arguments.c
  base/memory.c
  sn_base/dimensions.c
  sn_base/pointer.c
  sn/faces_kba.c
  sn/quantities.c
  sn/stepscheduler_kba.c
  sn/sweeper.c
  sn/sweeper_kernels.c
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
  STRING(REPLACE " " ";" CMAKE_C_FLAGS_ADJUSTED ${CMAKE_C_FLAGS})
  SET(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};${CMAKE_C_FLAGS_ADJUSTED}")
  CUDA_INCLUDE_DIRECTORIES(${INCLUDE_DIRS})
  CUDA_ADD_LIBRARY(sweeper STATIC ${CUDA_SOURCES})
  CUDA_ADD_EXECUTABLE(sweep driver/sweep.cu)
  TARGET_LINK_LIBRARIES(sweep sweeper)
  CUDA_ADD_EXECUTABLE(tester driver/tester.cu)
  TARGET_LINK_LIBRARIES(tester sweeper)
ELSE()
  INCLUDE_DIRECTORIES(${INCLUDE_DIRS})
  ADD_LIBRARY(sweeper STATIC ${SOURCES})
  ADD_EXECUTABLE(sweep driver/sweep.c)
  TARGET_LINK_LIBRARIES(sweep sweeper)
  ADD_EXECUTABLE(tester driver/tester.c)
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
