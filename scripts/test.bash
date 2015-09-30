#!/bin/bash -l
#==============================================================================
cat <<EOF >/dev/null
===============================================================================

Run tests on sweep miniapp.

This script builds multiple versions of the executable and executes
a battery of tests on each one.

Usage:

cd minisweep/scripts
bash ./test.bash

Environment variables:

VERBOSE - set to 1 for verbose build output.

NOTE: to run properly on Titan, this must be executed in a PBS job and
must be run from Lustre.

===============================================================================
EOF
#==============================================================================

set -eu

declare g_ntest=0
declare g_ntest_passed=0
declare g_verbose=1

#==============================================================================
# Run the executable once.
#==============================================================================
function perform_runs_impl
{
  local exec_config_args="$(echo "$1" | tr '\n' ' ')"
  local application_args="$(echo "$2" | tr '\n' ' ')"
  local environment="${3:-}"
  local executable="$4"
  if [ "$environment" = "" ] ; then
    local env_environment=""
  else
    local env_environment="env $environment"
  fi

  if [ "${PBS_NP:-}" != "" -a "${CRAY_MPICH2_ROOTDIR:-}" != "" ] ; then
    #---If running on Cray back-end node, must cd to Lustre to do the aprun.
    local wd="$PWD"
    pushd "$MEMBERWORK" >/dev/null
    pwd
    echo $env_environment aprun $exec_config_args "$wd/$executable"
    $env_environment aprun $exec_config_args "$wd/$executable"
    #assert $? = 0
    popd >/dev/null
  elif [ "${IMICROOT:-}" != "" ] ; then
    cp $executable $TMPDIR/mic0
    $env_environment micmpiexec -n 1 -wdir $TMPDIR -host ${HOSTNAME}-mic0 $TMPDIR/$executable $application_arg
    rm $executable $TMPDIR/mic0
  else
    #assert $exec_config_args = "-n1"
    $env_environment ./$executable
    #assert $? = 0
  fi
}
#==============================================================================

#==============================================================================
# Run the executable multiple times using stdin.
#==============================================================================
function perform_runs
{
  local executable
  for executable in "tester" "sweep" ; do
    perform_runs_impl "$1" "$2" "${3:-}" "$executable" | tee tmp_

    local ntest=$(       grep TESTS tmp_ | sed -e 's/.*TESTS *//'  -e 's/ .*//')
    local ntest_passed=$(grep TESTS tmp_ | sed -e 's/.*PASSED *//' -e 's/ .*//')
    rm -f tmp_

    g_ntest=$((        $g_ntest        + $ntest ))
    g_ntest_passed=$(( $g_ntest_passed + $ntest_passed ))
  done
}
#==============================================================================

#==============================================================================
# Initialize build/execution environment.
#==============================================================================
function initialize
{
  if [ "${CRAY_MPICH2_ROOTDIR:-}" != "" ] ; then
    if [ "${PE_ENV:-}" != "GNU" ] ; then
      module swap PrgEnv-pgi PrgEnv-gnu
    fi
    module load cmake

    if [ "${PE_ENV:-}" != "GNU" ] ; then
      echo "Error: GNU compiler required." 1>&2
      exit 1
    fi
  fi
}
#==============================================================================

#==============================================================================
# Output pass/fail results.
#==============================================================================

function print_results_summary
{
  local is_final="${1:-}"

  local ntest_failed=$(( $g_ntest - $g_ntest_passed ))

  echo
  if [ "$is_final" = "" ] ; then
    echo -n "RUNNING TOTAL: "
  else
    echo -n "FINAL TOTAL: "
  fi
  echo -n "TESTS $g_ntest"
  echo -n "    "
  echo -n "PASSED $g_ntest_passed"
  echo -n "    "
  echo -n "FAILED $(( $g_ntest - $g_ntest_passed ))"
  echo -n "."
  if [ "$is_final" != "" ] ; then
    if [ $ntest_failed = 0 ] ; then
      echo -n "          ========== SUCCESS =========="
    fi
  fi
  echo
  echo
}
#==============================================================================


#==============================================================================
function main
{
  initialize

  #---args to use below:
  #---  ncell_x ncell_y ncell_z ne nm na numiterations nproc_x nproc_y nblock_z 

  if [ "${VERBOSE:-}" != "" ] ; then
    VERBOSE_ARG="VERBOSE=$VERBOSE"
  else
    VERBOSE_ARG=""
  fi

  #BUILD_DIR=../../build_test
  BUILD_DIR=/ccs/home/$(whoami)/atlas/minisweep_work/build_test
  SCRIPTS_DIR=$PWD/../scripts

  if [ "${CRAY_MPICH2_ROOTDIR:-}" != "" ] ; then
    local IS_TITAN=1
  else
    local IS_TITAN=0
  fi

  rm -f tmp_

  #==============================
  # Serial
  #==============================

  echo "--------------------------------------------------------"
  echo "---Serial tests---"
  echo "--------------------------------------------------------"

  #make MPI_OPTION= NM_VALUE=4

  rm -rf $BUILD_DIR
  mkdir $BUILD_DIR
  pushd $BUILD_DIR

  env NM_VALUE=4 $SCRIPTS_DIR/cmake_serial.sh
  make $VERBOSE_ARG

  perform_runs "-n1" ""

  popd
  rm -rf $BUILD_DIR

  print_results_summary

  #==============================
  # OpenMP
  #==============================

  echo "--------------------------------------------------------"
  echo "---OpenMP tests---"
  echo "--------------------------------------------------------"

  #make MPI_OPTION= OPENMP_OPTION=THREADS NM_VALUE=4

  rm -rf $BUILD_DIR
  mkdir $BUILD_DIR
  pushd $BUILD_DIR

  env NM_VALUE=4 $SCRIPTS_DIR/cmake_openmp.sh
  make $VERBOSE_ARG

  perform_runs "-n1 -d16" ""

  popd
  rm -rf $BUILD_DIR

  print_results_summary

  #==============================
  # OpenMP/tasks
  #==============================

  echo "--------------------------------------------------------"
  echo "---OpenMP/tasks tests---"
  echo "--------------------------------------------------------"

  module swap gcc gcc/4.9.2

  #make MPI_OPTION= OPENMP_OPTION=THREADS NM_VALUE=4

  rm -rf $BUILD_DIR
  mkdir $BUILD_DIR
  pushd $BUILD_DIR

  env NM_VALUE=4 $SCRIPTS_DIR/cmake_openmp_tasks.sh
  make $VERBOSE_ARG

  local num_threads
  for num_threads in 1 2 4 8 16 32 ; do
    perform_runs "-n1 -d16" "" "OMP_NUM_THREADS=$num_threads"
  done

  popd
  rm -rf $BUILD_DIR

  module swap gcc/4.9.2 gcc

  print_results_summary

  #==============================
  # MPI.
  #==============================

  if [ "${PBS_NP:-}" != "" ] ; then

    echo "--------------------------------------------------------"
    echo "---MPI tests---"
    echo "--------------------------------------------------------"

    #make NM_VALUE=4

    rm -rf $BUILD_DIR
    mkdir $BUILD_DIR
    pushd $BUILD_DIR

    env NM_VALUE=4 $SCRIPTS_DIR/cmake_cray_xk7.sh
    make $VERBOSE_ARG

    perform_runs "-n16" ""

    popd
    rm -rf $BUILD_DIR

  fi #---PBS_NP

  print_results_summary

  #==============================
  # MPI + CUDA.
  #==============================

  if [ $IS_TITAN = 1 ] ; then
    module load cudatoolkit
  fi

  if [ "${PBS_NP:-}" != "" ] ; then
  if [ "${PBS_NP:-}" -ge 4 ] ; then

    echo "--------------------------------------------------------"
    echo "---MPI + CUDA tests---"
    echo "--------------------------------------------------------"

    #make CUDA_OPTION=1 NM_VALUE=4

    rm -rf $BUILD_DIR
    mkdir $BUILD_DIR
    pushd $BUILD_DIR

    env NM_VALUE=4 $SCRIPTS_DIR/cmake_cray_xk7_cuda.sh
    make $VERBOSE_ARG

    perform_runs "-n4" "" "CRAY_CUDA_MPS=1"

    popd
    rm -rf $BUILD_DIR

  fi #---PBS_NP
  fi #---PBS_NP

  if [ $IS_TITAN = 1 ] ; then
    module unload cudatoolkit
  fi

  print_results_summary

  #==============================
  # Variants.
  #==============================

  echo "--------------------------------------------------------"
  echo "---Tests of sweeper variants---"
  echo "--------------------------------------------------------"

  local alg_options

  for alg_options in -DSWEEPER_SIMPLE -DSWEEPER_TILEOCTANTS ; do

    #make MPI_OPTION= ALG_OPTIONS="$alg_options" NM_VALUE=16

    rm -rf $BUILD_DIR
    mkdir $BUILD_DIR
    pushd $BUILD_DIR

    env NM_VALUE=16 ALG_OPTIONS=$alg_options $SCRIPTS_DIR/cmake_serial.sh
    make $VERBOSE_ARG

    perform_runs "-n1" ""

    popd
    rm -rf $BUILD_DIR

  done #---alg_options

  print_results_summary

  #==============================
  # Finalize.
  #==============================

  print_results_summary FINAL

}
#==============================================================================

time main

#==============================================================================
