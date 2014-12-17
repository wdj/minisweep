#!/bin/bash -l
#==============================================================================
# Run tests on the code.
#==============================================================================

set -eu

declare g_ntest=0
declare g_ntest_passed=0
declare g_verbose=1

#==============================================================================
# Run the executable once.
#==============================================================================
function perform_run
{
  local exec_config_args="$(echo "$1" | tr '\n' ' ')"
  local application_args="$(echo "$2" | tr '\n' ' ')"

  if [ "${PBS_NP:-}" != "" -a "${CRAY_MPICH2_ROOTDIR:-}" != "" ] ; then
    #---If running on Cray back-end node, must cd to Lustre to do the aprun.
    local wd="$PWD"
    pushd "$MEMBERWORK" >/dev/null
    aprun $exec_config_args "$wd/sweep.x" $application_args
    #assert $? = 0
    popd >/dev/null
  else
    #assert $exec_config_args = "-n1"
    ./sweep.x $application_args
    #assert $? = 0
  fi
}
#==============================================================================

#==============================================================================
# Run the executable multiple times using stdin.
#==============================================================================
function perform_runs
{
  local exec_config_args="$(echo "$1" | tr '\n' ' ')"
  local application_args="$(echo "$2" | tr '\n' ' ')"

  perform_run "$1" "$2" | tee tmp_

  local ntest=$(        grep TESTS tmp_ | sed -e 's/.*TESTS *//'  -e 's/ .*//' )
  local ntest_passed=$( grep TESTS tmp_ | sed -e 's/.*PASSED *//' -e 's/ .*//' )
  rm -f tmp_

  g_ntest=$((        $g_ntest        + $ntest ))
  g_ntest_passed=$(( $g_ntest_passed + $ntest_passed ))
}
#==============================================================================

#==============================================================================
# Compare the output of two runs for match.
# This function is DEFUNCT.
#==============================================================================
function compare_runs
{
  local exec_config_args1="$(echo "$1" | tr '\n' ' ')"
  local exec_config_args2="$(echo "$3" | tr '\n' ' ')"

  local application_args1="$(echo "$2" | tr '\n' ' ')"
  local application_args2="$(echo "$4" | tr '\n' ' ')"

  local is_pass
  local is_pass1
  local is_pass2
  local normsq1
  local normsq2
  local time1
  local time2

  g_ntest=$(( $g_ntest + 1 ))

  #---Run 1.

  echo -n "$g_ntest // $exec_config_args1 / $application_args1 / "
  local result1="$( perform_run "$exec_config_args1" "$application_args1" )"
  normsq1=$( echo "$result1" | grep '^Normsq result: ' \
                             | sed -e 's/^Normsq result: *//' -e 's/ .*//' )
  time1=$(   echo "$result1" | grep '^Normsq result: ' \
                             | sed -e 's/.* time: *//' -e 's/ .*//' )
  echo -n "$time1 // "
  is_pass1=$([ "$result1" != "${result1/ PASS }" ] && echo 1 || echo 0 )
  [[ g_verbose -ge 2 ]] && echo "$result1"

  #---Run 2.

  echo -n "$exec_config_args2 / $application_args2 / "
  local result2="$( perform_run "$exec_config_args2" "$application_args2" )"
  normsq2=$( echo "$result2" | grep '^Normsq result: ' \
                             | sed -e 's/^Normsq result: *//' -e 's/ .*//' )
  time2=$(   echo "$result2" | grep '^Normsq result: ' \
                             | sed -e 's/.* time: *//' -e 's/ .*//' )
  echo -n "$time2 // "
  is_pass2=$([ "$result2" != "${result2/ PASS }" ] && echo 1 || echo 0 )
  [[ g_verbose -ge 2 ]] && echo "$result2"

  #---Final test of success.

  # Check whether each run passed and whether the results match each other.

  is_pass=$([ "$normsq1" = "$normsq2" ] && echo 1 || echo 0 )

  if [ $is_pass1 = 1 -a $is_pass2 = 1 -a $is_pass = 1 ] ; then
    echo "PASS"
    g_ntest_passed=$(( $g_ntest_passed + 1 ))
  else
    echo "FAIL"
    echo "$result1"
    echo "$result2"
  fi
}
#==============================================================================

#==============================================================================
# Strings for cases to run, MPI.
#==============================================================================
function argstrings_mpi
{
  local ARGS="--ncell_x  5 --ncell_y  4 --ncell_z  5 --ne 7 --na 10"

  echo "$ARGS --nproc_x 1 --nproc_y 1 --nblock_z 1"
  echo "$ARGS --nproc_x 2 --nproc_y 1 --nblock_z 1"

  echo "$ARGS --nproc_x 1 --nproc_y 1 --nblock_z 1"
  echo "$ARGS --nproc_x 1 --nproc_y 2 --nblock_z 1"

  local ARGS="--ncell_x  5 --ncell_y  4 --ncell_z  6 --ne 7 --na 10"

  echo "$ARGS --nproc_x 1 --nproc_y 1 --nblock_z 1"
  echo "$ARGS --nproc_x 4 --nproc_y 4 --nblock_z 2"

  local ARGS="--ncell_x  5 --ncell_y  4 --ncell_z  6 --ne 7 --na 10 --is_face_comm_async 0"

  echo "$ARGS --nproc_x 1 --nproc_y 1 --nblock_z 1"
  echo "$ARGS --nproc_x 4 --nproc_y 4 --nblock_z 2"

  local ARGS="--ncell_x 5 --ncell_y 8 --ncell_z 16 --ne 9 --na 12"

  echo "$ARGS --nproc_x 4 --nproc_y 4 --nblock_z 1"
  echo "$ARGS --nproc_x 4 --nproc_y 4 --nblock_z 2"

  echo "$ARGS --nproc_x 4 --nproc_y 4 --nblock_z 2"
  echo "$ARGS --nproc_x 4 --nproc_y 4 --nblock_z 4"
}
#==============================================================================

#==============================================================================
# Strings for cases to run, CUDA.
#==============================================================================
function argstrings_cuda
{
  local ARGS="--ncell_x  2 --ncell_y  3 --ncell_z  4 --ne 20 --na 5 --nblock_z 2"

  for nthread_octant in 1 2 4 8 ; do
    echo "$ARGS"
    echo "$ARGS --is_using_device 1 --nthread_e 1 \
                --nthread_octant $nthread_octant"
  done

  for nthread_e in 2 10 20 ; do
    echo "$ARGS"
    echo "$ARGS --is_using_device 1 --nthread_e $nthread_e --nthread_octant 8"
  done
}
#==============================================================================

#==============================================================================
# Strings for cases to run, MPI + CUDA.
#==============================================================================
function argstrings_mpi_cuda
{
  local ARGS="--ncell_x 3 --ncell_y 5 --ncell_z 6 --ne 2 --na 5 --nblock_z 2"

  local nproc_x=""
  local nproc_y=""

  for nproc_x in 1 2 ; do
  for nproc_y in 1 2 ; do
  for nthread_octant in 1 2 4 8 ; do
      echo "$ARGS"
      echo "$ARGS --is_using_device 1 --nproc_x $nproc_x --nproc_y $nproc_y \
                  --nthread_e 3 --nthread_octant $nthread_octant"
  done
  done
  done
}
#==============================================================================

#==============================================================================
# Strings for cases to run, serial case.
#==============================================================================
function argstrings_serial
{
  local key=""

  for key in {0..1} ; do

    local ARGS="--ncell_x $(( 2 + 3 * $key )) --ncell_y  $(( 1 + 3 * $key )) --ncell_z  $(( 2 + 3 * $key )) --ne 3 --na 7"

    local ncell_x_per_subblock=""
    local ncell_y_per_subblock=""
    local ncell_z_per_subblock=""

    for ncell_x_per_subblock in  {1..3} ; do
    for ncell_y_per_subblock in  {1..3} ; do
    for ncell_z_per_subblock in  {1..3} ; do
      echo "$ARGS"
      echo "$ARGS --ncell_x_per_subblock $ncell_x_per_subblock --ncell_y_per_subblock $ncell_y_per_subblock --ncell_z_per_subblock $ncell_z_per_subblock"
    done
    done
    done
  done
}
#==============================================================================

#==============================================================================
# Strings for cases to run, OpenMP.
#==============================================================================
function argstrings_openmp
{
#  local ARGS="--ncell_x 1 --ncell_y 2 --ncell_z  1 --ne 1 --na 1"
#
#  echo "$ARGS --ncell_x_per_subblock 1 --ncell_y_per_subblock 1 --ncell_z_per_subblock 1 --nthread_y 1"
#  echo "$ARGS --ncell_x_per_subblock 1 --ncell_y_per_subblock 1 --ncell_z_per_subblock 1 --nthread_y 2"

  local ARGS="--ncell_x  5 --ncell_y  4 --ncell_z  5 --ne 17 --na 10"

  echo "$ARGS --nthread_e 1"
  echo "$ARGS --nthread_e 2"
  echo "$ARGS --nthread_e 2"
  echo "$ARGS --nthread_e 3"
  echo "$ARGS --nthread_e 3"
  echo "$ARGS --nthread_e 4"

  echo "$ARGS --nthread_octant 1"
  echo "$ARGS --nthread_octant 2"
  echo "$ARGS --nthread_octant 2"
  echo "$ARGS --nthread_octant 4"
  echo "$ARGS --nthread_octant 4"
  echo "$ARGS --nthread_octant 8"

  echo "$ARGS --nthread_e 1 --nthread_octant 1"
  echo "$ARGS --nthread_e 2 --nthread_octant 1"
  echo "$ARGS --nthread_e 2 --nthread_octant 1"
  echo "$ARGS --nthread_e 2 --nthread_octant 2"
}
#==============================================================================

#==============================================================================
# Strings for cases to run, sweeper variants.
#==============================================================================
function argstrings_variants
{
  local ARG_NBLOCK_Z_1="$1"
  local ARG_NBLOCK_Z_5="$2"

  local ARGS="--ncell_x  4 --ncell_y  3 --ncell_z  5 --ne 11 --na 7"

  echo "$ARGS --niterations 1 $ARG_NBLOCK_Z_1"
  echo "$ARGS --niterations 2 $ARG_NBLOCK_Z_1"

  echo "$ARGS --niterations 1 $ARG_NBLOCK_Z_1"
  echo "$ARGS --niterations 1 $ARG_NBLOCK_Z_5"
}
#==============================================================================

#==============================================================================
# Initialize build/execution environment.
#==============================================================================
function initialize
{
  if [ "$PE_ENV" != "GNU" ] ; then
    module swap PrgEnv-pgi PrgEnv-gnu
  fi
  module load cudatoolkit

  if [ "$PE_ENV" != "GNU" ] ; then
    echo "Error: GNU compiler required." 1>&2
    exit 1
  fi
}
#==============================================================================

#==============================================================================
function main
{
  initialize

  #---args to use below:
  #---  ncell_x ncell_y ncell_z ne nm na numiterations nproc_x nproc_y nblock_z 

  cp /dev/null tmp_

  #==============================
  # MPI + CUDA.
  #==============================

  if [ "${PBS_NP:-}" != "" ] ; then
  if [ "${PBS_NP:-}" -ge 4 ] ; then

    echo "--------------------------------------------------------"
    echo "---MPI + CUDA tests---"
    echo "--------------------------------------------------------"

    make CUDA_OPTION=1 NM_VALUE=4

    perform_runs "-n4" "" <<EOF
$(argstrings_mpi_cuda)
EOF

  fi #---PBS_NP
  fi #---PBS_NP

  #==============================
  # CUDA.
  #==============================

  if [ "${PBS_NP:-}" != "" ] ; then

    echo "--------------------------------------------------------"
    echo "---CUDA tests---"
    echo "--------------------------------------------------------"

    make CUDA_OPTION=1 NM_VALUE=4

    perform_runs "-n1" "" <<EOF
$(argstrings_cuda)
EOF

  fi #---PBS_NP

  #==============================
  # MPI.
  #==============================

  if [ "${PBS_NP:-}" != "" ] ; then

    echo "--------------------------------------------------------"
    echo "---MPI tests---"
    echo "--------------------------------------------------------"

    make NM_VALUE=4

    perform_runs "-n16" "" <<EOF
$(argstrings_mpi)
EOF

  fi #---PBS_NP

  #==============================
  # OpenMP
  #==============================

  echo "--------------------------------------------------------"
  echo "---OpenMP tests---"
  echo "--------------------------------------------------------"

  make MPI_OPTION= OPENMP_OPTION=THREADS NM_VALUE=4

  perform_runs "-n1 -d8" "" <<EOF
$(argstrings_openmp)
EOF

  #==============================
  # Serial
  #==============================

  echo "--------------------------------------------------------"
  echo "---Serial tests---"
  echo "--------------------------------------------------------"

  make MPI_OPTION= NM_VALUE=4

  perform_runs "-n1" "" <<EOF
$(argstrings_serial)
EOF

  #==============================
  # Variants.
  #==============================

  echo "--------------------------------------------------------"
  echo "---Tests of sweeper variants---"
  echo "--------------------------------------------------------"

  local alg_options

  for alg_options in -DSWEEPER_KBA -DSWEEPER_SIMPLE -DSWEEPER_TILEOCTANTS ; do

    make MPI_OPTION= ALG_OPTIONS="$alg_options" NM_VALUE=16

    if [ $alg_options = "-DSWEEPER_KBA" ] ; then
      local ARG_NBLOCK_Z_1="--nblock_z 1"
      local ARG_NBLOCK_Z_5="--nblock_z 5"
    else
      local ARG_NBLOCK_Z_1=""
      local ARG_NBLOCK_Z_5=""
    fi

     perform_runs "-n1" "" <<EOF
$(argstrings_variants "$ARG_NBLOCK_Z_1" "$ARG_NBLOCK_Z_5")
EOF

  done #---alg_options

  #==============================
  # Finalize.
  #==============================

  local ntest_failed=$(( $g_ntest - $g_ntest_passed ))

  echo
  echo -n "TOTAL: "
  echo -n "TESTS $g_ntest"
  echo -n "    "
  echo -n "PASSED $g_ntest_passed"
  echo -n "    "
  echo -n "FAILED $(( $g_ntest - $g_ntest_passed ))"
  echo -n "."
  if [ $ntest_failed = 0 ] ; then
    echo -n "          ========== SUCCESS =========="
  fi
  echo


}
#==============================================================================

time main

#==============================================================================
