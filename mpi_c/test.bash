#!/bin/bash
#==============================================================================

set -eu

declare g_num_tests=0
declare g_num_passed=0
declare g_verbose=1

#==============================================================================
function perform_run
{
  local procs="$1"
  local args="$2"

  if [ "${PBS_NP:-}" != "" -a "${CRAY_MPICH2_ROOTDIR:-}" != "" ] ; then
    local wd="$PWD"
    pushd "$MEMBERWORK" >/dev/null
    aprun -n$procs "$wd/sweep.x" $args
    #assert $? = 0
    popd >/dev/null
  else
    #assert $procs = 1
    ./sweep.x $args
    #assert $? = 0
  fi

}
#==============================================================================

#==============================================================================
function compare_runs
{
  local procs1="$1"
  local args1="$2"
  local procs2="$3"
  local args2="$4"

  local pass pass1 pass2

  local result1="$( perform_run "$procs1" "$args1" )"
  [[ g_verbose -ge 2 ]] && echo "$result1"
  local result2="$( perform_run "$procs2" "$args2" )"
  [[ g_verbose -ge 2 ]] && echo "$result2"

  g_num_tests=$(( $g_num_tests + 1 ))

  pass1=$([ "$result1" != "${result1/ PASS }" ] && echo 1 || echo 0 )
  pass2=$([ "$result2" != "${result2/ PASS }" ] && echo 1 || echo 0 )

  local normsq1 normsq2
  normsq1=$( echo "$result1" | grep '^Normsq result: ' \
                             | sed -e 's/^Normsq result: *//' -e 's/ .*//' )
  normsq2=$( echo "$result2" | grep '^Normsq result: ' \
                             | sed -e 's/^Normsq result: *//' -e 's/ .*//' )

  pass=$([ "$normsq1" = "$normsq2" ] && echo 1 || echo 0 )

  if [ $pass1 = 1 -a $pass2 = 1 -a $pass = 1 ] ; then
    echo "$g_num_tests / $procs1 / $args1 / $procs2 / $args2 / PASS"
    g_num_passed=$(( $g_num_passed + 1 ))
  else
    echo "$g_num_tests / $procs1 / $args1 / $procs2 / $args2 / FAIL"
    echo "$result1"
    echo "$result2"
  fi


}
#==============================================================================

#==============================================================================
function main
{
# nx ny nz ne nm na numiterations nproc_x nproc_y nblock_z 

  if [ "${PBS_NP:-}" != "" ] ; then

    make

    local ARGS_SIZES="--nx  5 --ny  5 --nz  5 --ne 10 --nm 16 --na 20"
    compare_runs   "1"  "$ARGS_SIZES --nproc_x 1 --nproc_y 1 --nblock_z 1" \
                   "2"  "$ARGS_SIZES --nproc_x 2 --nproc_y 1 --nblock_z 1"
    compare_runs   "1"  "$ARGS_SIZES --nproc_x 1 --nproc_y 1 --nblock_z 1" \
                   "2"  "$ARGS_SIZES --nproc_x 1 --nproc_y 2 --nblock_z 1"

    local ARGS_SIZES="--nx  5 --ny  5 --nz  6 --ne 10 --nm 16 --na 20"
    compare_runs   "1"  "$ARGS_SIZES --nproc_x 1 --nproc_y 1 --nblock_z 1" \
                  "16"  "$ARGS_SIZES --nproc_x 4 --nproc_y 4 --nblock_z 2"

    local ARGS_SIZES="--nx 16 --ny 32 --nz 64 --ne 16 --nm 16 --na 32"
    compare_runs  "16"  "$ARGS_SIZES --nproc_x 4 --nproc_y 4 --nblock_z 1" \
                  "16"  "$ARGS_SIZES --nproc_x 4 --nproc_y 4 --nblock_z 2"
    compare_runs  "16"  "$ARGS_SIZES --nproc_x 4 --nproc_y 4 --nblock_z 2" \
                  "16"  "$ARGS_SIZES --nproc_x 4 --nproc_y 4 --nblock_z 4"

  fi

  local alg_options

  for alg_options in -DSWEEPER_KBA -DSWEEPER_SIMPLE -DSWEEPER_TILEOCTANTS ; do

    make MPI_OPTION= ALG_OPTIONS="$alg_options"

    if [ $alg_options = "-DSWEEPER_KBA" ] ; then
      local ARG_NBLOCK_Z_1="--nblock_z 1"
      local ARG_NBLOCK_Z_5="--nblock_z 5"
    else
      local ARG_NBLOCK_Z_1=""
      local ARG_NBLOCK_Z_5=""
    fi

    local ARGS_SIZES="--nx  5 --ny  5 --nz  5 --ne 10 --nm 16 --na 20"
    compare_runs  "1" "$ARGS_SIZES --niterations 1 $ARG_NBLOCK_Z_1" \
                  "1" "$ARGS_SIZES --niterations 2 $ARG_NBLOCK_Z_1"
    compare_runs  "1" "$ARGS_SIZES --niterations 1 $ARG_NBLOCK_Z_1" \
                  "1" "$ARGS_SIZES --niterations 1 $ARG_NBLOCK_Z_5"

  done

  echo -n "Total tests $g_num_tests PASSED $g_num_passed"
  echo " FAILED $(( $g_num_tests - $g_num_passed ))."

}
#==============================================================================

main

#==============================================================================
