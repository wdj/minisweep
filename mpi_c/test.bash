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

  compare_runs   "1"  "5 5 5 10 16 20  1  1 1  1" \
                 "1"  "5 5 5 10 16 20  2  1 1  1"
  compare_runs   "1"  "5 5 5 10 16 20  1  1 1  1" \
                 "1"  "5 5 5 10 16 20  1  1 1  5"

  if [ "${PBS_NP:-}" != "" ] ; then

    compare_runs   "1"  "5 5 5 10 16 20  1  1 1  1" \
                   "2"  "5 5 5 10 16 20  1  2 1  1"
    compare_runs   "1"  "5 5 5 10 16 20  1  1 1  1" \
                   "2"  "5 5 5 10 16 20  1  1 2  1"

    compare_runs   "1"  "5 5 6 10 16 20  1  1 1  1" \
                  "16"  "5 5 6 10 16 20  1  4 4  2"

    compare_runs  "16"  "16 32 64 16 16 32  1  4 4  1" \
                  "16"  "16 32 64 16 16 32  1  4 4  2"
    compare_runs  "16"  "16 32 64 16 16 32  1  4 4  2" \
                  "16"  "16 32 64 16 16 32  1  4 4  4"

  fi

  echo -n "Total tests $g_num_tests PASSED $g_num_passed"
  echo " FAILED $(( $g_num_tests - $g_num_passed ))."

}
#==============================================================================

main

#==============================================================================
