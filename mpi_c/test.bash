#!/bin/bash
#==============================================================================

set -eu

declare g_ntest=0
declare g_ntest_passed=0
declare g_verbose=1

#==============================================================================
function perform_run
{
  local exec_args="$1"
  local app_args="$2"

  if [ "${PBS_NP:-}" != "" -a "${CRAY_MPICH2_ROOTDIR:-}" != "" ] ; then
    local wd="$PWD"
    pushd "$MEMBERWORK" >/dev/null
    aprun $exec_args "$wd/sweep.x" $app_args
    #assert $? = 0
    popd >/dev/null
  else
    #assert $exec_args = "-n1"
    ./sweep.x $app_args
    #assert $? = 0
  fi

}
#==============================================================================

#==============================================================================
function compare_runs
{
  local app_args1="$2"
  local app_args2="$4"

  local exec_args1="$1"
  local exec_args2="$3"

  local is_pass is_pass1 is_pass2
  local normsq1 normsq2
  local time1 time2

  g_ntest=$(( $g_ntest + 1 ))

  #---Run 1.

  echo -n "$g_ntest // $exec_args1 / $app_args1 / "
  local result1="$( perform_run "$exec_args1" "$app_args1" )"
  normsq1=$( echo "$result1" | grep '^Normsq result: ' \
                             | sed -e 's/^Normsq result: *//' -e 's/ .*//' )
  time1=$(   echo "$result1" | grep '^Normsq result: ' \
                             | sed -e 's/.* time: *//' -e 's/ .*//' )
  echo -n "$time1 // "
  is_pass1=$([ "$result1" != "${result1/ PASS }" ] && echo 1 || echo 0 )
  [[ g_verbose -ge 2 ]] && echo "$result1"

  #---Run 2.

  echo -n "$exec_args2 / $app_args2 / "
  local result2="$( perform_run "$exec_args2" "$app_args2" )"
  normsq2=$( echo "$result2" | grep '^Normsq result: ' \
                             | sed -e 's/^Normsq result: *//' -e 's/ .*//' )
  time2=$(   echo "$result2" | grep '^Normsq result: ' \
                             | sed -e 's/.* time: *//' -e 's/ .*//' )
  echo -n "$time2 // "
  is_pass2=$([ "$result2" != "${result2/ PASS }" ] && echo 1 || echo 0 )
  [[ g_verbose -ge 2 ]] && echo "$result2"

  #---

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
function main
{
# nx ny nz ne nm na numiterations nproc_x nproc_y nblock_z 

  if [ "${PBS_NP:-}" != "" ] ; then

    echo "----------------"
    echo "---MPI + CUDA tests---"
    echo "----------------"

    if [ "$PBS_NUM_NODES" -ge 4 ] ; then

      make CUDA_OPTION=1

      local ARGS="--nx 3 --ny 5 --nz 2 --ne 2 --nm 4 --na 5 --nblock_z 2"

      compare_runs \
        "-n1 -N1" "$ARGS " \
        "-n2 -N1" "$ARGS --is_using_device 1 --nproc_x 2 --nproc_y 1"
      compare_runs \
        "-n1 -N1" "$ARGS " \
        "-n2 -N1" "$ARGS --is_using_device 1 --nproc_x 1 --nproc_y 2"
      compare_runs \
        "-n1 -N1" "$ARGS " \
        "-n4 -N1" "$ARGS --is_using_device 1 --nproc_x 2 --nproc_y 2"
      compare_runs \
        "-n1 -N1" "$ARGS " \
        "-n4 -N1" "$ARGS --is_using_device 1 --nproc_x 2 --nproc_y 2 --nthread_e 3"
      compare_runs \
        "-n1 -N1" "$ARGS " \
        "-n4 -N1" "$ARGS --is_using_device 1 --nproc_x 2 --nproc_y 2 --nthread_e 3 --nthread_octant 8"

    fi #---PBS_NUM_NODES

    echo "----------------"
    echo "---CUDA tests---"
    echo "----------------"

    make CUDA_OPTION=1

    local ARGS="--nx  2 --ny  3 --nz  4 --ne 20 --nm 4 --na 5 --nblock_z 2"

    compare_runs \
      "-n1"  "$ARGS " \
      "-n1"  "$ARGS --is_using_device 1 --nthread_e 1 --nthread_octant 1"
    compare_runs \
      "-n1"  "$ARGS " \
      "-n1"  "$ARGS --is_using_device 1 --nthread_e 1 --nthread_octant 2"
    compare_runs \
      "-n1"  "$ARGS " \
      "-n1"  "$ARGS --is_using_device 1 --nthread_e 1 --nthread_octant 4"
    compare_runs \
      "-n1"  "$ARGS " \
      "-n1"  "$ARGS --is_using_device 1 --nthread_e 1 --nthread_octant 8"
    compare_runs \
      "-n1"  "$ARGS " \
      "-n1"  "$ARGS --is_using_device 1 --nthread_e 2 --nthread_octant 8"
    compare_runs \
      "-n1"  "$ARGS " \
      "-n1"  "$ARGS --is_using_device 1 --nthread_e 10 --nthread_octant 8"
    compare_runs \
      "-n1"  "$ARGS " \
      "-n1"  "$ARGS --is_using_device 1 --nthread_e 20 --nthread_octant 8"

    echo "---------------"
    echo "---MPI tests---"
    echo "---------------"

    make

    local ARGS="--nx  5 --ny  4 --nz  5 --ne 7 --nm 4 --na 10"
    compare_runs   "-n1"  "$ARGS --nproc_x 1 --nproc_y 1 --nblock_z 1" \
                   "-n2"  "$ARGS --nproc_x 2 --nproc_y 1 --nblock_z 1"
    compare_runs   "-n1"  "$ARGS --nproc_x 1 --nproc_y 1 --nblock_z 1" \
                   "-n2"  "$ARGS --nproc_x 1 --nproc_y 2 --nblock_z 1"

    local ARGS="--nx  5 --ny  4 --nz  6 --ne 7 --nm 4 --na 10"
    compare_runs   "-n1"  "$ARGS --nproc_x 1 --nproc_y 1 --nblock_z 1" \
                  "-n16"  "$ARGS --nproc_x 4 --nproc_y 4 --nblock_z 2"

    local ARGS="--nx 5 --ny 8 --nz 16 --ne 9 --nm 1 --na 12"
    compare_runs  "-n16"  "$ARGS --nproc_x 4 --nproc_y 4 --nblock_z 1" \
                  "-n16"  "$ARGS --nproc_x 4 --nproc_y 4 --nblock_z 2"
    compare_runs  "-n16"  "$ARGS --nproc_x 4 --nproc_y 4 --nblock_z 2" \
                  "-n16"  "$ARGS --nproc_x 4 --nproc_y 4 --nblock_z 4"

    echo "------------------"
    echo "---OpenMP tests---"
    echo "------------------"

    make OPENMP_OPTION=THREADS

    local ARGS="--nx  5 --ny  4 --nz  5 --ne 200 --nm 4 --na 10"
    compare_runs   "-n1 -d1"  "$ARGS --nthread_e 1" \
                   "-n1 -d2"  "$ARGS --nthread_e 2"
    compare_runs   "-n1 -d2"  "$ARGS --nthread_e 2" \
                   "-n1 -d3"  "$ARGS --nthread_e 3"
    compare_runs   "-n1 -d3"  "$ARGS --nthread_e 3" \
                   "-n1 -d4"  "$ARGS --nthread_e 4"

    compare_runs   "-n1 -d1"  "$ARGS --nthread_octant 1" \
                   "-n1 -d2"  "$ARGS --nthread_octant 2"
    compare_runs   "-n1 -d2"  "$ARGS --nthread_octant 2" \
                   "-n1 -d4"  "$ARGS --nthread_octant 4"
    compare_runs   "-n1 -d4"  "$ARGS --nthread_octant 4" \
                   "-n1 -d8"  "$ARGS --nthread_octant 8"

    compare_runs   "-n1 -d1"  "$ARGS --nthread_e 1 --nthread_octant 1" \
                   "-n1 -d2"  "$ARGS --nthread_e 2 --nthread_octant 1"
    compare_runs   "-n1 -d2"  "$ARGS --nthread_e 2 --nthread_octant 1" \
                   "-n1 -d4"  "$ARGS --nthread_e 2 --nthread_octant 2"

  fi

  echo "---Tests of sweeper variants---"

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

    local ARGS="--nx  5 --ny  5 --nz  5 --ne 10 --nm 16 --na 20"
    compare_runs  "-n1" "$ARGS --niterations 1 $ARG_NBLOCK_Z_1" \
                  "-n1" "$ARGS --niterations 2 $ARG_NBLOCK_Z_1"
    compare_runs  "-n1" "$ARGS --niterations 1 $ARG_NBLOCK_Z_1" \
                  "-n1" "$ARGS --niterations 1 $ARG_NBLOCK_Z_5"

  done

  echo -n "Total tests $g_ntest PASSED $g_ntest_passed"
  echo " FAILED $(( $g_ntest - $g_ntest_passed ))."

}
#==============================================================================

main

#==============================================================================
