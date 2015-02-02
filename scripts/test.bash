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
    pwd
    echo aprun $exec_config_args "$wd/tester" $application_args
    aprun $exec_config_args "$wd/tester" $application_args
    #assert $? = 0
    popd >/dev/null
  elif [ "${IMICROOT:-}" != "" ] ; then
    cp tester $TMPDIR/mic0
    micmpiexec -n 1 -wdir $TMPDIR -host ${HOSTNAME}-mic0 $TMPDIR/tester $application_arg
    rm tester $TMPDIR/mic0
  else
    #assert $exec_config_args = "-n1"
    ./tester $application_args
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
# Initialize build/execution environment.
#==============================================================================
function initialize
{
  if [ "${CRAY_MPICH2_ROOTDIR:-}" != "" ] ; then
    if [ "${PE_ENV:-}" != "GNU" ] ; then
      module swap PrgEnv-pgi PrgEnv-gnu
    fi
    module load cudatoolkit

    if [ "${PE_ENV:-}" != "GNU" ] ; then
      echo "Error: GNU compiler required." 1>&2
      exit 1
    fi
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
  # Serial
  #==============================

  echo "--------------------------------------------------------"
  echo "---Serial tests---"
  echo "--------------------------------------------------------"

  make MPI_OPTION= NM_VALUE=4

  perform_runs "-n1" "" <<EOF
EOF

  #==============================
  # OpenMP
  #==============================

  echo "--------------------------------------------------------"
  echo "---OpenMP tests---"
  echo "--------------------------------------------------------"

  make MPI_OPTION= OPENMP_OPTION=THREADS NM_VALUE=4

  perform_runs "-n1 -d16" "" <<EOF
EOF

  #==============================
  # MPI.
  #==============================

  if [ "${PBS_NP:-}" != "" ] ; then

    echo "--------------------------------------------------------"
    echo "---MPI tests---"
    echo "--------------------------------------------------------"

    make NM_VALUE=4

    perform_runs "-n16" "" <<EOF
EOF

  fi #---PBS_NP

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
EOF

  fi #---PBS_NP
  fi #---PBS_NP

  #==============================
  # Variants.
  #==============================

  echo "--------------------------------------------------------"
  echo "---Tests of sweeper variants---"
  echo "--------------------------------------------------------"

  local alg_options

  for alg_options in -DSWEEPER_SIMPLE -DSWEEPER_TILEOCTANTS ; do

    make MPI_OPTION= ALG_OPTIONS="$alg_options" NM_VALUE=16

     perform_runs "-n1" "" <<EOF
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
