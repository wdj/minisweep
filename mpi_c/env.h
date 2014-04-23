/*---------------------------------------------------------------------------*/
/*!
 * \file   env.h
 * \author Wayne Joubert
 * \date   Wed Jan 15 16:06:28 EST 2014
 * \brief  Environment settings specific to this programming API.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _mpi_c__env_h_
#define _mpi_c__env_h_

#include <stdio.h>
#include <stdlib.h>

#include <sys/time.h>

#include "arguments.h"

/*===========================================================================*/
/*---Header file for assertions---*/

#include "env_assert.h"

/*===========================================================================*/
/*---Declarations relevant to environment---*/

#include "env_declarations.h"

/*===========================================================================*/
/*---MPI wrapper function definitions---*/

#include "env_mpi.h"

/*===========================================================================*/
/*---OpenMP wrapper function definitions---*/

#include "env_openmp.h"

/*===========================================================================*/
/*---CUDA wrapper function definitions---*/

#include "env_cuda.h"

/*===========================================================================*/

#ifdef __cplusplus
extern "C"
{
#endif

/*===========================================================================*/
/*---Initialize for execution---*/

static void Env_initialize( Env *env, int argc, char** argv )
{
  Env_mpi_initialize__(  env, argc, argv );
  Env_omp_initialize__(  env, argc, argv );
  Env_cuda_initialize__( env, argc, argv );
}

/*===========================================================================*/
/*---Set values from args---*/

static void Env_set_values( Env *env, Arguments* args )
{
  Env_mpi_set_values__(  env, args );
  Env_omp_set_values__(  env, args );
  Env_cuda_set_values__( env, args );
}

/*===========================================================================*/
/*---Finalize execution---*/

static void Env_finalize( Env* env )
{
  Env_cuda_finalize__( env );
  Env_omp_finalize__(  env );
  Env_mpi_finalize__(  env );
}

/*===========================================================================*/
/*---Indicate whether to do output---*/

static Bool_t Env_do_output( const Env* env )
{
  return ( Env_proc_this( env ) == 0 );
}

/*===========================================================================*/
/*---Timer type---*/

typedef double Timer_t;

/*===========================================================================*/
/*---Timer utilities---*/

static Timer_t Env_get_time()
{
    struct timeval tv;
    int i = gettimeofday( &tv, NULL );
    Timer_t result = ( (Timer_t) tv.tv_sec +
                       (Timer_t) tv.tv_usec * 1.e-6 );
    return result;
}

/*---------------------------------------------------------------------------*/

static Timer_t Env_get_synced_time()
{
  Env_mpi_barrier();
  return Env_get_time();
}

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

#endif /*---_mpi_c__env_h_---*/

/*---------------------------------------------------------------------------*/
