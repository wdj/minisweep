/*---------------------------------------------------------------------------*/
/*!
 * \file   env.h
 * \author Wayne Joubert
 * \date   Wed Jan 15 16:06:28 EST 2014
 * \brief  Environment settings specific to this programming API.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _env_h_
#define _env_h_

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
/*---Intel MIC definitions---*/

#include "env_mic.h"

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
/*---Reset values---*/

static void Env_reset_values( Env *env )
{
  Env_mpi_reset_values__( env );
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

static Bool_t Env_is_proc_master( const Env* env )
{
  return ( Env_is_proc_active( env ) && Env_proc_this( env ) == 0 );
}

/*===========================================================================*/
/*---Timer type---*/

typedef double Timer;

/*===========================================================================*/
/*---Timer utilities---*/

static Timer Env_get_time()
{
    struct timeval tv;
    int i = gettimeofday( &tv, NULL );
    Timer result = ( (Timer) tv.tv_sec +
                     (Timer) tv.tv_usec * 1.e-6 );
    return result;
}

/*---------------------------------------------------------------------------*/

static Timer Env_get_synced_time( Env* env )
{
  Env_mpi_barrier( env );
  return Env_get_time();
}

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

#endif /*---_env_h_---*/

/*---------------------------------------------------------------------------*/
