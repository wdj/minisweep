/*---------------------------------------------------------------------------*/
/*!
 * \file   env.h
 * \author Wayne Joubert
 * \date   Wed Jan 15 16:06:28 EST 2014
 * \brief  Environment settings specific to this programming API.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

/*=============================================================================

This is the main header file providing access to the specific parallel
environemnt being used.  Other header files are included for the candidate APIs
which are activated if being used.

=============================================================================*/

#ifndef _env_h_
#define _env_h_

#include <string.h>
#include "sys/time.h"

#include "arguments.h"

/*---Header file for assertions---*/
#include "env_assert.h"

/*---Data structure declarations relevant to env---*/
#include "env_types.h"

/*---Definitions relevant to specific parallel APIs---*/
#include "env_mpi.h"
#include "env_openmp.h"
#include "env_cuda.h"
#include "env_mic.h"

/*===========================================================================*/

#ifdef __cplusplus
extern "C"
{
#endif

/*===========================================================================*/
/*---Null object---*/

static Env Env_null()
{
  Env result;
  memset( (void*)&result, 0, sizeof(Env) );
  return result;
}

/*===========================================================================*/
/*---Initialize for execution---*/

static void Env_initialize( Env *env, int argc, char** argv )
{
  Env_mpi_initialize_(  env, argc, argv );
  Env_cuda_initialize_( env, argc, argv );
}

/*===========================================================================*/
/*---Set values from args---*/

static void Env_set_values( Env *env, Arguments* args )
{
  Env_mpi_set_values_(  env, args );
  Env_cuda_set_values_( env, args );
}

/*===========================================================================*/
/*---Finalize execution---*/

static void Env_finalize( Env* env )
{
  Env_cuda_finalize_( env );
  Env_mpi_finalize_(  env );
}

/*===========================================================================*/
/*---Indicate whether to do output---*/

static Bool_t Env_is_proc_master( Env* env )
{
  return ( Env_is_proc_active( env ) && Env_proc_this( env ) == 0 );
}

/*===========================================================================*/
/*---Timer type---*/

typedef double Timer;

/*===========================================================================*/
/*---Timer utilities---*/

static Timer Env_get_time( Env* env )
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
  return Env_get_time( env );
}

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

#endif /*---_env_h_---*/

/*---------------------------------------------------------------------------*/
