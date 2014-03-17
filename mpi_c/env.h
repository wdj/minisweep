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

#include "env_assert.h"

/*===========================================================================*/
/*---Basic types---*/

/*---Floating point type for sweep---*/

typedef double P;
enum{ P_IS_DOUBLE = 1 };

static inline P P_zero() { return (P)0.; }
static inline P P_one()  { return (P)1.; }

/*---Boolean type---*/

typedef int Bool_t;

enum{ Bool_true = 1, Bool_false = 0 };

/*===========================================================================*/
/*---Struct with environment information---*/

typedef struct
{
  int nproc_x;    /*---Number of procs along x axis---*/
  int nproc_y;    /*---Number of procs along y axis---*/
  int tag;        /*---Next free message tag---*/
} Env;

/*===========================================================================*/
/*---Include MPI wrapper function definitions---*/

#include "env_mpi.h"

/*===========================================================================*/
/*---Include OpenMP wrapper function definitions---*/

#include "env_openmp.h"

/*===========================================================================*/
/*---Initialize for execution---*/

static void Env_initialize( Env *env, int argc, char** argv )
{
  Env_initialize_mpi__( env, argc, argv );
  env->tag = 0;
}

/*===========================================================================*/
/*---Finalize execution---*/

static void Env_finalize( Env* env )
{
  Env_finalize_mpi__( env );
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
  Env_barrier();
  return Env_get_time();
}

/*===========================================================================*/

#endif /*---_mpi_c__env_h_---*/

/*---------------------------------------------------------------------------*/
