/*---------------------------------------------------------------------------*/
/*!
 * \file   env_openmp.h
 * \author Wayne Joubert
 * \date   Wed Jan 15 16:06:28 EST 2014
 * \brief  Environment settings for openmp.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _mpi_c__env_openmp_h_
#define _mpi_c__env_openmp_h_

#include "types.h"
#include "env_assert.h"
#include "arguments.h"

#ifdef USE_OPENMP
#include "omp.h"
#endif

#ifdef __cplusplus
extern "C"
{
#endif

/*===========================================================================*/
/*---Initialize OpenMP---*/

static void Env_omp_initialize__( Env *env, int argc, char** argv )
{
#ifdef USE_OPENMP
#endif
}

/*===========================================================================*/
/*---Finalize OpenMP---*/

static void Env_omp_finalize__( Env* env )
{
#ifdef USE_OPENMP
#endif
}

/*===========================================================================*/
/*---Set values from args---*/

static void Env_omp_set_values__( Env *env, Arguments* args )
{
#ifdef USE_OPENMP
#endif
}

/*===========================================================================*/
/*---Get openmp current thread number---*/

static inline int Env_omp_thread( const Env* env )
{
  int result = 0;
#ifdef USE_OPENMP
  result = omp_get_thread_num();  
#endif
  return result;
}

/*===========================================================================*/
/*---Are we in an openmp threaded region---*/

static inline Bool_t Env_omp_in_parallel( const Env* env )
{
  Bool_t result = Bool_false;
#ifdef USE_OPENMP
  result = omp_in_parallel();
#endif
  return result;
}

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

#endif /*---_mpi_c__env_openmp_h_---*/

/*---------------------------------------------------------------------------*/
