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

#include "env_assert.h"

#ifdef USE_OPENMP
#include "omp.h"
#endif

/*===========================================================================*/
/*---Set up enums---*/

#ifdef USE_OPENMP_OCTANT
enum{ IS_USING_OPENMP_OCTANT = 1 };
#else
enum{ IS_USING_OPENMP_OCTANT = 0 };
#endif

#ifdef USE_OPENMP_E
enum{ IS_USING_OPENMP_E = 1 };
#else
enum{ IS_USING_OPENMP_E = 0 };
#endif

/*===========================================================================*/
/*---Initialize OpenMP---*/

static void Env_initialize_openmp__( Env *env, int argc, char** argv )
{
#ifdef USE_OPENMP
  omp_set_nested( 1 );
#endif
}

/*===========================================================================*/
/*---Finalize OpenMP---*/

static void Env_finalize_openmp__( Env* env )
{
#ifdef USE_OPENMP
  omp_set_nested( 0 );
#endif
}


/*===========================================================================*/
/*---Get openmp current number of threads---*/

static inline int Env_num_threads( const Env* env )
{
  int result = 1;
#ifdef USE_OPENMP
  result = omp_get_num_threads();
#endif
  return result;
}

/*===========================================================================*/
/*---Get openmp current thread number---*/

static inline int Env_thread_this( const Env* env )
{
  int result = 0;
#ifdef USE_OPENMP
  result = omp_get_thread_num();  
#endif
  return result;
}

/*===========================================================================*/

#endif /*---_mpi_c__env_openmp_h_---*/

/*---------------------------------------------------------------------------*/
