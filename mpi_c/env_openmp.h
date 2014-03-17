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
