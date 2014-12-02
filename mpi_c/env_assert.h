/*---------------------------------------------------------------------------*/
/*!
 * \file   env_assert.h
 * \author Wayne Joubert
 * \date   Wed Jan 15 16:06:28 EST 2014
 * \brief  Environment settings for assertions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _mpi_c__env_assert_h_
#define _mpi_c__env_assert_h_

#include <stdio.h>
#include <assert.h>

#ifdef __cplusplus
extern "C"
{
#endif

/*===========================================================================*/
/*---Assertions---*/

#ifdef USE_CUDA
#define Assert(v) assert(v)
#else
#define Assert(v) assert(v)
#endif

#ifndef Insist
#define Insist( condition ) \
  (void)((condition) || (insist__ (#condition, __FILE__, __LINE__),0))
#endif

static void insist__( const char *condition_string, const char *file, int line )
{
  fprintf( stderr, "Insist error: \"%s\". At file %s, line %i.\n",
                   condition_string, file, line );
  exit( EXIT_FAILURE );
}

#ifndef NDEBUG
#define Static_Assert( condition ) { int a[ ( condition ) ? 1 : -1 ]; (void)a; }
#else
#define Static_Assert( condition )
#endif

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

#endif /*---_mpi_c__env_assert_h_---*/

/*---------------------------------------------------------------------------*/
