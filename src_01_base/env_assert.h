/*---------------------------------------------------------------------------*/
/*!
 * \file   env_assert.h
 * \author Wayne Joubert
 * \date   Wed Jan 15 16:06:28 EST 2014
 * \brief  Environment settings for assertions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _env_assert_h_
#define _env_assert_h_

#ifndef __CUDA_ARCH__

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C"
{
#endif

/*===========================================================================*/
/*---Assertions---*/

#define Assert(v) assert(v)

#ifndef Insist
#define Insist( condition ) \
  (void)((condition) || (insist_ (#condition, __FILE__, __LINE__),0))
#endif

static void insist_( const char *condition_string, const char *file, int line )
{
  fprintf( stderr, "Insist error: \"%s\". At file %s, line %i.\n",
                   condition_string, file, line );
  exit( EXIT_FAILURE );
}

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

#else /*---__CUDA_ARCH__---*/

/*---Ignore on device.---*/

#define Assert(v)
#define Insist(v)

#endif /*---__CUDA_ARCH__---*/

/*===========================================================================*/
/*---Static assertions---*/

#ifndef NDEBUG
/*---Fail compilation and (hopefully) give a filename/line number---*/
#define Static_Assert( condition ) { int a[ ( condition ) ? 1 : -1 ]; (void)a; }
#else
#define Static_Assert( condition )
#endif

#endif /*---_env_assert_h_---*/

/*---------------------------------------------------------------------------*/
