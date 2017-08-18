/*---------------------------------------------------------------------------*/
/*!
 * \file   env_mic.h
 * \author Wayne Joubert
 * \date   Wed Jun 11 09:33:15 EDT 2014
 * \brief  Environment settings for Intel MIC.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _env_mic_h_
#define _env_mic_h_

#include <stddef.h>

#include "types.h"
#include "env_mic_kernels.h"

/* TODO: move functions to .c file */

/*===========================================================================*/
/*---Memory management---*/

#ifdef __MIC__

static int* malloc_host_int( size_t n, Env* env )
{
  Insist( n+1 >= 1 );
  int* result = (int*)malloc( n * sizeof(int) );
  Insist( result );
  return result;
}

/*---------------------------------------------------------------------------*/

static Bool_t* malloc_host_bool( size_t n, Env* env )
{
  Insist( n+1 >= 1 );
  Bool_t* result = (Bool_t*)malloc( n * sizeof(*result) );
  Insist( result );
  return result;
}

/*---------------------------------------------------------------------------*/

static P* malloc_host_P( size_t n, Env* env )
{
  Insist( n+1 >= 1 );
  P* result = _mm_malloc( n * sizeof(P), VEC_LEN * sizeof(P) );
  Insist( result );
  return result;
}

/*---------------------------------------------------------------------------*/

static P* malloc_host_pinned_P( size_t n, Env* env )
{
  return malloc_host_P( n );
}

/*---------------------------------------------------------------------------*/

static P* malloc_device_P( size_t n, Env* env )
{
  Insist( n+1 >= 1 );
  P* result = NULL;
  return result;
}

/*---------------------------------------------------------------------------*/

static void free_host_int( int* p, size_t n, Env* env )
{
  Insist( p );
  free( (void*) p );
}

/*---------------------------------------------------------------------------*/

static void free_host_bool( Bool_t* p, size_t n, Env* env )
{
  Insist( p );
  free( (void*) p );
}

/*---------------------------------------------------------------------------*/

static void free_host_P( P* p, Env* env )
{
  Insist( p );
  _mm_free( p );
}

/*---------------------------------------------------------------------------*/

static void free_host_pinned_P( P* p, Env* env )
{
  free_host_P( p );
}

/*---------------------------------------------------------------------------*/

static void free_device_P( P* p, Env* env )
{
}

#endif /*---__MIC__---*/

/*===========================================================================*/

#endif /*---_env_mic_h_---*/

/*---------------------------------------------------------------------------*/
