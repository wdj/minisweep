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

#ifdef IS_MIC

static int* malloc_host_int( size_t n, Env* env )
{
  INSIST( n+1 >= 1 );

  int* p = (int*)malloc( n * sizeof(*p) );
  INSIST( p );
  env->cpu_mem += n * sizeof(*p);
  env->cpu_mem_max = env->cpu_mem > env->cpu_mem_max ?
                     env->cpu_mem : env->cpu_mem_max;
  return p;
}

/*---------------------------------------------------------------------------*/

static Bool_t* malloc_host_bool( size_t n, Env* env )
{
  INSIST( n+1 >= 1 );

  Bool_t* p = (Bool_t*)malloc( n * sizeof(*p) );
  INSIST( p );
  env->cpu_mem += n * sizeof(*p);
  env->cpu_mem_max = env->cpu_mem > env->cpu_mem_max ?
                     env->cpu_mem : env->cpu_mem_max;
  return p;
}

/*---------------------------------------------------------------------------*/

static P* malloc_host_P( size_t n, Env* env )
{
  INSIST( n+1 >= 1 );

  P* p = _mm_malloc( n * sizeof(*p), VEC_LEN * sizeof(*p) );
  INSIST( p );
  env->cpu_mem += n * sizeof(*p);
  env->cpu_mem_max = env->cpu_mem > env->cpu_mem_max ?
                     env->cpu_mem : env->cpu_mem_max;
  return p;
}

/*---------------------------------------------------------------------------*/

static P* malloc_host_pinned_P( size_t n, Env* env )
{
  return malloc_host_P( n, env );
}

/*---------------------------------------------------------------------------*/

static P* malloc_device_P( size_t n, Env* env )
{
  INSIST( n+1 >= 1 );
  P* p = NULL;
  return p;
}

/*---------------------------------------------------------------------------*/

static void free_host_int( int* p, size_t n, Env* env )
{
  INSIST( p );

  free( (void*) p );
  env->cpu_mem -= n * sizeof(*p);
}

/*---------------------------------------------------------------------------*/

static void free_host_bool( Bool_t* p, size_t n, Env* env )
{
  INSIST( p );

  free( (void*) p );
  env->cpu_mem -= n * sizeof(*p);
}

/*---------------------------------------------------------------------------*/

static void free_host_P( P* p, size_t n, Env* env )
{
  INSIST( p );

  _mm_free( p );
  env->cpu_mem -= n * sizeof(*p);
}

/*---------------------------------------------------------------------------*/

static void free_host_pinned_P( P* p, size_t n, Env* env )
{
  free_host_P( p, n, env );
  env->cpu_mem -= n * sizeof(*p);
}

/*---------------------------------------------------------------------------*/

static void free_device_P( P* p, size_t n, Env* env )
{
}

#endif /*---IS_MIC---*/

/*===========================================================================*/

#endif /*---_env_mic_h_---*/

/*---------------------------------------------------------------------------*/
