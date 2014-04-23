/*---------------------------------------------------------------------------*/
/*!
 * \file   env_cuda.h
 * \author Wayne Joubert
 * \date   Tue Apr 22 17:03:08 EDT 2014
 * \brief  Environment settings for cuda.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _mpi_c__env_cuda_h_
#define _mpi_c__env_cuda_h_

#ifdef USE_CUDA
#incude "cuda.h"
#endif

#include "types.h"
#include "env_assert.h"
#include "memory.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*===========================================================================*/
/*---Initialize CUDA---*/

static void Env_cuda_initialize__( Env *env, int argc, char** argv )
{
#ifdef __CUDACC__
#endif
}

/*===========================================================================*/
/*---Finalize CUDA---*/

static void Env_cuda_finalize__( Env* env )
{
#ifdef __CUDACC__
#endif
}

/*===========================================================================*/
/*---Set values from args---*/

static void Env_cuda_set_values__( Env *env, Arguments* args )
{
#ifdef __CUDACC__
  env->is_using_device__ = Arguments_consume_int_or_default( args,
                                             "--is_using_device", Bool_false );
  Insist( env->is_using_device__ == 0 ||
          env->is_using_device__ == 1 ? "Invalid is_using_device value." : 0 );
#endif
}

/*===========================================================================*/
/*---Determine whether using device---*/

static Bool_t Env_cuda_is_using_device( Env *env )
{
#ifdef __CUDACC__
  return env->is_using_device__;
#else
  return Bool_false;
#endif
}

/*===========================================================================*/

static P* Env_cuda_malloc_P( size_t n )
{
  assert( n+1 >= 1 );

  P* result = NULL;

#ifdef __CUDACC__
  cudaMalloc( &result, n==0 ? ((size_t)1) : n*sizeof(P) );
  assert( result );
#endif

  return result;
}

/*---------------------------------------------------------------------------*/

static P* Env_cuda_malloc_host_P( size_t n )
{
  assert( n+1 >= 1 );

  P* result = NULL;

#ifdef __CUDACC__
  cudaMallocHost( &result, n==0 ? ((size_t)1) : n*sizeof(P) );
#else
   result = malloc_P( n );
#endif
  assert( result );

  return result;
}

/*---------------------------------------------------------------------------*/

static void Env_cuda_free_P( P* p )
{
#ifdef __CUDACC__
  cudaFree( p );
#endif
}

/*---------------------------------------------------------------------------*/

static void Env_cuda_free_host_P( P* p )
{
#ifdef __CUDACC__
  cudaFreeHost( p );
#else
  free_P( p );
#endif
}

/*---------------------------------------------------------------------------*/

static void Env_cuda_copy_host_to_device_P( P*     p_d,
                                            P*     p_h,
                                            size_t n )
{
#ifdef __CUDACC__
  assert( p_d );
  assert( p_h );
  assert( n+1 >= 1 );

  cudaMemcpy( p_d, p_h, n*sizeof(P), cudaMemcpyHostToDevice );
#endif
}

/*---------------------------------------------------------------------------*/

static void Env_cuda_copy_device_to_host_P( P*     p_h,
                                            P*     p_d,
                                            size_t n )
{
#ifdef __CUDACC__
  assert( p_h );
  assert( p_d );
  assert( n+1 >= 1 );

  cudaMemcpy( p_h, p_d, n*sizeof(P), cudaMemcpyDeviceToHost );
#endif
}

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

#endif /*---_mpi_c__env_cuda_h_---*/

/*---------------------------------------------------------------------------*/
