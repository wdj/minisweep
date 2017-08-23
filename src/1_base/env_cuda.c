/*---------------------------------------------------------------------------*/
/*!
 * \file   env_cuda.c
 * \author Wayne Joubert
 * \date   Tue Apr 22 17:03:08 EDT 2014
 * \brief  Environment settings for cuda.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>

#include "types.h"
#include "env_types.h"
#include "env_assert.h"
#include "arguments.h"
#include "env_cuda.h"

/*===========================================================================*/
/*---Error handling---*/

Bool_t Env_cuda_last_call_succeeded()
{
  Bool_t result = Bool_true;

#ifdef USE_CUDA
  /*---NOTE: this read of the last error is a destructive read---*/
  cudaError_t error = cudaGetLastError();

  if ( error != cudaSuccess )
  {
      result = Bool_false;
      printf( "CUDA error detected: %s\n", cudaGetErrorString( error ) );
  }
#endif

  return result;
}

/*===========================================================================*/
/*---Initialize CUDA---*/

void Env_cuda_initialize_( Env *env, int argc, char** argv )
{
#ifdef USE_CUDA
  cudaStreamCreate( & env->stream_send_block_ );
  INSIST( Env_cuda_last_call_succeeded() );

  cudaStreamCreate( & env->stream_recv_block_ );
  INSIST( Env_cuda_last_call_succeeded() );

  cudaStreamCreate( & env->stream_kernel_faces_ );
  INSIST( Env_cuda_last_call_succeeded() );
#endif
}

/*===========================================================================*/
/*---Finalize CUDA---*/

void Env_cuda_finalize_( Env* env )
{
#ifdef USE_CUDA
  cudaStreamDestroy( env->stream_send_block_ );
  INSIST( Env_cuda_last_call_succeeded() );

  cudaStreamDestroy( env->stream_recv_block_ );
  INSIST( Env_cuda_last_call_succeeded() );

  cudaStreamDestroy( env->stream_kernel_faces_ );
  INSIST( Env_cuda_last_call_succeeded() );
#endif
}

/*===========================================================================*/
/*---Set values from args---*/

void Env_cuda_set_values_( Env *env, Arguments* args )
{
#ifdef USE_CUDA
  env->is_using_device_ = Arguments_consume_int_or_default( args,
                                             "--is_using_device", Bool_false );
  INSIST_UI(
          env->is_using_device_ == 0 ||
          env->is_using_device_ == 1 ? "Invalid is_using_device value." : 0 );
#endif
}

/*===========================================================================*/
/*---Determine whether using device---*/

Bool_t Env_cuda_is_using_device( const Env* const env )
{
#ifdef USE_CUDA
  return env->is_using_device_;
#else
  return Bool_false;
#endif
}

/*===========================================================================*/
/*---Memory management, for CUDA and all platforms ex. MIC---*/

#ifndef IS_MIC

int* malloc_host_int( size_t n, Env* env )
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

Bool_t* malloc_host_bool( size_t n, Env* env )
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

P* malloc_host_P( size_t n, Env* env )
{
  INSIST( n+1 >= 1 );

  P* p = (P*)malloc( n * sizeof(*p) );
  INSIST( p );
  env->cpu_mem += n * sizeof(*p);
  env->cpu_mem_max = env->cpu_mem > env->cpu_mem_max ?
                     env->cpu_mem : env->cpu_mem_max;
  return p;
}

/*---------------------------------------------------------------------------*/

P* malloc_host_pinned_P( size_t n, Env* env )
{
  INSIST( n+1 >= 1 );

  P* p = NULL;
#ifdef USE_CUDA
  cudaMallocHost( &p, n==0 ? ((size_t)1) : n*sizeof(*p) );
  INSIST( Env_cuda_last_call_succeeded() );
#else
  p = (P*)malloc( n * sizeof(P) );
#endif
  INSIST( p );
  env->cpu_mem += n * sizeof(*p);
  env->cpu_mem_max = env->cpu_mem > env->cpu_mem_max ?
                     env->cpu_mem : env->cpu_mem_max;
  return p;
}

/*---------------------------------------------------------------------------*/

P* malloc_device_P( size_t n, Env* env )
{
  INSIST( n+1 >= 1 );

  P* p = NULL;
#ifdef USE_CUDA
  cudaMalloc( &p, n==0 ? ((size_t)1) : n*sizeof(*p) );
  INSIST( Env_cuda_last_call_succeeded() );
  INSIST( p );
  env->gpu_mem += n * sizeof(*p);
  env->gpu_mem_max = env->gpu_mem > env->gpu_mem_max ?
                     env->gpu_mem : env->gpu_mem_max;
#endif
  return p;
}

/*---------------------------------------------------------------------------*/

void free_host_int( int* p, size_t n, Env* env )
{
  INSIST( p );
  INSIST( n+1 >= 1 );

  free( (void*) p );
  env->cpu_mem -= n * sizeof(*p);
}

/*---------------------------------------------------------------------------*/

void free_host_bool( Bool_t* p, size_t n, Env* env )
{
  INSIST( p );
  INSIST( n+1 >= 1 );

  free( (void*) p );
  env->cpu_mem -= n * sizeof(*p);
}

/*---------------------------------------------------------------------------*/

void free_host_P( P* p, size_t n, Env* env )
{
  INSIST( p );
  INSIST( n+1 >= 1 );

  free( (void*) p );
  env->cpu_mem -= n * sizeof(*p);
}

/*---------------------------------------------------------------------------*/

void free_host_pinned_P( P* p, size_t n, Env* env )
{
  INSIST( p );
  INSIST( n+1 >= 1 );

#ifdef USE_CUDA
  cudaFreeHost( p );
  INSIST( Env_cuda_last_call_succeeded() );
#else
  free( (void*) p );
#endif
  env->cpu_mem -= n * sizeof(*p);
}

/*---------------------------------------------------------------------------*/

void free_device_P( P* p, size_t n, Env* env )
{
  INSIST( p );
  INSIST( n+1 >= 1 );

#ifdef USE_CUDA
  cudaFree( p );
  INSIST( Env_cuda_last_call_succeeded() );
  env->gpu_mem -= n * sizeof(*p);
#endif
}

#endif /*---IS_MIC---*/

/*---------------------------------------------------------------------------*/

void cuda_copy_host_to_device_P( P*     p_d,
                                 P*     p_h,
                                 size_t n )
{
#ifdef USE_CUDA
  INSIST( p_d && p_h );
  INSIST( n+1 >= 1 );

  cudaMemcpy( p_d, p_h, n*sizeof(P), cudaMemcpyHostToDevice );
  INSIST( Env_cuda_last_call_succeeded() );
#endif
}

/*---------------------------------------------------------------------------*/

void cuda_copy_device_to_host_P( P*     p_h,
                                 P*     p_d,
                                 size_t n )
{
#ifdef USE_CUDA
  INSIST( p_d && p_h );
  INSIST( n+1 >= 1 );

  cudaMemcpy( p_h, p_d, n*sizeof(P), cudaMemcpyDeviceToHost );
  INSIST( Env_cuda_last_call_succeeded() );
#endif
}

/*---------------------------------------------------------------------------*/

void cuda_copy_host_to_device_stream_P( P*       p_d,
                                        P*       p_h,
                                        size_t   n,
                                        Stream_t stream )
{
#ifdef USE_CUDA
  INSIST( p_d && p_h );
  INSIST( n+1 >= 1 );

  cudaMemcpyAsync( p_d, p_h, n*sizeof(P), cudaMemcpyHostToDevice, stream );
  INSIST( Env_cuda_last_call_succeeded() );
#endif
}

/*---------------------------------------------------------------------------*/

void cuda_copy_device_to_host_stream_P( P*       p_h,
                                        P*       p_d,
                                        size_t   n,
                                        Stream_t stream )
{
#ifdef USE_CUDA
  INSIST( p_d && p_h );
  INSIST( n+1 >= 1 );

  cudaMemcpyAsync( p_h, p_d, n*sizeof(P), cudaMemcpyDeviceToHost, stream );
  INSIST( Env_cuda_last_call_succeeded() );
#endif
}

/*===========================================================================*/
/*---Stream management---*/

Stream_t Env_cuda_stream_send_block( Env* env )
{
#ifdef USE_CUDA
  return env->stream_send_block_;
#else
  return 0;
#endif
}

/*---------------------------------------------------------------------------*/

Stream_t Env_cuda_stream_recv_block( Env* env )
{
#ifdef USE_CUDA
  return env->stream_recv_block_;
#else
  return 0;
#endif
}

/*---------------------------------------------------------------------------*/

Stream_t Env_cuda_stream_kernel_faces( Env* env )
{
#ifdef USE_CUDA
  return env->stream_kernel_faces_;
#else
  return 0;
#endif
}

/*---------------------------------------------------------------------------*/

void Env_cuda_stream_wait( Env* env, Stream_t stream )
{
#ifdef USE_CUDA
  cudaStreamSynchronize( stream );
  INSIST( Env_cuda_last_call_succeeded() );
#endif
}

/*===========================================================================*/
