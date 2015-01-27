/*---------------------------------------------------------------------------*/
/*!
 * \file   env_cuda.h
 * \author Wayne Joubert
 * \date   Tue Apr 22 17:03:08 EDT 2014
 * \brief  Environment settings for cuda.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _env_cuda_h_
#define _env_cuda_h_

#ifdef USE_CUDA
#include "cuda.h"
#endif

#include "types.h"
#include "env_assert.h"
#include "env_declarations.h"
#include "memory.h"

#include "env_cuda_kernels.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*===========================================================================*/
/*---Error handling---*/

static Bool_t Env_cuda_last_call_succeeded()
{
  Bool_t result = Bool_true;

#ifdef __CUDACC__
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

static void Env_cuda_initialize_( Env *env, int argc, char** argv )
{
#ifdef __CUDACC__
  cudaStreamCreate( & env->stream_send_block_ );
  Assert( Env_cuda_last_call_succeeded() );
  cudaStreamCreate( & env->stream_recv_block_ );
  Assert( Env_cuda_last_call_succeeded() );
  cudaStreamCreate( & env->stream_kernel_faces_ );
  Assert( Env_cuda_last_call_succeeded() );
#endif
}

/*===========================================================================*/
/*---Finalize CUDA---*/

static void Env_cuda_finalize_( Env* env )
{
#ifdef __CUDACC__
  cudaStreamDestroy( env->stream_send_block_ );
  Assert( Env_cuda_last_call_succeeded() );
  cudaStreamDestroy( env->stream_recv_block_ );
  Assert( Env_cuda_last_call_succeeded() );
  cudaStreamDestroy( env->stream_kernel_faces_ );
  Assert( Env_cuda_last_call_succeeded() );
#endif
}

/*===========================================================================*/
/*---Set values from args---*/

static void Env_cuda_set_values_( Env *env, Arguments* args )
{
#ifdef __CUDACC__
  env->is_using_device_ = Arguments_consume_int_or_default( args,
                                             "--is_using_device", Bool_false );
  Insist( env->is_using_device_ == 0 ||
          env->is_using_device_ == 1 ? "Invalid is_using_device value." : 0 );
#endif
}

/*===========================================================================*/
/*---Determine whether using device---*/

static Bool_t Env_cuda_is_using_device( const Env *env )
{
#ifdef __CUDACC__
  return env->is_using_device_;
#else
  return Bool_false;
#endif
}

/*===========================================================================*/
/*---Memory management---*/

#ifndef __MIC__

static P* malloc_host_P( size_t n )
{
  Assert( n+1 >= 1 );
  P* result = malloc_P( n );
  Assert( result );
  return result;
}

/*---------------------------------------------------------------------------*/

static P* malloc_host_pinned_P( size_t n )
{
  Assert( n+1 >= 1 );

  P* result = NULL;

#ifdef __CUDACC__
  cudaMallocHost( &result, n==0 ? ((size_t)1) : n*sizeof(P) );
  Assert( Env_cuda_last_call_succeeded() );
#else
  result = malloc_P( n );
#endif
  Assert( result );

  return result;
}

/*---------------------------------------------------------------------------*/

static P* malloc_device_P( size_t n )
{
  Assert( n+1 >= 1 );

  P* result = NULL;

#ifdef __CUDACC__
  cudaMalloc( &result, n==0 ? ((size_t)1) : n*sizeof(P) );
  Assert( Env_cuda_last_call_succeeded() );
  Assert( result );
#endif

  return result;
}

/*---------------------------------------------------------------------------*/

static void free_host_P( P* p )
{
  Assert( p );
  free_P( p );
}

/*---------------------------------------------------------------------------*/

static void free_host_pinned_P( P* p )
{
  Assert( p );
#ifdef __CUDACC__
  cudaFreeHost( p );
  Assert( Env_cuda_last_call_succeeded() );
#else
  free_P( p );
#endif
}

/*---------------------------------------------------------------------------*/

static void free_device_P( P* p )
{
#ifdef __CUDACC__
  cudaFree( p );
  Assert( Env_cuda_last_call_succeeded() );
#endif
}

#endif /*---__MIC__---*/

/*---------------------------------------------------------------------------*/

static void Env_cuda_copy_host_to_device_P( P*     p_d,
                                            P*     p_h,
                                            size_t n )
{
#ifdef __CUDACC__
  Assert( p_d );
  Assert( p_h );
  Assert( n+1 >= 1 );

  cudaMemcpy( p_d, p_h, n*sizeof(P), cudaMemcpyHostToDevice );
  Assert( Env_cuda_last_call_succeeded() );
#endif
}

/*---------------------------------------------------------------------------*/

static void Env_cuda_copy_device_to_host_P( P*     p_h,
                                            P*     p_d,
                                            size_t n )
{
#ifdef __CUDACC__
  Assert( p_h );
  Assert( p_d );
  Assert( n+1 >= 1 );

  cudaMemcpy( p_h, p_d, n*sizeof(P), cudaMemcpyDeviceToHost );
  Assert( Env_cuda_last_call_succeeded() );
#endif
}

/*---------------------------------------------------------------------------*/

static void Env_cuda_copy_host_to_device_stream_P( P*       p_d,
                                                   P*       p_h,
                                                   size_t   n,
                                                   Stream_t stream )
{
#ifdef __CUDACC__
  Assert( p_d );
  Assert( p_h );
  Assert( n+1 >= 1 );

  cudaMemcpyAsync( p_d, p_h, n*sizeof(P), cudaMemcpyHostToDevice, stream );
  Assert( Env_cuda_last_call_succeeded() );
#endif
}

/*---------------------------------------------------------------------------*/

static void Env_cuda_copy_device_to_host_stream_P( P*       p_h,
                                                   P*       p_d,
                                                   size_t   n,
                                                   Stream_t stream )
{
#ifdef __CUDACC__
  Assert( p_h );
  Assert( p_d );
  Assert( n+1 >= 1 );

  cudaMemcpyAsync( p_h, p_d, n*sizeof(P), cudaMemcpyDeviceToHost, stream );
  Assert( Env_cuda_last_call_succeeded() );
#endif
}

/*===========================================================================*/
/*---Stream management---*/

static Stream_t Env_cuda_stream_send_block( Env* env )
{
#ifdef __CUDACC__
  return env->stream_send_block_;
#else
  return 0;
#endif
}

/*---------------------------------------------------------------------------*/

static Stream_t Env_cuda_stream_recv_block( Env* env )
{
#ifdef __CUDACC__
  return env->stream_recv_block_;
#else
  return 0;
#endif
}

/*---------------------------------------------------------------------------*/

static Stream_t Env_cuda_stream_kernel_faces( Env* env )
{
#ifdef __CUDACC__
  return env->stream_kernel_faces_;
#else
  return 0;
#endif
}

/*---------------------------------------------------------------------------*/

static void Env_cuda_stream_wait( Stream_t stream )
{
#ifdef __CUDACC__
  cudaStreamSynchronize( stream );
  Assert( Env_cuda_last_call_succeeded() );
#endif
}

/*---------------------------------------------------------------------------*/

static void Env_cuda_stream_wait_all()
{
#ifdef __CUDACC__
  cudaThreadSynchronize();
  Assert( Env_cuda_last_call_succeeded() );
#endif
}

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

#endif /*---_env_cuda_h_---*/

/*---------------------------------------------------------------------------*/
