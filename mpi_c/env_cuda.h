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
#include "cuda.h"
#endif

#include "types.h"
#include "env_assert.h"
#include "memory.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*===========================================================================*/
/*---Types---*/

#ifdef __CUDACC__
typedef cudaStream_t Stream_t;
#else
typedef int Stream_t;
#endif

/*===========================================================================*/
/*---Pointer to device shared memory---*/

#ifdef __CUDA_ARCH__
__shared__ extern char cuda_shared_memory[];
#endif

/*===========================================================================*/
/*---Initialize CUDA---*/

static void Env_cuda_initialize__( Env *env, int argc, char** argv )
{
#ifdef __CUDACC__
  cudaStreamCreate( & env->stream_send_block__ );
  cudaStreamCreate( & env->stream_recv_block__ );
  cudaStreamCreate( & env->stream_kernel_faces__ );
#endif
}

/*===========================================================================*/
/*---Finalize CUDA---*/

static void Env_cuda_finalize__( Env* env )
{
#ifdef __CUDACC__
  cudaStreamDestroy( env->stream_send_block__ );
  cudaStreamDestroy( env->stream_recv_block__ );
  cudaStreamDestroy( env->stream_kernel_faces__ );
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
/*---Memory management---*/

static P* Env_cuda_malloc_P( size_t n )
{
  Assert( n+1 >= 1 );

  P* result = NULL;

#ifdef __CUDACC__
  cudaMalloc( &result, n==0 ? ((size_t)1) : n*sizeof(P) );
  Assert( result );
#endif

  return result;
}

/*---------------------------------------------------------------------------*/

static P* Env_cuda_malloc_host_P( size_t n )
{
  Assert( n+1 >= 1 );

  P* result = NULL;

#ifdef __CUDACC__
  cudaMallocHost( &result, n==0 ? ((size_t)1) : n*sizeof(P) );
#else
   result = malloc_P( n );
#endif
  Assert( result );

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
  Assert( p_d );
  Assert( p_h );
  Assert( n+1 >= 1 );

  cudaMemcpy( p_d, p_h, n*sizeof(P), cudaMemcpyHostToDevice );
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
#endif
}

/*===========================================================================*/
/*---Device thread management---*/

TARGET_HD static int Env_cuda_threadblock( int axis )
{
  Assert( axis >= 0 && axis < 2 );

#ifdef __CUDA_ARCH__
  return axis==0 ? blockIdx.x :
         axis==1 ? blockIdx.y :
                   blockIdx.z;
#else
  return 0;
#endif
}

/*---------------------------------------------------------------------------*/

TARGET_HD static int Env_cuda_thread_in_threadblock( int axis )
{
  Assert( axis >= 0 && axis < 2 );

#ifdef __CUDA_ARCH__
  return axis==0 ? threadIdx.x :
         axis==1 ? threadIdx.y :
                   threadIdx.z;
#else
  return 0;
#endif
}

/*---------------------------------------------------------------------------*/

TARGET_HD static void Env_cuda_sync_threadblock()
{
#ifdef __CUDA_ARCH__
  __syncthreads();
#endif
}

/*===========================================================================*/
/*---Stream management---*/

static Stream_t Env_cuda_stream_send_block( Env* env )
{
#ifdef __CUDACC__
  return env->stream_send_block__;
#else
  return 0;
#endif
}

/*---------------------------------------------------------------------------*/

static Stream_t Env_cuda_stream_recv_block( Env* env )
{
#ifdef __CUDACC__
  return env->stream_recv_block__;
#else
  return 0;
#endif
}

/*---------------------------------------------------------------------------*/

static Stream_t Env_cuda_stream_kernel_faces( Env* env )
{
#ifdef __CUDACC__
  return env->stream_kernel_faces__;
#else
  return 0;
#endif
}

/*---------------------------------------------------------------------------*/

static void Env_cuda_stream_wait( Stream_t stream )
{
#ifdef __CUDACC__
  cudaStreamSynchronize( stream );
#endif
}

/*---------------------------------------------------------------------------*/

static void Env_cuda_stream_wait_all()
{
#ifdef __CUDACC__
  cudaThreadSynchronize();
#endif
}

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

#endif /*---_mpi_c__env_cuda_h_---*/

/*---------------------------------------------------------------------------*/
