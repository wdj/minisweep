/*---------------------------------------------------------------------------*/
/*!
 * \file   env_mic.h
 * \author Wayne Joubert
 * \date   Wed Jun 11 09:33:15 EDT 2014
 * \brief  Environment settings for Intel MIC.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _mpi_c__env_mic_h_
#define _mpi_c__env_mic_h_

#include "types.h"

/*===========================================================================*/
/*---Enums---*/

#ifdef __MIC__
enum{ IS_USING_MIC = Bool_true };
#else
enum{ IS_USING_MIC = Bool_true };
#endif

#ifdef __MIC__
enum{ VEC_LEN = P_IS_DOUBLE ? 8 : 16 };
#endif

/*===========================================================================*/
/*---Define dummy function/macros for non-MIC case---*/

#ifndef __MIC__
#define __assume( a )
#define __assume_aligned( a, b )
#define __declspec( a )
#endif

/*===========================================================================*/
/*---Memory management---*/

#ifdef __MIC__

static P* malloc_host_P( size_t n )
{
  Assert( n+1 >= 1 );
  P* result = _mm_malloc( n * sizeof(P), VEC_LEN * sizeof(P) );
  Assert( result );
  return result;
}

/*---------------------------------------------------------------------------*/

static P* malloc_host_pinned_P( size_t n )
{
  Assert( n+1 >= 1 );
  P* result = _mm_malloc( n * sizeof(P), VEC_LEN * sizeof(P) );
  Assert( result );
  return result;
}

/*---------------------------------------------------------------------------*/

static P* malloc_device_P( size_t n )
{
  Assert( n+1 >= 1 );
  P* result = NULL;
  return result;
}

/*---------------------------------------------------------------------------*/

static void free_host_P( P* p )
{
  Assert( p );
  _mm_free( p );
}

/*---------------------------------------------------------------------------*/

static void free_host_pinned_P( P* p )
{
  Assert( p );
  _mm_free( p );
}

/*---------------------------------------------------------------------------*/

static void free_device_P( P* p )
{
}

#endif /*---__MIC__---*/

/*===========================================================================*/

#endif /*---_mpi_c__env_cuda_h_---*/

/*---------------------------------------------------------------------------*/
