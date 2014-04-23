/*---------------------------------------------------------------------------*/
/*!
 * \file   memory.c
 * \author Wayne Joubert
 * \date   Wed Jan 15 16:06:28 EST 2014
 * \brief  Utilities for memory management.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include <stdlib.h>

#include "memory.h"
#include "types.h"
#include "env_assert.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*===========================================================================*/
/*---Allocate array of type int---*/

int* malloc_i( size_t n )
{
  int* p = (int*) malloc( n * sizeof( int ) );
  assert( p );
  return p;
}

/*===========================================================================*/
/*---Deallocate array of type int---*/

void free_i( int* p )
{
  assert( p );
  free( (void*) p );
}

/*===========================================================================*/
/*---Allocate array of type P---*/

P* malloc_P( size_t n )
{
  P* p = (P*) malloc( n * sizeof( P ) );
  assert( p );
  return p;
}

/*===========================================================================*/
/*---Deallocate array of type P---*/

void free_P( P* p )
{
  assert( p );
  free( (void*) p );
}

/*===========================================================================*/
  
#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*---------------------------------------------------------------------------*/
