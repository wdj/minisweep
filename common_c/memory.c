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
#include "env.h"
#include "definitions.h"
#include "memory.h"

/*===========================================================================*/
/*---Allocate array of type "P"---*/

P* pmalloc( size_t n )
{
  P* p = (P*) malloc( n * sizeof( P ) );
  assert( p );
  return p;
}

/*===========================================================================*/
/*---Deallocate array of type "P"---*/

void pfree( P* p )
{
  assert( p );
  free( (void*) p );
}

/*---------------------------------------------------------------------------*/
