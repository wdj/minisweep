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

#include "defs.h"

/*---------------------------------------------------------------------------*/
/*---Allocate array of type "P"---*/

P* pmalloc( size_t n )
{
  return (P*) malloc( n * sizeof( P ) );
}

/*---------------------------------------------------------------------------*/
/*---Deallocate array of type "P"---*/

void pfree( P* p )
{
  free( (void*) p );
}

/*---------------------------------------------------------------------------*/
