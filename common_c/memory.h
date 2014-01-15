/*---------------------------------------------------------------------------*/
/*!
 * \file   memory.h
 * \author Wayne Joubert
 * \date   Wed Jan 15 16:06:28 EST 2014
 * \brief  Utilities for memory management.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _common_c__memory_h_
#define _common_c__memory_h_

#include "defs.h"

/*---------------------------------------------------------------------------*/
/*---Allocate array of type "P"---*/

P* pmalloc( size_t n );

/*---------------------------------------------------------------------------*/
/*---Deallocate array of type "P"---*/

void pfree( P* p );

/*---------------------------------------------------------------------------*/

#endif /*---_common_c__memory_h_---*/

/*---------------------------------------------------------------------------*/
