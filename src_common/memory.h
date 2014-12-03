/*---------------------------------------------------------------------------*/
/*!
 * \file   memory.h
 * \author Wayne Joubert
 * \date   Wed Jan 15 16:06:28 EST 2014
 * \brief  Utilities for memory management.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _memory_h_
#define _memory_h_

#include <stddef.h>

#include "types.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*===========================================================================*/
/*---Allocate array of type int---*/

int* malloc_i( size_t n );

/*===========================================================================*/
/*---Deallocate array of type int---*/

void free_i( int* p );

/*===========================================================================*/
/*---Allocate array of type P---*/

P* malloc_P( size_t n );

/*===========================================================================*/
/*---Deallocate array of type P---*/

void free_P( P* p );

/*===========================================================================*/
  
#ifdef __cplusplus
} /*---extern "C"---*/
#endif

#endif /*---_memory_h_---*/

/*---------------------------------------------------------------------------*/
