/*---------------------------------------------------------------------------*/
/*!
 * \file   pointer.c
 * \author Wayne Joubert
 * \date   Tue Apr 22 14:57:52 EDT 2014
 * \brief  Utilities for memory management.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include <stdlib.h>

#include "pointer.h"
#include "types.h"
#include "env_assert.h"
#include "memory.h"

/*===========================================================================*/
/*---Pseudo-constructor---*/

void Pointer_ctor( Pointer* p,
                   size_t   n,
                   Bool_t   is_using_device )
{
  assert( p );
  p->h__ = NULL;
  p->d__ = NULL;
  p->n__ = n;
  p->is_using_device__ = is_using_device;
  p->is_pinned__ = Bool_false;
}

/*---------------------------------------------------------------------------*/

void Pointer_set_pinned( Pointer* p,
                         Bool_t   is_pinned )
{
  assert( p );
  assert( ! p->h__ );
  p->is_pinned__ = is_pinned;
}

/*===========================================================================*/
/*---Pseudo-destructor---*/

void Pointer_dtor( Pointer* p )
{
  assert( p );

  if( p->h__ )
  {
    if( p->is_pinned__ )
    {
    }
    else
    {
      free_P( p->h__ );
    }
  }

  if( p->d__ )
  {
  }

  p->h__ = NULL;
  p->d__ = NULL;
  p->n__ = -1;
  p->is_using_device__ = Bool_false;
  p->is_pinned__       = Bool_false;
}

/*===========================================================================*/
/*---De/allocate memory---*/

void Pointer_create_h( Pointer* p )
{
  assert( p );
  assert( ! p->h__ );

  if( p->is_pinned__ )
  {
  }
  else
  {
    p->h__ = malloc_P( p->n__ );
    assert( p->h__ );
  }
}

/*---------------------------------------------------------------------------*/

void Pointer_create_d( Pointer* p )
{
  assert( p );
  assert( ! p->d__ );



  assert( p->d__ );
}

/*---------------------------------------------------------------------------*/

void Pointer_delete_h( Pointer* p )
{
  assert( p );
  assert( p->h__ );

  if( p->is_pinned__ )
  {
  }
  else
  {
    free_P( p->h__ );
    p->h__ = NULL;
  }
}

/*---------------------------------------------------------------------------*/

void Pointer_delete_d( Pointer* p )
{
  assert( p );
  assert( p->d__ );



  p->h__ = NULL;
}

/*===========================================================================*/
/*---Copy memory---*/

void Pointer_update_h( Pointer* p )
{
}

/*---------------------------------------------------------------------------*/

void Pointer_update_d( Pointer* p )
{
}

/*===========================================================================*/
/*---Accessors---*/

P* __restrict__ Pointer_h( Pointer* p )
{
  assert( p );
  assert( p->h__ );
  return p->h__;
}

/*---------------------------------------------------------------------------*/

P* __restrict__ Pointer_d( Pointer* p )
{
  assert( p );
  assert( p->d__ );
  assert( p->is_using_device__ );
  return p->d__;
}

/*---------------------------------------------------------------------------*/

P* __restrict__ Pointer_a( Pointer* p )
{
  assert( p );
  P* __restrict__ result = p->is_using_device__ ? p->d__ : p->h__;
  assert( result );
  return result;
}

/*---------------------------------------------------------------------------*/
