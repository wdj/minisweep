/*---------------------------------------------------------------------------*/
/*!
 * \file   pointer.c
 * \author Wayne Joubert
 * \date   Tue Apr 22 14:57:52 EDT 2014
 * \brief  Tools for handling mirrored host/device arrays.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include <stdlib.h>

#include "env.h"
#include "pointer.h"
#include "types.h"
#include "memory.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*===========================================================================*/
/*---Pseudo-constructors---*/

void Pointer_ctor( Pointer* p,
                   size_t   n,
                   Bool_t   is_using_device )
{
  Assert( p );
  Assert( n+1 >= 1 );

  p->h__ = NULL;
  p->d__ = NULL;
  p->n__ = n;
  p->is_using_device__ = is_using_device;
  p->is_pinned__ = Bool_false;
  p->is_alias__  = Bool_false;
}

/*---------------------------------------------------------------------------*/

void Pointer_ctor_alias( Pointer* p,
                         Pointer* source,
                         size_t   base,
                         size_t   n )
{
  Assert( p );
  Assert( source );
  Assert( base+1 >= 1 );
  Assert( n+1 >= 1 );
  Assert( base+n <= source->n__ );

  p->h__ = NULL;
  if( source->h__ )
  {
    p->h__ = source->h__ + base;
  }

  p->d__ = NULL;
  if( source->h__ )
  {
    p->d__ = source->d__ + base;
  }

  p->n__               = n;
  p->is_using_device__ = source->is_using_device__;
  p->is_pinned__       = source->is_pinned__;
  p->is_alias__        = Bool_true;
}

/*---------------------------------------------------------------------------*/

void Pointer_set_pinned( Pointer* p,
                         Bool_t   is_pinned )
{
  Assert( p );
  Assert( ! p->h__
              ? "Currently cannot change pinnedness of allocated array" : 0 );
  Assert( ! p->is_alias__ );

  p->is_pinned__ = is_pinned;
}

/*===========================================================================*/
/*---Pseudo-destructor---*/

void Pointer_dtor( Pointer* p )
{
  Assert( p );

  if( p->h__ && ! p->is_alias__ )
  {
    if( p->is_pinned__ && p->is_using_device__ )
    {
      free_host_pinned_P( p->h__ );
    }
    else
    {
      free_host_P( p->h__ );
    }
  }
  p->h__ = NULL;

  if( p->d__ && ! p->is_alias__ )
  {
    free_device_P( p->d__ );
  }
  p->d__ = NULL;

  p->n__ = 0;
  p->is_using_device__ = Bool_false;
  p->is_pinned__       = Bool_false;
}

/*===========================================================================*/
/*---De/allocate memory---*/

void Pointer_create_h( Pointer* p )
{
  Assert( p );
  Assert( ! p->is_alias__ );
  Assert( ! p->h__ );

  if( p->is_pinned__ && p->is_using_device__ )
  {
    p->h__ = malloc_host_pinned_P( p->n__ );
  }
  else
  {
    p->h__ = malloc_host_P( p->n__ );
  }
  Assert( p->h__ );
}

/*---------------------------------------------------------------------------*/

void Pointer_create_d( Pointer* p )
{
  Assert( p );
  Assert( ! p->is_alias__ );
  Assert( ! p->d__ );

  if( p->is_using_device__ )
  {
    p->d__ = malloc_device_P( p->n__ );

    Assert( p->d__ );
  }
}

/*---------------------------------------------------------------------------*/

void Pointer_create( Pointer* p )
{
  Assert( p );
  Assert( ! p->is_alias__ );

  Pointer_create_h( p );
  Pointer_create_d( p );
}

/*---------------------------------------------------------------------------*/

void Pointer_delete_h( Pointer* p )
{
  Assert( p );
  Assert( ! p->is_alias__ );
  Assert( p->h__ );

  if( p->is_pinned__ && p->is_using_device__ )
  {
    free_host_pinned_P( p->h__ );
  }
  else
  {
    free_host_P( p->h__ );
  }
  p->h__ = NULL;
}

/*---------------------------------------------------------------------------*/

void Pointer_delete_d( Pointer* p )
{
  Assert( p );
  Assert( ! p->is_alias__ );

  if( p->is_using_device__ )
  {
    Assert( p->d__ );

    free_device_P( p->d__ );

    p->d__ = NULL;
  }
}

/*---------------------------------------------------------------------------*/

void Pointer_delete( Pointer* p )
{
  Assert( p );
  Assert( ! p->is_alias__ );

  Pointer_delete_h( p );
  Pointer_delete_d( p );
}

/*===========================================================================*/
/*---Copy memory---*/

void Pointer_update_h( Pointer* p )
{
  Assert( p );

  if( p->is_using_device__ )
  {
    Env_cuda_copy_device_to_host_P( p->h__, p->d__, p->n__ );
  }
}

/*---------------------------------------------------------------------------*/

void Pointer_update_d( Pointer* p )
{
  Assert( p );

  if( p->is_using_device__ )
  {
    Env_cuda_copy_host_to_device_P( p->d__, p->h__, p->n__ );
  }
}

/*---------------------------------------------------------------------------*/

void Pointer_update_h_stream( Pointer* p, Stream_t stream )
{
  Assert( p );

  if( p->is_using_device__ )
  {
    Env_cuda_copy_device_to_host_stream_P( p->h__, p->d__, p->n__, stream );
  }
}

/*---------------------------------------------------------------------------*/

void Pointer_update_d_stream( Pointer* p, Stream_t stream )
{
  Assert( p );

  if( p->is_using_device__ )
  {
    Env_cuda_copy_host_to_device_stream_P( p->d__, p->h__, p->n__, stream );
  }
}

/*===========================================================================*/
  
#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*---------------------------------------------------------------------------*/
