/*---------------------------------------------------------------------------*/
/*!
 * \file   pointer.h
 * \author Wayne Joubert
 * \date   Tue Apr 22 14:57:52 EDT 2014
 * \brief  Tools for handling mirrored host/device arrays.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _common_c__pointer_h_
#define _common_c__pointer_h_

#include <stddef.h>

#include "types.h"
#include "env_assert.h"

#include "pointer_kernels.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*===========================================================================*/
/*---Null object---*/

static inline Pointer Pointer_null()
{
  Pointer p;

  p.n__ = 0;
  p.h__ = NULL;
  p.d__ = NULL;
  p.is_using_device__ = Bool_false;
  p.is_pinned__       = Bool_false;
  p.is_alias__        = Bool_false;

  return p;
}

/*===========================================================================*/
/*---Pseudo-constructors---*/

void Pointer_ctor( Pointer* p,
                   size_t   n,
                   Bool_t   is_using_device );

/*---------------------------------------------------------------------------*/

void Pointer_ctor_alias( Pointer* p,
                         Pointer* source,
                         size_t   base,
                         size_t   n );

/*---------------------------------------------------------------------------*/

void Pointer_set_pinned( Pointer* p,
                         Bool_t   is_pinned );

/*===========================================================================*/
/*---Pseudo-destructor---*/

void Pointer_dtor( Pointer* p );

/*===========================================================================*/
/*---De/allocate memory---*/

void Pointer_create_h( Pointer* p );

/*---------------------------------------------------------------------------*/

void Pointer_create_d( Pointer* p );

/*---------------------------------------------------------------------------*/

void Pointer_create( Pointer* p );

/*---------------------------------------------------------------------------*/

void Pointer_delete_h( Pointer* p );

/*---------------------------------------------------------------------------*/

void Pointer_delete_d( Pointer* p );

/*---------------------------------------------------------------------------*/

void Pointer_delete( Pointer* p );

/*===========================================================================*/
/*---Copy memory---*/

void Pointer_update_h( Pointer* p );

/*---------------------------------------------------------------------------*/

void Pointer_update_d( Pointer* p );

/*---------------------------------------------------------------------------*/

void Pointer_update_h_stream( Pointer* p, Stream_t stream );

/*---------------------------------------------------------------------------*/

void Pointer_update_d_stream( Pointer* p, Stream_t stream );

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

#endif /*---_common_c__pointer_h_---*/

/*---------------------------------------------------------------------------*/
