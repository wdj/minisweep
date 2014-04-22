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

#include <stdlib.h>

#include "types.h"
#include "env_assert.h"

/*===========================================================================*/
/*---Enums---*/

enum{ IS_USING_DEVICE = Bool_true, IS_NOT_USING_DEVICE = Bool_false };
enum{ IS_PINNED = Bool_true, IS_NOT_PINNED = Bool_false };

/*===========================================================================*/
/*---Pointer struct---*/

typedef struct
{
  size_t          n__;
  P* __restrict__ h__;
  P* __restrict__ d__;
  Bool_t          is_using_device__;
  Bool_t          is_pinned__;
} Pointer;

/*===========================================================================*/
/*---Pseudo-constructor---*/

void Pointer_ctor( Pointer* p,
                   size_t   n,
                   Bool_t   is_using_device );

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

void Pointer_delete_h( Pointer* p );

/*---------------------------------------------------------------------------*/

void Pointer_delete_d( Pointer* p );

/*===========================================================================*/
/*---Copy memory---*/

void Pointer_update_h( Pointer* p );

/*---------------------------------------------------------------------------*/

void Pointer_update_d( Pointer* p );

/*===========================================================================*/
/*---Accessors---*/

P* __restrict__ Pointer_h( Pointer* p );

/*---------------------------------------------------------------------------*/

P* __restrict__ Pointer_d( Pointer* p );

/*---------------------------------------------------------------------------*/

P* __restrict__ Pointer_a( Pointer* p );

/*===========================================================================*/

#endif /*---_common_c__pointer_h_---*/

/*---------------------------------------------------------------------------*/
