/*---------------------------------------------------------------------------*/
/*!
 * \file   pointer_kernels.h
 * \author Wayne Joubert
 * \date   Tue Apr 22 14:57:52 EDT 2014
 * \brief  pointer code for device kernels.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _pointer_kernels_h_
#define _pointer_kernels_h_

#include <stddef.h>

#include "types_kernels.h"

#ifndef Assert
#define Assert(v)
#endif

#ifdef __cplusplus
extern "C"
{
#endif

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
  Bool_t          is_alias__;
} Pointer;

/*===========================================================================*/
/*---Accessors---*/

TARGET_HD static inline P* __restrict__ Pointer_h( Pointer* p )
{
  Assert( p );
  Assert( p->h__ );
  return p->h__;
}

/*---------------------------------------------------------------------------*/

TARGET_HD static inline const P* __restrict__ Pointer_const_h( const Pointer* p )
{
  Assert( p );
  Assert( p->h__ );
  return p->h__;
}

/*---------------------------------------------------------------------------*/

TARGET_HD static inline P* __restrict__ Pointer_d( Pointer* p )
{
  Assert( p );
  Assert( p->d__ );
  Assert( p->is_using_device__ );
  return p->d__;
}

/*---------------------------------------------------------------------------*/

TARGET_HD static inline const P* __restrict__ Pointer_const_d( const Pointer* p )
{
  Assert( p );
  Assert( p->d__ );
  Assert( p->is_using_device__ );
  return p->d__;
}

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

#endif /*---_pointer_kernels_h_---*/

/*---------------------------------------------------------------------------*/
