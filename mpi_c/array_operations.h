/*---------------------------------------------------------------------------*/
/*!
 * \file   array_operations.h
 * \author Wayne Joubert
 * \date   Thu Jan 16 15:39:53 EST 2014
 * \brief  Functions for operating on special-purpose multidimensional arrays.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _serial_c__array_operations_h_
#define _serial_c__array_operations_h_

#include "env.h"
#include "definitions.h"
#include "dimensions.h"
#include "array_accessors.h"
#include "quantities.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*===========================================================================*/
/*---Initialize state vector to required input value---*/

static void initialize_state( P* __restrict__    v,
                              const Dimensions   dims,
                              int                nu,
                              const Quantities*  quan )
{
  int ix = 0;
  int iy = 0;
  int iz = 0;
  int ie = 0;
  int im = 0;
  int iu = 0;

  for( iz=0; iz<dims.nz; ++iz )
  for( iy=0; iy<dims.ny; ++iy )
  for( ix=0; ix<dims.nx; ++ix )
  for( ie=0; ie<dims.ne; ++ie )
  for( im=0; im<dims.nm; ++im )
  for( iu=0; iu<nu; ++iu )
  {
    *ref_state( v, dims, nu, ix, iy, iz, ie, im, iu )
                = Quantities_init_state( quan, ix, iy, iz, ie, im, iu, dims );
  }
}

/*===========================================================================*/
/*---Initialize state vector to zero---*/

static void initialize_state_zero( P* __restrict__   v,
                                   const Dimensions  dims,
                                   int               nu )
{
  size_t i = 0;
  size_t n = Dimensions_size_state( dims, nu );

  for( i=0; i<n; ++i )
  {
    v[i] = P_zero();
  }
}

/*===========================================================================*/
/*---Compute vector norm info for state vector---*/

static void get_state_norms( const P* __restrict__  vi,
                             const P* __restrict__  vo,
                             const Dimensions       dims,
                             int                    nu,
                             P*                     normsqp,
                             P*                     normsqdiffp )
{
  Assert( normsqp     != NULL ? "Null pointer encountered" : 0 );
  Assert( normsqdiffp != NULL ? "Null pointer encountered" : 0);

  int ix = 0;
  int iy = 0;
  int iz = 0;
  int ie = 0;
  int im = 0;
  int iu = 0;

  P normsq     = P_zero();
  P normsqdiff = P_zero();

  for( iz=0; iz<dims.nz; ++iz )
  for( iy=0; iy<dims.ny; ++iy )
  for( ix=0; ix<dims.nx; ++ix )
  for( ie=0; ie<dims.ne; ++ie )
  for( im=0; im<dims.nm; ++im )
  for( iu=0; iu<nu; ++iu )
  {
    const P val_vi = *const_ref_state( vi, dims, nu, ix, iy, iz, ie, im, iu );
    const P val_vo = *const_ref_state( vo, dims, nu, ix, iy, iz, ie, im, iu );
    const P diff   = val_vi - val_vo;
    normsq        += val_vo * val_vo;
    normsqdiff    += diff   * diff;
  }
  Assert( normsq     >= P_zero() );
  Assert( normsqdiff >= P_zero() );
  normsq     = Env_sum_P( normsq );
  normsqdiff = Env_sum_P( normsqdiff );

  *normsqp     = normsq;
  *normsqdiffp = normsqdiff;
}

/*===========================================================================*/
/*---Copy vector---*/

static void copy_vector(       P* __restrict__  vo,
                         const P* __restrict__  vi,
                         size_t                 n )
{
  size_t i = 0;

  for( i=0; i<n; ++i )
  {
    vo[i] = vi[i];
  }
}

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

#endif /*---_serial_c__array_operations_h_---*/

/*---------------------------------------------------------------------------*/
