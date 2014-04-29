/*---------------------------------------------------------------------------*/
/*!
 * \file   array_accessors.h
 * \author Wayne Joubert
 * \date   Wed Jan 15 16:06:28 EST 2014
 * \brief  Functions for referencing special-purpose multidimensional arrays.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _serial_c__array_accessors_h_
#define _serial_c__array_accessors_h_

#include "function_attributes.h"
#include "env.h"
#include "definitions.h"
#include "dimensions.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*===========================================================================*/
/*---Multidimensional indexing function---*/

TARGET_HD static inline size_t ind_state(
    Dimensions       dims,
    int              nu,
    int              ix,
    int              iy,
    int              iz,
    int              ie,
    int              im,
    int              iu )
{
  Assert( nu > 0 );
  Assert( ix >= 0 && ix < dims.nx );
  Assert( iy >= 0 && iy < dims.ny );
  Assert( iz >= 0 && iz < dims.nz );
  Assert( ie >= 0 && ie < dims.ne );
  Assert( im >= 0 && im < dims.nm );
  Assert( iu >= 0 && iu < nu );

  return  im + dims.nm * (
          iu + nu      * (
          ix + dims.nx * (
          iy + dims.ny * (
          ie + dims.ne * (
          iz + dims.nz * ( /*---NOTE: This axis MUST be slowest-varying---*/
          0 ))))));
}

/*===========================================================================*/
/*---Multidimensional array accessor function---*/

TARGET_HD static inline P* ref_state(
    P* __restrict__  v,
    Dimensions       dims,
    int              nu,
    int              ix,
    int              iy,
    int              iz,
    int              ie,
    int              im,
    int              iu )
{
  Assert( v != NULL );
  Assert( nu > 0 );
  Assert( ix >= 0 && ix < dims.nx );
  Assert( iy >= 0 && iy < dims.ny );
  Assert( iz >= 0 && iz < dims.nz );
  Assert( ie >= 0 && ie < dims.ne );
  Assert( im >= 0 && im < dims.nm );
  Assert( iu >= 0 && iu < nu );

  return & v[ ind_state( dims, nu, ix, iy, iz, ie, im, iu ) ];
}

/*===========================================================================*/
/*---Multidimensional array accessor function---*/

TARGET_HD static inline const P* const_ref_state(
    const P* __restrict__  v,
    Dimensions             dims,
    int                    nu,
    int                    ix,
    int                    iy,
    int                    iz,
    int                    ie,
    int                    im,
    int                    iu )
{
  Assert( v != NULL );
  Assert( nu > 0 );
  Assert( ix >= 0 && ix < dims.nx );
  Assert( iy >= 0 && iy < dims.ny );
  Assert( iz >= 0 && iz < dims.nz );
  Assert( ie >= 0 && ie < dims.ne );
  Assert( im >= 0 && im < dims.nm );
  Assert( iu >= 0 && iu < nu );

  return & v[ ind_state( dims, nu, ix, iy, iz, ie, im, iu ) ];
}

/*===========================================================================*/
/*---Multidimensional array accessor function---*/

TARGET_HD static inline P* ref_v_local(
    P* __restrict__  v,
    Dimensions       dims,
    int              nu,
    int              ia,
    int              iu )
{
  Assert( v != NULL );
  Assert( nu > 0 );
  Assert( ia >= 0 && ia < dims.na );
  Assert( iu >= 0 && iu < nu );

  return & v[ ia + dims.na * (
              iu + nu      * (
              0 )) ];
}

/*===========================================================================*/
/*---Multidimensional array accessor function---*/

TARGET_HD static inline const P* const_ref_v_local(
    const P* __restrict__  v,
    Dimensions             dims,
    int                    nu,
    int                    ia,
    int                    iu )
{
  Assert( v != NULL );
  Assert( nu > 0 );
  Assert( ia >= 0 && ia < dims.na );
  Assert( iu >= 0 && iu < nu );

  return & v[ ia + dims.na * (
              iu + nu      * (
              0 )) ];
}

/*===========================================================================*/
/*---Multidimensional array accessor function---*/

TARGET_HD static inline P* ref_a_from_m(
    P* __restrict__  v,
    Dimensions       dims,
    int              im,
    int              ia,
    int              octant )
{
  Assert( v != NULL );
  Assert( im >= 0 && im < dims.nm );
  Assert( ia >= 0 && ia < dims.na );
  Assert( octant >= 0 && octant < NOCTANT );

  return & v[ im     + dims.nm * (
              ia     + dims.na * (
              octant + NOCTANT * (
              0 ))) ];
}

/*===========================================================================*/
/*---Multidimensional array accessor function---*/

TARGET_HD static inline const P* const_ref_a_from_m(
    const P* __restrict__  v,
    Dimensions             dims,
    int                    im,
    int                    ia,
    int                    octant )
{
  Assert( v != NULL );
  Assert( im >= 0 && im < dims.nm );
  Assert( ia >= 0 && ia < dims.na );
  Assert( octant >= 0 && octant < NOCTANT );

  return & v[ im     + dims.nm * (
              ia     + dims.na * (
              octant + NOCTANT * (
              0 ))) ];
}

/*===========================================================================*/
/*---Multidimensional array accessor function---*/

TARGET_HD static inline P* ref_m_from_a(
    P* __restrict__  v,
    Dimensions       dims,
    int              im,
    int              ia,
    int              octant )
{
  Assert( v != NULL );
  Assert( im >= 0 && im < dims.nm );
  Assert( ia >= 0 && ia < dims.na );
  Assert( octant >= 0 && octant < NOCTANT );

  return & v[ im     + dims.nm * (
              ia     + dims.na * (
              octant + NOCTANT * (
              0 ))) ];
}

/*===========================================================================*/
/*---Multidimensional array accessor function---*/

TARGET_HD static inline const P* const_ref_m_from_a(
    const P* __restrict__  v,
    Dimensions             dims,
    int                    im,
    int                    ia,
    int                    octant )
{
  Assert( v != NULL );
  Assert( im >= 0 && im < dims.nm );
  Assert( ia >= 0 && ia < dims.na );
  Assert( octant >= 0 && octant < NOCTANT );

  return & v[ im     + dims.nm * (
              ia     + dims.na * (
              octant + NOCTANT * (
              0 ))) ];
}

/*===========================================================================*/
/*---Multidimensional array accessor function---*/

TARGET_HD static inline P* ref_facexy(
    P* __restrict__  v,
    Dimensions       dims,
    int              nu,
    int              noctant_per_block,
    int              ix,
    int              iy,
    int              ie,
    int              ia,
    int              iu,
    int              octant_in_block )
{
  Assert( v != NULL );
  Assert( nu > 0 );
  Assert( ix >= 0 && ix < dims.nx );
  Assert( iy >= 0 && iy < dims.ny );
  Assert( ie >= 0 && ie < dims.ne );
  Assert( ia >= 0 && ia < dims.na );
  Assert( iu >= 0 && iu < nu );
  Assert( octant_in_block >= 0 && octant_in_block < noctant_per_block );

  return & v[ ia + dims.na * (
              iu + nu      * (
              ix + dims.nx * (
              iy + dims.ny * (
              ie + dims.ne * (
              octant_in_block ))))) ];
}

/*===========================================================================*/
/*---Multidimensional array accessor function---*/

TARGET_HD static inline const P* const_ref_facexy(
    const P* __restrict__  v,
    Dimensions             dims,
    int                    nu,
    int                    noctant_per_block,
    int                    ix,
    int                    iy,
    int                    ie,
    int                    ia,
    int                    iu,
    int                    octant_in_block )
{
  Assert( v != NULL );
  Assert( nu > 0 );
  Assert( ix >= 0 && ix < dims.nx );
  Assert( iy >= 0 && iy < dims.ny );
  Assert( ie >= 0 && ie < dims.ne );
  Assert( ia >= 0 && ia < dims.na );
  Assert( iu >= 0 && iu < nu );
  Assert( octant_in_block >= 0 && octant_in_block < noctant_per_block );

  return & v[ ia + dims.na * (
              iu + nu      * (
              ix + dims.nx * (
              iy + dims.ny * (
              ie + dims.ne * (
              octant_in_block ))))) ];
}

/*===========================================================================*/
/*---Multidimensional array accessor function---*/

TARGET_HD static inline P* ref_facexz(
    P* __restrict__  v,
    Dimensions       dims,
    int              nu,
    int              noctant_per_block,
    int              ix,
    int              iz,
    int              ie,
    int              ia,
    int              iu,
    int              octant_in_block )
{
  Assert( v != NULL );
  Assert( nu > 0 );
  Assert( ix >= 0 && ix < dims.nx );
  Assert( iz >= 0 && iz < dims.nz );
  Assert( ie >= 0 && ie < dims.ne );
  Assert( ia >= 0 && ia < dims.na );
  Assert( iu >= 0 && iu < nu );
  Assert( octant_in_block >= 0 && octant_in_block < noctant_per_block );

  return & v[ ia + dims.na * (
              iu + nu      * (
              ix + dims.nx * (
              iz + dims.nz * (
              ie + dims.ne * (
              octant_in_block ))))) ];
}

/*===========================================================================*/
/*---Multidimensional array accessor function---*/

TARGET_HD static inline const P* const_ref_facexz(
    const P* __restrict__  v,
    Dimensions             dims,
    int                    nu,
    int                    noctant_per_block,
    int                    ix,
    int                    iz,
    int                    ie,
    int                    ia,
    int                    iu,
    int                    octant_in_block )
{
  Assert( v != NULL );
  Assert( nu > 0 );
  Assert( ix >= 0 && ix < dims.nx );
  Assert( iz >= 0 && iz < dims.nz );
  Assert( ie >= 0 && ie < dims.ne );
  Assert( ia >= 0 && ia < dims.na );
  Assert( iu >= 0 && iu < nu );
  Assert( octant_in_block >= 0 && octant_in_block < noctant_per_block );

  return & v[ ia + dims.na * (
              iu + nu      * (
              ix + dims.nx * (
              iz + dims.nz * (
              ie + dims.ne * (
              octant_in_block ))))) ];
}

/*===========================================================================*/
/*---Multidimensional array accessor function---*/

TARGET_HD static inline P* ref_faceyz(
    P* __restrict__  v,
    Dimensions       dims,
    int              nu,
    int              noctant_per_block,
    int              iy,
    int              iz,
    int              ie,
    int              ia,
    int              iu,
    int              octant_in_block )
{
  Assert( v != NULL );
  Assert( nu > 0 );
  Assert( iy >= 0 && iy < dims.ny );
  Assert( iz >= 0 && iz < dims.nz );
  Assert( ie >= 0 && ie < dims.ne );
  Assert( ia >= 0 && ia < dims.na );
  Assert( iu >= 0 && iu < nu );
  Assert( octant_in_block >= 0 && octant_in_block < noctant_per_block );

  return & v[ ia + dims.na * (
              iu + nu      * (
              iy + dims.ny * (
              iz + dims.nz * (
              ie + dims.ne * (
              octant_in_block ))))) ];
}

/*===========================================================================*/
/*---Multidimensional array accessor function---*/

TARGET_HD static inline const P* const_ref_faceyz(
    const P* __restrict__  v,
    Dimensions             dims,
    int                    nu,
    int                    noctant_per_block,
    int                    iy,
    int                    iz,
    int                    ie,
    int                    ia,
    int                    iu,
    int                    octant_in_block )
{
  Assert( v != NULL );
  Assert( nu > 0 );
  Assert( iy >= 0 && iy < dims.ny );
  Assert( iz >= 0 && iz < dims.nz );
  Assert( ie >= 0 && ie < dims.ne );
  Assert( ia >= 0 && ia < dims.na );
  Assert( iu >= 0 && iu < nu );
  Assert( octant_in_block >= 0 && octant_in_block < noctant_per_block );

  return & v[ ia + dims.na * (
              iu + nu      * (
              iy + dims.ny * (
              iz + dims.nz * (
              ie + dims.ne * (
              octant_in_block ))))) ];
}

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

#endif /*---_serial_c__array_accessors_h_---*/

/*---------------------------------------------------------------------------*/
