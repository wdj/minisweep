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
/*---Multidimensional indexing function---*/

TARGET_HD static inline size_t ind_state_flat(
    int              dims_nx,
    int              dims_ny,
    int              dims_nz,
    int              dims_ne,
    int              dims_nm,
    int              nu,
    int              ix,
    int              iy,
    int              iz,
    int              ie,
    int              im,
    int              iu )
{
  Assert( nu > 0 );
  Assert( ix >= 0 && ix < dims_nx );
  Assert( iy >= 0 && iy < dims_ny );
  Assert( iz >= 0 && iz < dims_nz );
  Assert( ie >= 0 && ie < dims_ne );
  Assert( im >= 0 && im < dims_nm );
  Assert( iu >= 0 && iu < nu );

  return  im + dims_nm * (
          iu + nu      * (
          ix + dims_nx * (
          iy + dims_ny * (
          ie + dims_ne * (
          iz + dims_nz * ( /*---NOTE: This axis MUST be slowest-varying---*/
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

TARGET_HD static inline const P* const_ref_state_flat(
    const P* __restrict__  v,
    int                    dims_nx,
    int                    dims_ny,
    int                    dims_nz,
    int                    dims_ne,
    int                    dims_nm,
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
  Assert( ix >= 0 && ix < dims_nx );
  Assert( iy >= 0 && iy < dims_ny );
  Assert( iz >= 0 && iz < dims_nz );
  Assert( ie >= 0 && ie < dims_ne );
  Assert( im >= 0 && im < dims_nm );
  Assert( iu >= 0 && iu < nu );

  return & v[ ind_state_flat( dims_nx, dims_ny, dims_nz, dims_ne, dims_nm, nu,
                              ix, iy, iz, ie, im, iu ) ];
}

/*===========================================================================*/
/*---Multidimensional array accessor function---*/

TARGET_HD static inline P* ref_vilocal(
    P* const __restrict__  v,
    const Dimensions       dims,
    const int              nu,
    const int              immax,
    const int              im,
    const int              iu )
{
  Assert( v != NULL );
  Assert( nu > 0 );
  Assert( im >= 0 && im < immax );
  Assert( iu >= 0 && iu < nu );

  return & v[ im + immax * (
              iu + nu    * (
              0 )) ];
}

/*===========================================================================*/
/*---Multidimensional array accessor function---*/

TARGET_HD static inline const P* const_ref_vilocal(
    const P* const __restrict__  v,
    const Dimensions             dims,
    const int                    nu,
    const int                    immax,
    const int                    im,
    const int                    iu )
{
  Assert( v != NULL );
  Assert( nu > 0 );
  Assert( im >= 0 && im < immax );
  Assert( iu >= 0 && iu < nu );

  return & v[ im + immax * (
              iu + nu    * (
              0 )) ];
}

/*===========================================================================*/
/*---Multidimensional array indexing function---*/

TARGET_HD static inline int ind_vslocal(
    const Dimensions dims,
    const int        nu,
    const int        iamax,
    const int        ia,
    const int        iu )
{
  Assert( nu > 0 );
  Assert( ia >= 0 && ia < iamax );
  Assert( iu >= 0 && iu < nu );

  return ia + iamax * (
         iu + nu    * (
         0 ));
}

/*===========================================================================*/
/*---Multidimensional array accessor function---*/

TARGET_HD static inline P* ref_vslocal(
    P* const __restrict__  v,
    const Dimensions       dims,
    const int              nu,
    const int              iamax,
    const int              ia,
    const int              iu )
{
  Assert( v != NULL );
  Assert( nu > 0 );
  Assert( ia >= 0 && ia < iamax );
  Assert( iu >= 0 && iu < nu );

  return & v[ ia + iamax * (
              iu + nu    * (
              0 )) ];
}

/*===========================================================================*/
/*---Multidimensional array accessor function---*/

TARGET_HD static inline const P* const_ref_vslocal(
    const P* const __restrict__  v,
    const Dimensions             dims,
    const int                    nu,
    const int                    iamax,
    const int                    ia,
    const int                    iu )
{
  Assert( v != NULL );
  Assert( nu > 0 );
  Assert( ia >= 0 && ia < iamax );
  Assert( iu >= 0 && iu < nu );

  return & v[ ia + iamax * (
              iu + nu    * (
              0 )) ];
}

/*===========================================================================*/
/*---Multidimensional array accessor function---*/

TARGET_HD static inline P* ref_volocal(
    P* const __restrict__  v,
    const Dimensions       dims,
    const int              nu,
    const int              immax,
    const int              im,
    const int              iu )
{
  Assert( v != NULL );
  Assert( nu > 0 );
  Assert( im >= 0 && im < immax );
  Assert( iu >= 0 && iu < nu );

  return & v[ im + immax * (
              iu + nu    * (
              0 )) ];
}

/*===========================================================================*/
/*---Multidimensional array accessor function---*/

TARGET_HD static inline P* ref_a_from_m(
    P* const __restrict__  v,
    const Dimensions       dims,
    const int              im,
    const int              ia,
    const int              octant )
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
    const P* const __restrict__  v,
    const Dimensions             dims,
    const int                    im,
    const int                    ia,
    const int                    octant )
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

TARGET_HD static inline const P* const_ref_a_from_m_flat(
    const P* const __restrict__  v,
    const int                    dims_nm,
    const int                    dims_na,
    const int                    im,
    const int                    ia,
    const int                    octant )
{
  Assert( v != NULL );
  Assert( im >= 0 && im < dims_nm );
  Assert( ia >= 0 && ia < dims_na );
  Assert( octant >= 0 && octant < NOCTANT );

  return & v[ im     + dims_nm * (
              ia     + dims_na * (
              octant + NOCTANT * (
              0 ))) ];
}

/*===========================================================================*/
/*---Multidimensional array accessor function---*/

TARGET_HD static inline P* ref_m_from_a(
    P* const __restrict__  v,
    const Dimensions       dims,
    const int              im,
    const int              ia,
    const int              octant )
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
    const P* const __restrict__  v,
    const Dimensions             dims,
    const int                    im,
    const int                    ia,
    const int                    octant )
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
/*---Multidimensional array indexing function---*/

TARGET_HD static inline int ind_m_from_a_flat(
    const int dims_nm,
    const int dims_na,
    const int im,
    const int ia,
    const int octant )
{
  Assert( im >= 0 && im < dims_nm );
  Assert( ia >= 0 && ia < dims_na );
  Assert( octant >= 0 && octant < NOCTANT );

  return im     + dims_nm * (
         ia     + dims_na * (
         octant + NOCTANT * (
         0 )));
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
