/*---------------------------------------------------------------------------*/
/*!
 * \file   array_accessors_kernels.h
 * \author Wayne Joubert
 * \date   Wed Jan 15 16:06:28 EST 2014
 * \brief  Functions to reference multi-dim arrays, code for comp. kernel.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _array_accessors_kernels_h_
#define _array_accessors_kernels_h_

#include <stddef.h>

#include "env_assert_kernels.h"

#include "types_kernels.h"
#include "definitions_kernels.h"
#include "dimensions_kernels.h"

/*===========================================================================*/
/*---Multidimensional indexing function---*/

TARGET_HD static inline size_t ind_state_flat(
    const int dims_ncell_x,
    const int dims_ncell_y,
    const int dims_ncell_z,
    const int dims_ne,
    const int dims_nm,
    const int nu,
    const int ix,
    const int iy,
    const int iz,
    const int ie,
    const int im,
    const int iu )
{
  ASSERT( nu > 0 );
  ASSERT( ix >= 0 && ix < dims_ncell_x );
  ASSERT( iy >= 0 && iy < dims_ncell_y );
  ASSERT( iz >= 0 && iz < dims_ncell_z );
  ASSERT( ie >= 0 && ie < dims_ne );
  ASSERT( im >= 0 && im < dims_nm );
  ASSERT( iu >= 0 && iu < nu );

  return  im + dims_nm      * (
          iu + nu           * (
          ix + dims_ncell_x * (
          iy + dims_ncell_y * (
          ie + dims_ne      * (
          iz + dims_ncell_z * ( /*---NOTE: This axis MUST be slowest-varying---*/
          0 ))))));
}

/*===========================================================================*/
/*---Multidimensional indexing function---*/

TARGET_HD static inline size_t ind_state(
    const Dimensions dims,
    const int        nu,
    const int        ix,
    const int        iy,
    const int        iz,
    const int        ie,
    const int        im,
    const int        iu )
{
  ASSERT( nu > 0 );
  ASSERT( ix >= 0 && ix < dims.ncell_x );
  ASSERT( iy >= 0 && iy < dims.ncell_y );
  ASSERT( iz >= 0 && iz < dims.ncell_z );
  ASSERT( ie >= 0 && ie < dims.ne );
  ASSERT( im >= 0 && im < dims.nm );
  ASSERT( iu >= 0 && iu < nu );

  return  im + dims.nm      * (
          iu + nu           * (
          ix + dims.ncell_x * (
          iy + dims.ncell_y * (
          ie + dims.ne      * (
          iz + dims.ncell_z * ( /*---NOTE: This axis MUST be slowest-varying---*/
          0 ))))));
}

/*===========================================================================*/
/*---Multidimensional array accessor function---*/

TARGET_HD static inline P* ref_state_flat(
    P* const __restrict__  v,
    const int              dims_ncell_x,
    const int              dims_ncell_y,
    const int              dims_ncell_z,
    const int              dims_ne,
    const int              dims_nm,
    const int              nu,
    const int              ix,
    const int              iy,
    const int              iz,
    const int              ie,
    const int              im,
    const int              iu )
{
  ASSERT( v );
  ASSERT( nu > 0 );
  ASSERT( ix >= 0 && ix < dims_ncell_x );
  ASSERT( iy >= 0 && iy < dims_ncell_y );
  ASSERT( iz >= 0 && iz < dims_ncell_z );
  ASSERT( ie >= 0 && ie < dims_ne );
  ASSERT( im >= 0 && im < dims_nm );
  ASSERT( iu >= 0 && iu < nu );

  return & v[ ind_state_flat( dims_ncell_x, dims_ncell_y, dims_ncell_z,
                              dims_ne, dims_nm, nu,
                              ix, iy, iz, ie, im, iu ) ];
}

/*===========================================================================*/
/*---Multidimensional array accessor function---*/

TARGET_HD static inline P* ref_state(
    P* const __restrict__  v,
    const Dimensions       dims,
    const int              nu,
    const int              ix,
    const int              iy,
    const int              iz,
    const int              ie,
    const int              im,
    const int              iu )
{
  ASSERT( v );
  ASSERT( nu > 0 );
  ASSERT( ix >= 0 && ix < dims.ncell_x );
  ASSERT( iy >= 0 && iy < dims.ncell_y );
  ASSERT( iz >= 0 && iz < dims.ncell_z );
  ASSERT( ie >= 0 && ie < dims.ne );
  ASSERT( im >= 0 && im < dims.nm );
  ASSERT( iu >= 0 && iu < nu );

  return & v[ ind_state( dims, nu, ix, iy, iz, ie, im, iu ) ];
}

/*===========================================================================*/
/*---Multidimensional array accessor function---*/

TARGET_HD static inline const P* const_ref_state_flat(
    const P* const __restrict__  v,
    const int                    dims_ncell_x,
    const int                    dims_ncell_y,
    const int                    dims_ncell_z,
    const int                    dims_ne,
    const int                    dims_nm,
    const int                    nu,
    const int                    ix,
    const int                    iy,
    const int                    iz,
    const int                    ie,
    const int                    im,
    const int                    iu )
{
  ASSERT( v );
  ASSERT( nu > 0 );
  ASSERT( ix >= 0 && ix < dims_ncell_x );
  ASSERT( iy >= 0 && iy < dims_ncell_y );
  ASSERT( iz >= 0 && iz < dims_ncell_z );
  ASSERT( ie >= 0 && ie < dims_ne );
  ASSERT( im >= 0 && im < dims_nm );
  ASSERT( iu >= 0 && iu < nu );

  return & v[ ind_state_flat( dims_ncell_x, dims_ncell_y, dims_ncell_z,
                              dims_ne, dims_nm, nu,
                              ix, iy, iz, ie, im, iu ) ];
}

/*===========================================================================*/
/*---Multidimensional array accessor function---*/

TARGET_HD static inline const P* const_ref_state(
    const P* const __restrict__  v,
    const Dimensions             dims,
    const int                    nu,
    const int                    ix,
    const int                    iy,
    const int                    iz,
    const int                    ie,
    const int                    im,
    const int                    iu )
{
  ASSERT( v );
  ASSERT( nu > 0 );
  ASSERT( ix >= 0 && ix < dims.ncell_x );
  ASSERT( iy >= 0 && iy < dims.ncell_y );
  ASSERT( iz >= 0 && iz < dims.ncell_z );
  ASSERT( ie >= 0 && ie < dims.ne );
  ASSERT( im >= 0 && im < dims.nm );
  ASSERT( iu >= 0 && iu < nu );

  return & v[ ind_state( dims, nu, ix, iy, iz, ie, im, iu ) ];
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
  ASSERT( v );
  ASSERT( nu > 0 );
  ASSERT( im >= 0 && im < immax );
  ASSERT( iu >= 0 && iu < nu );

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
  ASSERT( v );
  ASSERT( nu > 0 );
  ASSERT( im >= 0 && im < immax );
  ASSERT( iu >= 0 && iu < nu );

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
  ASSERT( nu > 0 );
  ASSERT( ia >= 0 && ia < iamax );
  ASSERT( iu >= 0 && iu < nu );

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
  ASSERT( v );
  ASSERT( nu > 0 );
  ASSERT( ia >= 0 && ia < iamax );
  ASSERT( iu >= 0 && iu < nu );

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
  ASSERT( v );
  ASSERT( nu > 0 );
  ASSERT( ia >= 0 && ia < iamax );
  ASSERT( iu >= 0 && iu < nu );

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
  ASSERT( v );
  ASSERT( nu > 0 );
  ASSERT( im >= 0 && im < immax );
  ASSERT( iu >= 0 && iu < nu );

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
  ASSERT( v );
  ASSERT( im >= 0 && im < dims.nm );
  ASSERT( ia >= 0 && ia < dims.na );
  ASSERT( octant >= 0 && octant < NOCTANT );

  return & v[ ia     + dims.na * (
              im     +      NM * (
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
  ASSERT( v );
  ASSERT( im >= 0 && im < dims.nm );
  ASSERT( ia >= 0 && ia < dims.na );
  ASSERT( octant >= 0 && octant < NOCTANT );

  return & v[ ia     + dims.na * (
              im     +      NM * (
              octant + NOCTANT * (
              0 ))) ];
}

/*===========================================================================*/
/*---Multidimensional array accessor function---*/

TARGET_HD static inline const P* __restrict__ const_ref_a_from_m_flat(
    const P* const __restrict__  v,
    const int                    dims_nm,
    const int                    dims_na,
    const int                    im,
    const int                    ia,
    const int                    octant )
{
  ASSERT( v );
  ASSERT( im >= 0 && im < dims_nm );
  ASSERT( ia >= 0 && ia < dims_na );
  ASSERT( octant >= 0 && octant < NOCTANT );

  return & v[ ia     + dims_na * (
              im     +      NM * (
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
  ASSERT( v );
  ASSERT( im >= 0 && im < dims.nm );
  ASSERT( ia >= 0 && ia < dims.na );
  ASSERT( octant >= 0 && octant < NOCTANT );

  return & v[ im     +      NM * (
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
  ASSERT( v );
  ASSERT( im >= 0 && im < dims.nm );
  ASSERT( ia >= 0 && ia < dims.na );
  ASSERT( octant >= 0 && octant < NOCTANT );

  return & v[ im     +      NM * (
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
  ASSERT( im >= 0 && im < dims_nm );
  ASSERT( ia >= 0 && ia < dims_na );
  ASSERT( octant >= 0 && octant < NOCTANT );

  return im     +      NM * (
         ia     + dims_na * (
         octant + NOCTANT * (
         0 )));
}

/*===========================================================================*/
/*---Multidimensional array accessor function---*/

TARGET_HD static inline P* ref_facexy(
    P* const __restrict__  v,
    const Dimensions       dims,
    const int              nu,
    const int              noctant_per_block,
    const int              ix,
    const int              iy,
    const int              ie,
    const int              ia,
    const int              iu,
    const int              octant_in_block )
{
  ASSERT( v );
  ASSERT( nu > 0 );
  ASSERT( ix >= 0 && ix < dims.ncell_x );
  ASSERT( iy >= 0 && iy < dims.ncell_y );
  ASSERT( ie >= 0 && ie < dims.ne );
  ASSERT( ia >= 0 && ia < dims.na );
  ASSERT( iu >= 0 && iu < nu );
  ASSERT( octant_in_block >= 0 && octant_in_block < noctant_per_block );

  return & v[ ia + dims.na      * (
              iu + nu           * (
              ix + dims.ncell_x * (
              iy + dims.ncell_y * (
              ie + dims.ne      * (
              octant_in_block ))))) ];
}

/*===========================================================================*/
/*---Multidimensional array accessor function---*/

TARGET_HD static inline const P* const_ref_facexy(
    const P* const __restrict__  v,
    const Dimensions             dims,
    const int                    nu,
    const int                    noctant_per_block,
    const int                    ix,
    const int                    iy,
    const int                    ie,
    const int                    ia,
    const int                    iu,
    const int                    octant_in_block )
{
  ASSERT( v );
  ASSERT( nu > 0 );
  ASSERT( ix >= 0 && ix < dims.ncell_x );
  ASSERT( iy >= 0 && iy < dims.ncell_y );
  ASSERT( ie >= 0 && ie < dims.ne );
  ASSERT( ia >= 0 && ia < dims.na );
  ASSERT( iu >= 0 && iu < nu );
  ASSERT( octant_in_block >= 0 && octant_in_block < noctant_per_block );

  return & v[ ia + dims.na      * (
              iu + nu           * (
              ix + dims.ncell_x * (
              iy + dims.ncell_y * (
              ie + dims.ne      * (
              octant_in_block  ))))) ];
}

/*===========================================================================*/
/*---Multidimensional array accessor function---*/

TARGET_HD static inline P* ref_facexz(
    P* const __restrict__  v,
    const Dimensions       dims,
    const int              nu,
    const int              noctant_per_block,
    const int              ix,
    const int              iz,
    const int              ie,
    const int              ia,
    const int              iu,
    const int              octant_in_block )
{
  ASSERT( v );
  ASSERT( nu > 0 );
  ASSERT( ix >= 0 && ix < dims.ncell_x );
  ASSERT( iz >= 0 && iz < dims.ncell_z );
  ASSERT( ie >= 0 && ie < dims.ne );
  ASSERT( ia >= 0 && ia < dims.na );
  ASSERT( iu >= 0 && iu < nu );
  ASSERT( octant_in_block >= 0 && octant_in_block < noctant_per_block );

  return & v[ ia + dims.na      * (
              iu + nu           * (
              ix + dims.ncell_x * (
              iz + dims.ncell_z * (
              ie + dims.ne      * (
              octant_in_block ))))) ];
}

/*===========================================================================*/
/*---Multidimensional array accessor function---*/

TARGET_HD static inline const P* const_ref_facexz(
    const P* const __restrict__  v,
    const Dimensions             dims,
    const int                    nu,
    const int                    noctant_per_block,
    const int                    ix,
    const int                    iz,
    const int                    ie,
    const int                    ia,
    const int                    iu,
    const int                    octant_in_block )
{
  ASSERT( v );
  ASSERT( nu > 0 );
  ASSERT( ix >= 0 && ix < dims.ncell_x );
  ASSERT( iz >= 0 && iz < dims.ncell_z );
  ASSERT( ie >= 0 && ie < dims.ne );
  ASSERT( ia >= 0 && ia < dims.na );
  ASSERT( iu >= 0 && iu < nu );
  ASSERT( octant_in_block >= 0 && octant_in_block < noctant_per_block );

  return & v[ ia + dims.na      * (
              iu + nu           * (
              ix + dims.ncell_x * (
              iz + dims.ncell_z * (
              ie + dims.ne      * (
              octant_in_block ))))) ];
}

/*===========================================================================*/
/*---Multidimensional array accessor function---*/

TARGET_HD static inline P* ref_faceyz(
    P* const __restrict__ v,
    const Dimensions      dims,
    const int             nu,
    const int             noctant_per_block,
    const int             iy,
    const int             iz,
    const int             ie,
    const int             ia,
    const int             iu,
    const int             octant_in_block )
{
  ASSERT( v );
  ASSERT( nu > 0 );
  ASSERT( iy >= 0 && iy < dims.ncell_y );
  ASSERT( iz >= 0 && iz < dims.ncell_z );
  ASSERT( ie >= 0 && ie < dims.ne );
  ASSERT( ia >= 0 && ia < dims.na );
  ASSERT( iu >= 0 && iu < nu );
  ASSERT( octant_in_block >= 0 && octant_in_block < noctant_per_block );

  return & v[ ia + dims.na      * (
              iu + nu           * (
              iy + dims.ncell_y * (
              iz + dims.ncell_z * (
              ie + dims.ne      * (
              octant_in_block ))))) ];
}

/*===========================================================================*/
/*---Multidimensional array accessor function---*/

TARGET_HD static inline const P* const_ref_faceyz(
    const P* const __restrict__  v,
    const Dimensions             dims,
    const int                    nu,
    const int                    noctant_per_block,
    const int                    iy,
    const int                    iz,
    const int                    ie,
    const int                    ia,
    const int                    iu,
    const int                    octant_in_block )
{
  ASSERT( v );
  ASSERT( nu > 0 );
  ASSERT( iy >= 0 && iy < dims.ncell_y );
  ASSERT( iz >= 0 && iz < dims.ncell_z );
  ASSERT( ie >= 0 && ie < dims.ne );
  ASSERT( ia >= 0 && ia < dims.na );
  ASSERT( iu >= 0 && iu < nu );
  ASSERT( octant_in_block >= 0 && octant_in_block < noctant_per_block );

  return & v[ ia + dims.na      * (
              iu + nu           * (
              iy + dims.ncell_y * (
              iz + dims.ncell_z * (
              ie + dims.ne      * (
              octant_in_block ))))) ];
}

/*===========================================================================*/

#endif /*---_array_accessors_kernels_h_---*/

/*---------------------------------------------------------------------------*/
