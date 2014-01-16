/*---------------------------------------------------------------------------*/
/*!
 * \file   quantities_testing.h
 * \author Wayne Joubert
 * \date   Wed Jan 15 16:06:28 EST 2014
 * \brief  Declarations for physical quantities, testing case.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _serial_c__quantities_testing_h_
#define _serial_c__quantities_testing_h_

#include <stdio.h>
#include "dimensions.h"
#include "array_accessors.h"

/*---------------------------------------------------------------------------*/
/*---Struct to hold pointers to arrays associated with physical quantities---*/

typedef struct
{
  P* __restrict__ a_from_m;
  P* __restrict__ m_from_a;
} Quantities;

/*---------------------------------------------------------------------------*/
/*---Pseudo-constructor for Quantities struct---*/

void Quantities_ctor( Quantities* quan, Dimensions dims );

/*---------------------------------------------------------------------------*/
/*---Pseudo-destructor for Quantities struct---*/

void Quantities_dtor( Quantities* quan );

/*---------------------------------------------------------------------------*/
/*---Scale factor for energy---*/
/*---pseudo-private member function---*/

static inline int Quantities_scalefactor_energy__( int ie )
{
  /*---Power-of-two multiplier for each energy group, to help catch errors
       regarding indexing of energy groups.
  ---*/

  return 1 << ( ie & ( (1<<2) - 1 ) );
}

/*---------------------------------------------------------------------------*/
/*---Scale factor for space---*/
/*---pseudo-private member function---*/

static inline int Quantities_scalefactor_space__( int ix, int iy, int iz )
{
  /*---Red/black checkerboard pattern in space, either "1" or "2",
       to help catch spatial indexing errors.
  --*/

  return 1. + 1. * ( ( (ix+2) + (iy+2) + (iz+2) ) % 2 );
}

/*---------------------------------------------------------------------------*/
/*---Scale factor for angles---*/
/*---pseudo-private member function---*/

static inline int Quantities_scalefactor_angle__( int ia )
{
  /*---Create a "random" power of 2. Limit the size by taking only
       the low order bits of ia
  ---*/

  return 1 << ( ia & ( (1<<3) - 1 ) );
}

/*---------------------------------------------------------------------------*/
/*---Accessor for weights---*/
/*---pseudo-private member function---*/

static inline P Quantities_xfluxweight__( int ia )
{
  return (P) ( 1 / (P) 2 );
}

/*---------------------------------------------------------------------------*/
/*---Accessor for weights---*/
/*---pseudo-private member function---*/

static inline P Quantities_yfluxweight__( int ia )
{
  return (P) ( 1 / (P) 4 );
}

/*---------------------------------------------------------------------------*/
/*---Affine function used as basis for vector of linear values.
---pseudo-private member function---*/

static inline P Quantities_affinefunction__( int i )
{
  /*---NOTE: to insure that the mappings from moments to angles to moments
       restores the original vector, this must be just like this ---*/
  return (P) ( 1 + i );
}

/*---------------------------------------------------------------------------*/
/*---Accessor for weights---*/
/*---pseudo-private member function---*/

static inline P Quantities_zfluxweight__( int ia )
{
  /*---NOTE: The flux weights are calculated so that each flux "out" is
       identical to any flux "in", under assumptions on the incoming fluxes
       and the state vector at the cell.  The problem boundaries are
       initialized to this required flux value to make this work.
       Note that this value is exactly Quantities_scalefactor_angle__( ia )
       times the state vector value (in angle space) at this cell.
       Powers of 2 are used so that the divides are exact in
       floating point arithmetic.
  ---*/

  return (P) ( 1 / (P) 4 - 1 / (P) Quantities_scalefactor_angle__( ia ) );
}

/*---------------------------------------------------------------------------*/
/*---Initial values for boundary array---*/

static inline P Quantities_init_facexy(
  int ix,
  int iy,
  int iz,
  int ie,
  int ia,
  int iu)
{
  /*---NOTE: this is constructed to be affine in ia (except for scale factor)
       and independent of ix, iy, iz, to facilitate calculating the
       result analytically---*/

  return (P) Quantities_affinefunction__( ia )
           * Quantities_scalefactor_angle__( ia )
           * Quantities_scalefactor_space__( ix, iy, iz )
           * Quantities_scalefactor_energy__( ie );
}

/*---------------------------------------------------------------------------*/
/*---Initial values for boundary array---*/

static inline P Quantities_init_facexz(
  int ix,
  int iy,
  int iz,
  int ie,
  int ia,
  int iu)
{
  return (P) Quantities_affinefunction__( ia )
           * Quantities_scalefactor_angle__( ia )
           * Quantities_scalefactor_space__( ix, iy, iz )
           * Quantities_scalefactor_energy__( ie );
}

/*---------------------------------------------------------------------------*/
/*---Initial values for boundary array---*/

static inline P Quantities_init_faceyz(
  int ix,
  int iy,
  int iz,
  int ie,
  int ia,
  int iu)
{
  return (P) Quantities_affinefunction__( ia )
           * Quantities_scalefactor_angle__( ia )
           * Quantities_scalefactor_space__( ix, iy, iz )
           * Quantities_scalefactor_energy__( ie );
}

/*---------------------------------------------------------------------------*/
/*---Initial values for state vector---*/

static inline P Quantities_init_state(
  int ix,
  int iy,
  int iz,
  int ie,
  int im,
  int iu )
{
  return (P) Quantities_affinefunction__( im )
           * Quantities_scalefactor_space__( ix, iy, iz )
           * Quantities_scalefactor_energy__( ie );
}

/*---------------------------------------------------------------------------*/
/*---Perform equation solve at a gridcell---*/

static inline void Quantities_solve(
  P* __restrict__ v_local,
  P* __restrict__ facexy,
  P* __restrict__ facexz,
  P* __restrict__ faceyz,
  int ix,
  int iy,
  int iz,
  int ie,
  int ioctant,
  int my_rank,
  int my_group,
  int my_octant,
  Quantities quan,
  Dimensions dims )
{
  int iu = 0;
  int ia = 0;

  /*---Average the face values and accumulate---*/

  /*---The state value and incoming face values are first adjusted to
       normalized values by removing the spatial scaling.
       They are then combined using a weighted average chosen in a special
       way to give just the expected result.
       Finally, spatial scaling is applied to the result which is then
       stored.
  ---*/

  for( iu=0; iu<NU; ++iu )
  for( ia=0; ia<dims.na; ++ia )
  {
    const P result = (
          *ref_v_local( v_local, dims, ia, iu, my_group, my_octant, my_rank )
                             / Quantities_scalefactor_space__( ix, iy, iz )
        + *ref_facexy( facexy, dims, ix, iy, ie, ia, iu, ioctant )
                             * Quantities_xfluxweight__( ia )
                             / Quantities_scalefactor_space__( ix-idirx, iy, iz )
        + *ref_facexz( facexz, dims, ix, iz, ie, ia, iu, ioctant )
                             * Quantities_yfluxweight__( ia )
                             / Quantities_scalefactor_space__( ix, iy-idiry, iz )
        + *ref_faceyz( faceyz, dims, iy, iz, ie, ia, iu, ioctant )
                             * Quantities_zfluxweight__( ia )
                             / Quantities_scalefactor_space__( ix, iy, iz-idirz )
      )                      * Quantities_scalefactor_space__( ix, iy, iz );

    *ref_v_local( v_local, dims, ia, iu, my_group, my_octant, my_rank )
                                                                      = result;
    *ref_facexy( facexy, dims, ix, iy, ie, ia, iu, ioctant ) = result;
    *ref_facexz( facexz, dims, ix, iz, ie, ia, iu, ioctant ) = result;
    *ref_faceyz( faceyz, dims, iy, iz, ie, ia, iu, ioctant ) = result;
  }

} /*---Quantities_solve---*/

/*---------------------------------------------------------------------------*/

#endif /*---_serial_c__quantities_testing_h_---*/

/*---------------------------------------------------------------------------*/