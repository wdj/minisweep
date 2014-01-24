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

#include "dimensions.h"
#include "array_accessors.h"

/*===========================================================================*/
/*---Number of unknowns per gridcell---*/

enum{ NU = 1 };

/*===========================================================================*/
/*---Struct to hold pointers to arrays associated with physical quantities---*/

typedef struct
{
  P* __restrict__ a_from_m;
  P* __restrict__ m_from_a;
} Quantities;

/*===========================================================================*/
/*---Pseudo-constructor for Quantities struct---*/

void Quantities_ctor( Quantities* quan, Dimensions dims );

/*===========================================================================*/
/*---Pseudo-destructor for Quantities struct---*/

void Quantities_dtor( Quantities* quan );

/*===========================================================================*/
/*---Flops cost of solve per element---*/

double Quantities_flops_per_solve( const Dimensions dims );

/*===========================================================================*/
/*---Scale factor for energy---*/
/*---pseudo-private member function---*/

static inline int Quantities_scalefactor_energy__( int ie )
{
  /*---Power-of-two multiplier for each energy group, to help catch errors
       regarding indexing of energy groups.
  ---*/
  assert( ie >= 0 );

  return 1 << ( ie & ( (1<<2) - 1 ) );
}

/*===========================================================================*/
/*---Scale factor for space---*/
/*---pseudo-private member function---*/

static inline int Quantities_scalefactor_space__( int ix, int iy, int iz )
{
  /*---Red/black checkerboard pattern in space, either "1" or "2",
       to help catch spatial indexing errors.
  --*/
  assert( ix >= -1 );
  assert( iy >= -1 );
  assert( iz >= -1 );

  return 1. + 1. * ( ( (ix+2) + (iy+2) + (iz+2) ) % 2 );
}

/*===========================================================================*/
/*---Scale factor for angles---*/
/*---pseudo-private member function---*/

static inline int Quantities_scalefactor_angle__( int ia )
{
  /*---Create a "random" power of 2. Limit the size by taking only
       the low order bits of ia
  ---*/
  assert( ia >= 0 );

  return 1 << ( ia & ( (1<<3) - 1 ) );
}

/*===========================================================================*/
/*---Accessor for weights---*/
/*---pseudo-private member function---*/

static inline P Quantities_xfluxweight__( int ia )
{
  assert( ia >= 0 );

  return (P) ( 1 / (P) 2 );
}

/*===========================================================================*/
/*---Accessor for weights---*/
/*---pseudo-private member function---*/

static inline P Quantities_yfluxweight__( int ia )
{
  assert( ia >= 0 );

  return (P) ( 1 / (P) 4 );
}

/*===========================================================================*/
/*---Affine function used as basis for vector of linear values.
---pseudo-private member function---*/

static inline P Quantities_affinefunction__( int i )
{
  /*---NOTE: to insure that the mappings from moments to angles to moments
       restores the original vector, this must be just like this ---*/
  return (P) ( 1 + i );
}

/*===========================================================================*/
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
  assert( ia >= 0 );

  return (P) ( 1 / (P) 4 - 1 / (P) Quantities_scalefactor_angle__( ia ) );
}

/*===========================================================================*/
/*---Initial values for boundary array---*/

static inline P Quantities_init_facexy(
  int        ix,
  int        iy,
  int        iz,
  int        ie,
  int        ia,
  int        iu,
  Dimensions dims )
{
  assert( ix >=  0 && ix < dims.nx );
  assert( iy >=  0 && iy < dims.ny );
  assert( iz >= -1 && iz < dims.nz+1 );
  assert( ie >=  0 && ie < dims.ne );
  assert( ia >=  0 && ia < dims.na );
  assert( iu >=  0 && iu < NU );

  /*---NOTE: this is constructed to be affine in ia (except for scale factor)
       and independent of ix, iy, iz, to facilitate calculating the
       result analytically---*/

  return (P) Quantities_affinefunction__( ia )
           * Quantities_scalefactor_angle__( ia )
           * Quantities_scalefactor_space__( ix, iy, iz )
           * Quantities_scalefactor_energy__( ie );
}

/*===========================================================================*/
/*---Initial values for boundary array---*/

static inline P Quantities_init_facexz(
  int        ix,
  int        iy,
  int        iz,
  int        ie,
  int        ia,
  int        iu,
  Dimensions dims )
{
  assert( ix >=  0 && ix < dims.nx );
  assert( iy >= -1 && iy < dims.ny+1 );
  assert( iz >=  0 && iz < dims.nz );
  assert( ie >=  0 && ie < dims.ne );
  assert( ia >=  0 && ia < dims.na );
  assert( iu >=  0 && iu < NU );

  return (P) Quantities_affinefunction__( ia )
           * Quantities_scalefactor_angle__( ia )
           * Quantities_scalefactor_space__( ix, iy, iz )
           * Quantities_scalefactor_energy__( ie );
}

/*===========================================================================*/
/*---Initial values for boundary array---*/

static inline P Quantities_init_faceyz(
  int        ix,
  int        iy,
  int        iz,
  int        ie,
  int        ia,
  int        iu,
  Dimensions dims )
{
  assert( ix >= -1 && ix < dims.nx+1 );
  assert( iy >=  0 && iy < dims.ny );
  assert( iz >=  0 && iz < dims.nz );
  assert( ie >=  0 && ie < dims.ne );
  assert( ia >=  0 && ia < dims.na );
  assert( iu >=  0 && iu < NU );

  return (P) Quantities_affinefunction__( ia )
           * Quantities_scalefactor_angle__( ia )
           * Quantities_scalefactor_space__( ix, iy, iz )
           * Quantities_scalefactor_energy__( ie );
}

/*===========================================================================*/
/*---Initial values for state vector---*/

static inline P Quantities_init_state(
  int        ix,
  int        iy,
  int        iz,
  int        ie,
  int        im,
  int        iu,
  Dimensions dims )
{
  assert( ix >= 0 && ix < dims.nx );
  assert( iy >= 0 && iy < dims.ny );
  assert( iz >= 0 && iz < dims.nz );
  assert( ie >= 0 && ie < dims.ne );
  assert( im >= 0 && im < dims.nm );
  assert( iu >= 0 && iu < NU );

  return (P) Quantities_affinefunction__( im )
           * Quantities_scalefactor_space__( ix, iy, iz )
           * Quantities_scalefactor_energy__( ie );
}

/*===========================================================================*/
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
  int octant,
  int octant_ind,
  Quantities quan,
  Dimensions dims )
{
  assert( v_local );
  assert( facexy );
  assert( facexz );
  assert( faceyz );
  assert( ix >= 0 && ix < dims.nx );
  assert( iy >= 0 && iy < dims.ny );
  assert( iz >= 0 && iz < dims.nz );
  assert( ie >= 0 && ie < dims.ne );
  assert( octant >= 0 && octant < NOCTANT );
  /*---NOTE: the following may not be tight---*/
  assert( octant_ind >= 0 && octant_ind < NOCTANT );

  int iu = 0;
  int ia = 0;

  int idirx = Dir_x( octant );
  int idiry = Dir_y( octant );
  int idirz = Dir_z( octant );

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
          *ref_v_local( v_local, dims, NU, ia, iu )
                  / Quantities_scalefactor_space__( ix, iy, iz )
        + *ref_facexy( facexy, dims, NU, ix, iy, ie, ia, iu, octant_ind )
                  * Quantities_xfluxweight__( ia )
                  / Quantities_scalefactor_space__( ix-Dir_inc(idirx), iy, iz )
        + *ref_facexz( facexz, dims, NU, ix, iz, ie, ia, iu, octant_ind )
                  * Quantities_yfluxweight__( ia )
                  / Quantities_scalefactor_space__( ix, iy-Dir_inc(idiry), iz )
        + *ref_faceyz( faceyz, dims, NU, iy, iz, ie, ia, iu, octant_ind )
                  * Quantities_zfluxweight__( ia )
                  / Quantities_scalefactor_space__( ix, iy, iz-Dir_inc(idirz) )
      )           * Quantities_scalefactor_space__( ix, iy, iz );

    *ref_v_local( v_local, dims, NU, ia, iu ) = result;
    *ref_facexy( facexy, dims, NU, ix, iy, ie, ia, iu, octant_ind ) = result;
    *ref_facexz( facexz, dims, NU, ix, iz, ie, ia, iu, octant_ind ) = result;
    *ref_faceyz( faceyz, dims, NU, iy, iz, ie, ia, iu, octant_ind ) = result;
  }

} /*---Quantities_solve---*/

/*===========================================================================*/

#endif /*---_serial_c__quantities_testing_h_---*/

/*---------------------------------------------------------------------------*/
