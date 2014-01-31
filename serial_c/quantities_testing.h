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

enum{ NU = 4 };

/*===========================================================================*/
/*---Struct to hold pointers to arrays associated with physical quantities---*/

typedef struct
{
  P* __restrict__  a_from_m;
  P* __restrict__  m_from_a;
  int*             ix_base_vals;
  int*             iy_base_vals;
  int              ix_base;
  int              iy_base;
  int              nx_g;
  int              ny_g;
  int              nz_g;
} Quantities;

/*===========================================================================*/
/*---Pseudo-constructor for Quantities struct---*/

void Quantities_ctor( Quantities*       quan,
                      const Dimensions  dims,
                      Env*              env );

/*===========================================================================*/
/*---Pseudo-destructor for Quantities struct---*/

void Quantities_dtor( Quantities* quan );

/*===========================================================================*/
/*---Initialize Quantities a_from_m, m_from_a matrices---*/
/*---pseudo-private member function---*/

void Quantities_init_am_matrices__( Quantities*       quan,
                                    const Dimensions  dims );

/*===========================================================================*/
/*---Initialize Quantities subgrid decomp info---*/
/*---pseudo-private member function---*/

void Quantities_init_decomp__( Quantities*       quan,
                               const Dimensions  dims,
                               Env*              env );

/*===========================================================================*/
/*---Flops cost of solve per element---*/

double Quantities_flops_per_solve( const Dimensions dims );

/*===========================================================================*/
/*---Scale factor for energy---*/
/*---pseudo-private member function---*/

static inline int Quantities_scalefactor_energy__( int ie,
                                                   Dimensions dims )
{
  /*---Random power-of-two multiplier for each energy group,
       to help catch errors regarding indexing of energy groups.
  ---*/
  assert( ie >= 0 && ie < dims.ne );

  const int im = 714025;
  const int ia = 1366;
  const int ic = 150889;

  int result = ( (ie)*ia + ic ) % im;
  result = result & ( (1<<2) - 1 );
  result = 1 << result;

  return result;
}

/*===========================================================================*/
/*---Scale factor for unknown---*/
/*---pseudo-private member function---*/

static inline int Quantities_scalefactor_unknown__( int iu )
{
  /*---Random power-of-two multiplier for each cell unknown,
       to help catch errors regarding indexing of cell unknowns.
  ---*/
  assert( iu >= 0 && iu < NU );

  const int im = 312500;
  const int ia = 741;
  const int ic = 66037;

  int result = ( (iu)*ia + ic ) % im;
  result = result & ( (1<<2) - 1 );
  result = 1 << result;

  return result;
}

/*===========================================================================*/
/*---Scale factor for space---*/
/*---pseudo-private member function---*/

static inline int Quantities_scalefactor_space__( const Quantities* quan,
                                                  int ix_g,
                                                  int iy_g,
                                                  int iz_g )
{
  /*---Create power of 2 based on hash of the spatial location.
  ---*/
  assert( ix_g >= -1 && ix_g <= quan->nx_g );
  assert( iy_g >= -1 && iy_g <= quan->ny_g );
  assert( iz_g >= -1 && iz_g <= quan->nz_g );

  const int im = 134456;
  const int ia = 8121;
  const int ic = 28411;

  int result = 0;
  result = ( (result+(ix_g+2))*ia + ic ) % im;
  result = ( (result+(iy_g+2))*ia + ic ) % im;
  result = ( (result+(iz_g+2))*ia + ic ) % im;
  result = result & ( (1<<2) - 1 );
  result = 1 << result;

  return result;
}

/*===========================================================================*/
/*---Scale factor for angles---*/
/*---pseudo-private member function---*/

static inline int Quantities_scalefactor_angle__( Dimensions dims, int ia )
{
  /*---Create a "random" power of 2. Limit the size by taking only
       the low order bits of ia
  ---*/
  assert( ia >= 0 && ia < dims.na );

  return 1 << ( ia & ( (1<<3) - 1 ) );
}

/*===========================================================================*/
/*---Accessor for weights---*/
/*---pseudo-private member function---*/

static inline P Quantities_xfluxweight__( Dimensions dims, int ia )
{
  assert( ia >= 0 && ia < dims.na );

  return (P) ( 1 / (P) 2 );
}

/*===========================================================================*/
/*---Accessor for weights---*/
/*---pseudo-private member function---*/

static inline P Quantities_yfluxweight__( Dimensions dims, int ia )
{
  assert( ia >= 0 && ia < dims.na );

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

static inline P Quantities_zfluxweight__( Dimensions dims,  int ia )
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
  assert( ia >= 0 && ia < dims.na );

  return (P) ( 1 / (P) 4 - 1 / (P) Quantities_scalefactor_angle__( dims, ia ) );
}

/*===========================================================================*/
/*---Initial values for boundary array---*/

static inline P Quantities_init_facexy(
  const Quantities*  quan,
  int                ix_g,
  int                iy_g,
  int                iz_g,
  int                ie,
  int                ia,
  int                iu,
  int                octant,
  const Dimensions   dims_g )
{
  assert( ix_g >=  0 && ix_g <  dims_g.nx );
  assert( iy_g >=  0 && iy_g <  dims_g.ny );
  assert( ( iz_g == -1        && Dir_z(octant)==Dir_up() ) ||
          ( iz_g == dims_g.nz && Dir_z(octant)==Dir_dn() ) );
  assert( ie >=  0 && ie < dims_g.ne );
  assert( ia >=  0 && ia < dims_g.na );
  assert( iu >=  0 && iu < NU );
  assert( octant >= 0 && octant < NOCTANT );

  /*---NOTE: this is constructed to be affine in ia (except for scale factor)
       and independent of ix, iy, iz, to facilitate calculating the
       result analytically---*/

  return   ( (P) Quantities_affinefunction__( ia ) )
         * ( (P) Quantities_scalefactor_angle__( dims_g, ia ) )
         * ( (P) Quantities_scalefactor_space__( quan, ix_g, iy_g, iz_g ) )
         * ( (P) Quantities_scalefactor_energy__( ie, dims_g ) )
         * ( (P) Quantities_scalefactor_unknown__( iu ) );
}

/*===========================================================================*/
/*---Initial values for boundary array---*/

static inline P Quantities_init_facexz(
  const Quantities*  quan,
  int                ix_g,
  int                iy_g,
  int                iz_g,
  int                ie,
  int                ia,
  int                iu,
  int                octant,
  const Dimensions   dims_g )
{
  assert( ix_g >=  0 && ix_g < dims_g.nx );
  assert( ( iy_g == -1        && Dir_y(octant)==Dir_up() ) ||
          ( iy_g == dims_g.ny && Dir_y(octant)==Dir_dn() ) );
  assert( iz_g >=  0 && iz_g < dims_g.nz );
  assert( ie >=  0 && ie < dims_g.ne );
  assert( ia >=  0 && ia < dims_g.na );
  assert( iu >=  0 && iu < NU );
  assert( octant >= 0 && octant < NOCTANT );

  return   ( (P) Quantities_affinefunction__( ia ) )
         * ( (P) Quantities_scalefactor_angle__( dims_g, ia ) )
         * ( (P) Quantities_scalefactor_space__( quan, ix_g, iy_g, iz_g ) )
         * ( (P) Quantities_scalefactor_energy__( ie, dims_g ) )
         * ( (P) Quantities_scalefactor_unknown__( iu ) );
}

/*===========================================================================*/
/*---Initial values for boundary array---*/

static inline P Quantities_init_faceyz(
  const Quantities*  quan,
  int                ix_g,
  int                iy_g,
  int                iz_g,
  int                ie,
  int                ia,
  int                iu,
  int                octant,
  const Dimensions   dims_g )
{
  assert( ( ix_g == -1        && Dir_x(octant)==Dir_up() ) ||
          ( ix_g == dims_g.nx && Dir_x(octant)==Dir_dn() ) );
  assert( iy_g >=  0 && iy_g < dims_g.ny );
  assert( iz_g >=  0 && iz_g < dims_g.nz );
  assert( ie >=  0 && ie < dims_g.ne );
  assert( ia >=  0 && ia < dims_g.na );
  assert( iu >=  0 && iu < NU );
  assert( octant >= 0 && octant < NOCTANT );

  return   ( (P) Quantities_affinefunction__( ia ) )
         * ( (P) Quantities_scalefactor_angle__( dims_g, ia ) )
         * ( (P) Quantities_scalefactor_space__( quan, ix_g, iy_g, iz_g ) )
         * ( (P) Quantities_scalefactor_energy__( ie, dims_g ) )
         * ( (P) Quantities_scalefactor_unknown__( iu ) );
}

/*===========================================================================*/
/*---Initial values for state vector---*/

static inline P Quantities_init_state(
  const Quantities*  quan,
  int                ix,
  int                iy,
  int                iz,
  int                ie,
  int                im,
  int                iu,
  const Dimensions   dims )
{
  assert( ix >= 0 && ix < dims.nx);
  assert( iy >= 0 && iy < dims.ny );
  assert( iz >= 0 && iz < dims.nz );
  assert( ie >= 0 && ie < dims.ne );
  assert( im >= 0 && im < dims.nm );
  assert( iu >= 0 && iu < NU );

  return   ( (P) Quantities_affinefunction__( im ) )
         * ( (P) Quantities_scalefactor_space__( quan,
                                                 ix+quan->ix_base,
                                                 iy+quan->iy_base, iz ) )
         * ( (P) Quantities_scalefactor_energy__( ie, dims ) )
         * ( (P) Quantities_scalefactor_unknown__( iu ) );
}

/*===========================================================================*/
/*---Perform equation solve at a gridcell---*/

static inline void Quantities_solve(
  const Quantities*  quan,
  P* __restrict__    v_local,
  P* __restrict__    facexy,
  P* __restrict__    facexz,
  P* __restrict__    faceyz,
  int                ix_b,
  int                iy_b,
  int                iz_b,
  int                ie,
  int                ix_g,
  int                iy_g,
  int                iz_g,
  int                octant,
  int                octant_ind,
  int                num_face_octants_allocated,
  const Dimensions   dims_b,
  const Dimensions   dims_g )
{
  assert( v_local );
  assert( facexy );
  assert( facexz );
  assert( faceyz );
  assert( ix_b >= 0 && ix_b < dims_b.nx );
  assert( iy_b >= 0 && iy_b < dims_b.ny );
  assert( iz_b >= 0 && iz_b < dims_b.nz );
  assert( ie   >= 0 && ie   < dims_b.ne );
  assert( ix_g >= 0 && ix_g < dims_g.nx );
  assert( iy_g >= 0 && iy_g < dims_g.ny );
  assert( iz_g >= 0 && iz_g < dims_g.nz );
  assert( octant >= 0 && octant < NOCTANT );
  assert( octant_ind >= 0 && octant_ind < num_face_octants_allocated );

  int iu = 0;
  int ia = 0;

  int dir_x = Dir_x( octant );
  int dir_y = Dir_y( octant );
  int dir_z = Dir_z( octant );

  /*---Average the face values and accumulate---*/

  /*---The state value and incoming face values are first adjusted to
       normalized values by removing the spatial scaling.
       They are then combined using a weighted average chosen in a special
       way to give just the expected result.
       Finally, spatial scaling is applied to the result which is then
       stored.
  ---*/

  for( iu=0; iu<NU; ++iu )
  for( ia=0; ia<dims_b.na; ++ia )
  {
    const P result = (
          *const_ref_v_local( v_local, dims_b, NU, ia, iu )
             / Quantities_scalefactor_space__( quan, ix_g, iy_g, iz_g )
        + *const_ref_facexy( facexy, dims_b, NU, num_face_octants_allocated,
                                            ix_b, iy_b, ie, ia, iu, octant_ind )
             * Quantities_xfluxweight__( dims_g, ia )
             / Quantities_scalefactor_space__( quan,
                                               ix_g, iy_g, iz_g-Dir_inc(dir_z) )
        + *const_ref_facexz( facexz, dims_b, NU, num_face_octants_allocated,
                                            ix_b, iz_b, ie, ia, iu, octant_ind )
             * Quantities_yfluxweight__( dims_g, ia )
             / Quantities_scalefactor_space__( quan,
                                               ix_g, iy_g-Dir_inc(dir_y), iz_g )
        + *const_ref_faceyz( faceyz, dims_b, NU, num_face_octants_allocated,
                                            iy_b, iz_b, ie, ia, iu, octant_ind )
             * Quantities_zfluxweight__( dims_g, ia )
             / Quantities_scalefactor_space__( quan,
                                               ix_g-Dir_inc(dir_x), iy_g, iz_g )
      )      * Quantities_scalefactor_space__( quan, ix_g, iy_g, iz_g );

    *ref_v_local( v_local, dims_b, NU, ia, iu ) =
    *ref_facexy( facexy, dims_b, NU, num_face_octants_allocated,
                 ix_b, iy_b, ie, ia, iu, octant_ind ) =
    *ref_facexz( facexz, dims_b, NU, num_face_octants_allocated,
                 ix_b, iz_b, ie, ia, iu, octant_ind ) =
    *ref_faceyz( faceyz, dims_b, NU, num_face_octants_allocated,
                 iy_b, iz_b, ie, ia, iu, octant_ind ) =
        result;
  } /*---for---*/

} /*---Quantities_solve---*/

/*===========================================================================*/

#endif /*---_serial_c__quantities_testing_h_---*/

/*---------------------------------------------------------------------------*/
