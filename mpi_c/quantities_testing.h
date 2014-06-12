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

#include "function_attributes.h"
#include "dimensions.h"
#include "array_accessors.h"
#include "pointer.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*===========================================================================*/
/*---Type of boundary conditions---*/

TARGET_HD static inline Bool_t Quantities_bc_vacuum()
{
  return Bool_false;
}

/*===========================================================================*/
/*---Struct to hold pointers to arrays associated with physical quantities---*/

typedef struct
{
  Pointer  a_from_m;
  Pointer  m_from_a;
  int*     ix_base_vals;
  int*     iy_base_vals;
  int      ix_base;
  int      iy_base;
  int      nx_g;
  int      ny_g;
  int      nz_g;
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
                                    const Dimensions  dims,
                                    Env*              env );

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

TARGET_HD static inline int Quantities_scalefactor_energy__( int ie,
                                                             Dimensions dims )
{
  /*---Random power-of-two multiplier for each energy group,
       to help catch errors regarding indexing of energy groups.
  ---*/
  Assert( ie >= 0 && ie < dims.ne );

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

TARGET_HD static inline int Quantities_scalefactor_unknown__( int iu )
{
  /*---Random power-of-two multiplier for each cell unknown,
       to help catch errors regarding indexing of cell unknowns.
  ---*/
  Assert( iu >= 0 && iu < NU );

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

TARGET_HD static inline int Quantities_scalefactor_space__(
                                                  const Quantities* quan,
                                                  int ix_g,
                                                  int iy_g,
                                                  int iz_g )
{
  /*---Create power of 2 based on hash of the spatial location.
  ---*/
  Assert( ix_g >= -1 && ix_g <= quan->nx_g );
  Assert( iy_g >= -1 && iy_g <= quan->ny_g );
  Assert( iz_g >= -1 && iz_g <= quan->nz_g );

  int result = 0;

#ifndef RELAXED_TESTING
  const int im = 134456;
  const int ia = 8121;
  const int ic = 28411;

  result = ( (result+(ix_g+2))*ia + ic ) % im;
  result = ( (result+(iy_g+2))*ia + ic ) % im;
  result = ( (result+(iz_g+2))*ia + ic ) % im;
  result = ( (result+(ix_g+3*iy_g+7*iz_g+2))*ia + ic ) % im;
  result = ix_g+3*iy_g+7*iz_g+2;
  result = result & ( (1<<2) - 1 );
#endif
  result = 1 << result;

  return result;
}

/*===========================================================================*/
/*---Scale factor for angles---*/
/*---pseudo-private member function---*/

TARGET_HD static inline int Quantities_scalefactor_angle__( Dimensions dims,
                                                            int ia )
{
  /*---Create a "random" power of 2. Limit the size by taking only
       the low order bits of ia
  ---*/
  Assert( ia >= 0 && ia < dims.na );

  return 1 << ( ia & ( (1<<3) - 1 ) );
}

/*===========================================================================*/
/*---Accessor for weights---*/
/*---pseudo-private member function---*/

TARGET_HD static inline P Quantities_xfluxweight__( Dimensions dims,
                                                    int ia )
{
  Assert( ia >= 0 && ia < dims.na );

  return (P) ( 1 / (P) 2 );
}

/*===========================================================================*/
/*---Accessor for weights---*/
/*---pseudo-private member function---*/

TARGET_HD static inline P Quantities_yfluxweight__( Dimensions dims,
                                                    int ia )
{
  Assert( ia >= 0 && ia < dims.na );

  return (P) ( 1 / (P) 4 );
}

/*===========================================================================*/
/*---Affine function used as basis for vector of linear values.
---pseudo-private member function---*/

TARGET_HD static inline P Quantities_affinefunction__( int i )
{
  /*---NOTE: to insure that the mappings from moments to angles to moments
       restores the original vector, this must be just like this ---*/
  return (P) ( 1 + i );
}

/*===========================================================================*/
/*---Accessor for weights---*/
/*---pseudo-private member function---*/

TARGET_HD static inline P Quantities_zfluxweight__( Dimensions dims,
                                                    int ia )
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
  Assert( ia >= 0 && ia < dims.na );

  return (P) ( 1 / (P) 4 - 1 / (P) Quantities_scalefactor_angle__( dims, ia ) );
}

/*===========================================================================*/
/*---Scale factor for octants---*/
/*---pseudo-private member function---*/

TARGET_HD static inline int Quantities_scalefactor_octant__( int octant )
{
  Assert( octant>=0 && octant<NOCTANT );

#ifndef RELAXED_TESTING
  const int result = 1 + octant;
#else
  const int result = 1;
#endif

  return result;
}

/*===========================================================================*/
/*---Initial values for boundary array---*/

TARGET_HD static inline P Quantities_init_facexy(
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
  Assert( ix_g >=  0 && ix_g <  dims_g.nx );
  Assert( iy_g >=  0 && iy_g <  dims_g.ny );
  Assert( ( iz_g == -1        && Dir_z(octant)==DIR_UP ) ||
          ( iz_g == dims_g.nz && Dir_z(octant)==DIR_DN ) );
  Assert( ie >=  0 && ie < dims_g.ne );
  Assert( ia >=  0 && ia < dims_g.na );
  Assert( iu >=  0 && iu < NU );
  Assert( octant >= 0 && octant < NOCTANT );

  /*---NOTE: this is constructed to be affine in ia (except for scale factor)
       and independent of ix, iy, iz, to facilitate calculating the
       result analytically---*/

  if( Quantities_bc_vacuum() )
  {
    return ((P)0);
  }
  else
  {
    return   ( (P) Quantities_affinefunction__( ia ) )
           * ( (P) Quantities_scalefactor_angle__( dims_g, ia ) )
           * ( (P) Quantities_scalefactor_space__( quan, ix_g, iy_g, iz_g ) )
           * ( (P) Quantities_scalefactor_energy__( ie, dims_g ) )
           * ( (P) Quantities_scalefactor_unknown__( iu ) )
           * ( (P) Quantities_scalefactor_octant__( octant ) );
  }
}

/*===========================================================================*/
/*---Initial values for boundary array---*/

TARGET_HD static inline P Quantities_init_facexz(
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
  Assert( ix_g >=  0 && ix_g < dims_g.nx );
  Assert( ( iy_g == -1        && Dir_y(octant)==DIR_UP ) ||
          ( iy_g == dims_g.ny && Dir_y(octant)==DIR_DN ) );
  Assert( iz_g >=  0 && iz_g < dims_g.nz );
  Assert( ie >=  0 && ie < dims_g.ne );
  Assert( ia >=  0 && ia < dims_g.na );
  Assert( iu >=  0 && iu < NU );
  Assert( octant >= 0 && octant < NOCTANT );

  if( Quantities_bc_vacuum() )
  {
    return ((P)0);
  }
  else
  {
    return   ( (P) Quantities_affinefunction__( ia ) )
           * ( (P) Quantities_scalefactor_angle__( dims_g, ia ) )
           * ( (P) Quantities_scalefactor_space__( quan, ix_g, iy_g, iz_g ) )
           * ( (P) Quantities_scalefactor_energy__( ie, dims_g ) )
           * ( (P) Quantities_scalefactor_unknown__( iu ) )
           * ( (P) Quantities_scalefactor_octant__( octant ) );
  }
}

/*===========================================================================*/
/*---Initial values for boundary array---*/

TARGET_HD static inline P Quantities_init_faceyz(
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
  Assert( ( ix_g == -1        && Dir_x(octant)==DIR_UP ) ||
          ( ix_g == dims_g.nx && Dir_x(octant)==DIR_DN ) );
  Assert( iy_g >=  0 && iy_g < dims_g.ny );
  Assert( iz_g >=  0 && iz_g < dims_g.nz );
  Assert( ie >=  0 && ie < dims_g.ne );
  Assert( ia >=  0 && ia < dims_g.na );
  Assert( iu >=  0 && iu < NU );
  Assert( octant >= 0 && octant < NOCTANT );

  if( Quantities_bc_vacuum() )
  {
    return ((P)0);
  }
  else
  {
    return   ( (P) Quantities_affinefunction__( ia ) )
           * ( (P) Quantities_scalefactor_angle__( dims_g, ia ) )
           * ( (P) Quantities_scalefactor_space__( quan, ix_g, iy_g, iz_g ) )
           * ( (P) Quantities_scalefactor_energy__( ie, dims_g ) )
           * ( (P) Quantities_scalefactor_unknown__( iu ) )
           * ( (P) Quantities_scalefactor_octant__( octant ) );
  }
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
  Assert( ix >= 0 && ix < dims.nx);
  Assert( iy >= 0 && iy < dims.ny );
  Assert( iz >= 0 && iz < dims.nz );
  Assert( ie >= 0 && ie < dims.ne );
  Assert( im >= 0 && im < dims.nm );
  Assert( iu >= 0 && iu < NU );

  if( Quantities_bc_vacuum() )
  {
    return ((P)0);
  }
  else
  {
    return   ( (P) Quantities_affinefunction__( im ) )
           * ( (P) Quantities_scalefactor_space__( quan,
                                                   ix+quan->ix_base,
                                                   iy+quan->iy_base, iz ) )
           * ( (P) Quantities_scalefactor_energy__( ie, dims ) )
           * ( (P) Quantities_scalefactor_unknown__( iu ) );
  }
}

/*===========================================================================*/
/*---Perform equation solve at a gridcell---*/

TARGET_HD static inline void Quantities_solve(
  const Quantities*  quan,
  P* __restrict__    vslocal,
  int                ia,
  int                iaind,
  int                iamax,
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
  int                octant_in_block,
  int                noctant_per_block,
  const Dimensions   dims_b,
  const Dimensions   dims_g,
  Bool_t             is_cell_active )
{
  Assert( vslocal );
  Assert( ia >= 0 && ia < dims_b.na );
  Assert( iaind >= 0 && iaind < iamax );
  Assert( iamax >= 0 );
  Assert( facexy );
  Assert( facexz );
  Assert( faceyz );
  Assert( ix_b >= 0 && ix_b < dims_b.nx );
  Assert( iy_b >= 0 && iy_b < dims_b.ny );
  Assert( iz_b >= 0 && iz_b < dims_b.nz );
  Assert( ie   >= 0 && ie   < dims_b.ne );
  Assert( ix_g >= 0 && ix_g < dims_g.nx );
  Assert( iy_g >= 0 && iy_g < dims_g.ny );
  Assert( iz_g >= 0 && iz_g < dims_g.nz );
  Assert( octant >= 0 && octant < NOCTANT );
  Assert( octant_in_block >= 0 && octant_in_block < noctant_per_block );

  const int dir_x = Dir_x( octant );
  const int dir_y = Dir_y( octant );
  const int dir_z = Dir_z( octant );

  int iu = 0;

  /*---Average the face values and accumulate---*/

  /*---The state value and incoming face values are first adjusted to
       normalized values by removing the spatial scaling.
       They are then combined using a weighted average chosen in a special
       way to give just the expected result.
       Finally, spatial scaling is applied to the result which is then
       stored.
  ---*/

#pragma unroll
  for( iu=0; iu<NU; ++iu )
  {
    if( ia < dims_b.na && is_cell_active )
    {
    const P result = (
          *const_ref_vslocal( vslocal, dims_b, NU, iamax, iaind, iu )
             / Quantities_scalefactor_space__( quan, ix_g, iy_g, iz_g )
        + *const_ref_facexy( facexy, dims_b, NU, noctant_per_block,
                                     ix_b, iy_b, ie, ia, iu, octant_in_block )
           * ( Quantities_xfluxweight__( dims_g, ia ) * ((P)1)
             / Quantities_scalefactor_octant__( octant ) )
             / Quantities_scalefactor_space__( quan,
                                               ix_g, iy_g, iz_g-Dir_inc(dir_z) )
        + *const_ref_facexz( facexz, dims_b, NU, noctant_per_block,
                                     ix_b, iz_b, ie, ia, iu, octant_in_block )
           * ( Quantities_yfluxweight__( dims_g, ia ) * ((P)1)
             / Quantities_scalefactor_octant__( octant ) )
             / Quantities_scalefactor_space__( quan,
                                               ix_g, iy_g-Dir_inc(dir_y), iz_g )
        + *const_ref_faceyz( faceyz, dims_b, NU, noctant_per_block,
                                     iy_b, iz_b, ie, ia, iu, octant_in_block )
           * ( Quantities_zfluxweight__( dims_g, ia ) * ((P)1)
             / Quantities_scalefactor_octant__( octant ) )
             / Quantities_scalefactor_space__( quan,
                                               ix_g-Dir_inc(dir_x), iy_g, iz_g )
      )      * Quantities_scalefactor_space__( quan, ix_g, iy_g, iz_g );

    *ref_vslocal( vslocal, dims_b, NU, iamax, iaind, iu ) = result;
    *ref_facexy( facexy, dims_b, NU, noctant_per_block,
                 ix_b, iy_b, ie, ia, iu, octant_in_block ) = result *
                 Quantities_scalefactor_octant__( octant );
    *ref_facexz( facexz, dims_b, NU, noctant_per_block,
                 ix_b, iz_b, ie, ia, iu, octant_in_block ) = result *
                 Quantities_scalefactor_octant__( octant );
    *ref_faceyz( faceyz, dims_b, NU, noctant_per_block,
                 iy_b, iz_b, ie, ia, iu, octant_in_block ) = result *
                 Quantities_scalefactor_octant__( octant );
    }
  } /*---for---*/

} /*---Quantities_solve---*/

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

#endif /*---_serial_c__quantities_testing_h_---*/

/*---------------------------------------------------------------------------*/
