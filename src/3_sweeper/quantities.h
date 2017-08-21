/*---------------------------------------------------------------------------*/
/*!
 * \file   quantities.h
 * \author Wayne Joubert
 * \date   Wed Jan 15 16:06:28 EST 2014
 * \brief  Declarations for physical quantities, testing case.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _quantities_h_
#define _quantities_h_

#include "types.h"
#include "dimensions.h"

#include "quantities_kernels.h"

/*===========================================================================*/
/*---Pseudo-constructor for Quantities struct---*/

void Quantities_create( Quantities*       quan,
                        const Dimensions  dims,
                        Env*              env );

/*===========================================================================*/
/*---Pseudo-destructor for Quantities struct---*/

void Quantities_destroy( Quantities* quan, Env* env );

/*===========================================================================*/
/*---Initialize Quantities a_from_m, m_from_a matrices---*/
/*---pseudo-private member function---*/

void Quantities_init_am_matrices_( Quantities*       quan,
                                   const Dimensions  dims,
                                   Env*              env );

/*===========================================================================*/
/*---Initialize Quantities subgrid decomp info---*/
/*---pseudo-private member function---*/

void Quantities_init_decomp_( Quantities*       quan,
                              const Dimensions  dims,
                              Env*              env );

/*===========================================================================*/
/*---Flops cost of solve per element---*/

double Quantities_flops_per_solve( const Dimensions dims );

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
  ASSERT( ix >= 0 && ix < dims.ncell_x);
  ASSERT( iy >= 0 && iy < dims.ncell_y );
  ASSERT( iz >= 0 && iz < dims.ncell_z );
  ASSERT( ie >= 0 && ie < dims.ne );
  ASSERT( im >= 0 && im < dims.nm );
  ASSERT( iu >= 0 && iu < NU );

  if( Quantities_bc_vacuum() )
  {
    return ((P)0);
  }
  else
  {
    return   ( (P) Quantities_affinefunction_( im ) )
           * ( (P) Quantities_scalefactor_space_( quan,
                                                   ix+quan->ix_base,
                                                   iy+quan->iy_base, iz ) )
           * ( (P) Quantities_scalefactor_energy_( ie, dims ) )
           * ( (P) Quantities_scalefactor_unknown_( iu ) );
  }
}

/*===========================================================================*/

#endif /*---_quantities_h_---*/

/*---------------------------------------------------------------------------*/
