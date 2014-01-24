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

#include "definitions.h"
#include "dimensions.h"
#include "array_accessors.h"
#include "quantities.h"

/*===========================================================================*/
/*---Initialize state vector to required input value---*/

static void initialize_state( P* __restrict__ v, Dimensions dims, int nu )
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
                       = Quantities_init_state( ix, iy, iz, ie, im, iu, dims );
  }
}

/*===========================================================================*/
/*---Initialize state vector to zero---*/

static void initialize_state_zero( P* __restrict__ v, Dimensions dims, int nu )
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
    *ref_state( v, dims, nu, ix, iy, iz, ie, im, iu ) = P_zero();
  }
}

/*===========================================================================*/
/*---Compute vector norm info for state vector---*/

static void get_state_norms( P* __restrict__  vi,
                             P* __restrict__  vo,
                             Dimensions       dims,
                             int              nu,
                             P*               normsqp,
                             P*               normsqdiffp )
{
  assert( normsqp );
  assert( normsqdiffp );

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
    normsq += *ref_state( vo, dims, nu, ix, iy, iz, ie, im, iu ) *
              *ref_state( vo, dims, nu, ix, iy, iz, ie, im, iu );
    normsqdiff +=
             ( *ref_state( vo, dims, nu, ix, iy, iz, ie, im, iu ) -
               *ref_state( vi, dims, nu, ix, iy, iz, ie, im, iu ) ) *
             ( *ref_state( vo, dims, nu, ix, iy, iz, ie, im, iu ) -
               *ref_state( vi, dims, nu, ix, iy, iz, ie, im, iu ) );
  }
  assert( normsq     >= P_zero() );
  assert( normsqdiff >= P_zero() );

  *normsqp     = normsq;
  *normsqdiffp = normsqdiff;
}

/*===========================================================================*/

#endif /*---_serial_c__array_operations_h_---*/

/*---------------------------------------------------------------------------*/
