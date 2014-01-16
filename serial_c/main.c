/*---------------------------------------------------------------------------*/
/*!
 * \file   main.c
 * \author Wayne Joubert
 * \date   Wed May 22 11:22:14 EDT 2013
 * \brief  Main driver for KBA sweep miniapp.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "environment.h"
#include "dimensions.h"
#include "memory.h"
#include "sweeper.h"

/*===========================================================================*/
/*---Initialize state vector to required input value---*/

void initialize_state( P* __restrict__ v, Dimensions dims )
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
  for( iu=0; iu<NU; ++iu )
  {
    *ref_state( v, dims, ix, iy, iz, ie, im, iu )
                            = Quantities_init_state( ix, iy, iz, ie, im, iu );

  }
}

/*===========================================================================*/
/*---Initialize state vector to zero---*/

void initialize_state_zero( P* __restrict__ v, Dimensions dims )
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
  for( iu=0; iu<NU; ++iu )
  {
    *ref_state( v, dims, ix, iy, iz, ie, im, iu ) = P_zero();
  }
}

/*===========================================================================*/
/*---Compute vector norm info for state vector---*/

void get_state_norms( P* __restrict__ vi,
                      P* __restrict__ vo,
                      Dimensions dims,
                      P* normsqp,
                      P* normsqdiffp )
{
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
  for( iu=0; iu<NU; ++iu )
  {
    normsq += *ref_state( vo, dims, ix, iy, iz, ie, im, iu ) *
              *ref_state( vo, dims, ix, iy, iz, ie, im, iu );
    normsqdiff +=
             ( *ref_state( vo, dims, ix, iy, iz, ie, im, iu ) -
               *ref_state( vi, dims, ix, iy, iz, ie, im, iu ) ) *
             ( *ref_state( vo, dims, ix, iy, iz, ie, im, iu ) -
               *ref_state( vi, dims, ix, iy, iz, ie, im, iu ) );
  }

  *normsqp     = normsq;
  *normsqdiffp = normsqdiff;
}

/*===========================================================================*/
/*---Main---*/

int main( int argc, char** argv )
{
  /*---Declarations---*/

  Dimensions dims;
  Quantities quan;
  Sweeper sweeper;

  P* vi = 0;
  P* vo = 0;

  P normsq     = P_zero();
  P normsqdiff = P_zero();

  int iteration = 0;
  int numiterations = 0;
  int tile_octants = 0;

  double t1 = 0.;
  double t2 = 0.;
  double time = 0.;
  double flops = 0.;
  double floprate = 0.;

  /*---Initialize for execution---*/

  initialize( argc, argv );

  /*---Set problem size---*/

  dims.nx       = ( argc>1 && argv[1]!="" ) ? atoi(argv[1]) : 10;
  dims.ny       = ( argc>2 && argv[2]!="" ) ? atoi(argv[2]) : 10;
  dims.nz       = ( argc>3 && argv[3]!="" ) ? atoi(argv[3]) : 10;
  dims.ne       = ( argc>4 && argv[4]!="" ) ? atoi(argv[4]) : 30;
  dims.nm       = ( argc>5 && argv[5]!="" ) ? atoi(argv[5]) : 16;
  dims.na       = ( argc>6 && argv[6]!="" ) ? atoi(argv[6]) : 33;
  tile_octants  = ( argc>7 && argv[7]!="" ) ? atoi(argv[7]) : 0;
  numiterations = ( argc>8 && argv[8]!="" ) ? atoi(argv[8]) : 1;

  /*---Initialize quantities---*/

  Quantities_ctor( &quan, dims );

  /*---Allocate arrays---*/

  vi = pmalloc( dims.nx * dims.ny * dims.nz * dims.ne * dims.nm * NU );
  vo = pmalloc( dims.nx * dims.ny * dims.nz * dims.ne * dims.nm * NU );

  /*---Initialize input state array---*/

  initialize_state( vi, dims );

  /*---Initialize output state array---*/
  /*---This is not strictly required for the output vector but might
       have a performance effect from pre-touching pages.
  ---*/

  initialize_state_zero( vo, dims );

  /*---Initialize sweeper---*/

  Sweeper_ctor( &sweeper, dims, tile_octants );

  /*---Call sweeper---*/

  t1 = get_time();

  for( iteration=0; iteration<numiterations; ++iteration )
  {
    Sweeper_sweep( &sweeper,
                   iteration%2==0 ? vo : vi,
                   iteration%2==0 ? vi : vo,
                   quan,
                   dims );
  }

  t2 = get_time();
  time = t2 - t1;

  /*---Compute flops used---*/

  flops = 1. * NOCTANT *
          1. * dims.ne *
          1. * dims.nx *
          1. * dims.ny *
          1. * dims.nz *
          1. * NU *
          ( 2    * dims.nm * dims.na +
            NDIM * dims.na +
            2    * dims.nm * dims.na ) *
          numiterations;

  floprate = time == 0. ? 0. : flops / time / 1e9;

  /*---Compute, print norm squared of result---*/

  get_state_norms( vi, vo, dims, &normsq, &normsqdiff );

  if( do_output )
  {
    printf( "Normsq result: %e  diff: %e  %s  time: %.3f  GF/s %.4f\n",
            (double)normsq, (double)normsqdiff,
            normsqdiff==P_zero() ? "PASS" : "FAIL",
            time, floprate );
  }

  /*---Deallocations---*/

  pfree( vi );
  pfree( vo );
  Sweeper_dtor( &sweeper );
  Quantities_dtor( &quan );

  /*---Finalize execution---*/

  finalize();

} /*---main---*/

/*===========================================================================*/
