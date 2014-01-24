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

#include "environment.h"
#include "definitions.h"
#include "dimensions.h"
#include "memory.h"
#include "quantities.h"
#include "array_operations.h"
#include "sweeper.h"

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

  double t1 = 0.;
  double t2 = 0.;
  double time = 0.;
  double flops = 0.;
  double floprate = 0.;

  /*---Initialize for execution---*/

  Env_initialize( argc, argv );

  /*---Set problem size---*/

  dims.nx       = ( argc>1 && argv[1]!="" ) ? atoi(argv[1]) : 5;
  dims.ny       = ( argc>2 && argv[2]!="" ) ? atoi(argv[2]) : 5;
  dims.nz       = ( argc>3 && argv[3]!="" ) ? atoi(argv[3]) : 5;
  dims.ne       = ( argc>4 && argv[4]!="" ) ? atoi(argv[4]) : 30;
  dims.nm       = ( argc>5 && argv[5]!="" ) ? atoi(argv[5]) : 16;
  dims.na       = ( argc>6 && argv[6]!="" ) ? atoi(argv[6]) : 33;
  numiterations = ( argc>7 && argv[7]!="" ) ? atoi(argv[7]) : 1;

  Insist( dims.nx > 0, "Invalid nx supplied." );
  Insist( dims.ny > 0, "Invalid ny supplied." );
  Insist( dims.nz > 0, "Invalid nz supplied." );
  Insist( dims.ne > 0, "Invalid ne supplied." );
  Insist( dims.nm > 0, "Invalid nm supplied." );
  Insist( dims.na > 0, "Invalid na supplied." );
  Insist( numiterations >= 0, "Invalid iteration count supplied." );

  /*---Initialize quantities---*/

  Quantities_ctor( &quan, dims );

  /*---Allocate arrays---*/

  vi = pmalloc( Dimensions_size_state( dims, NU ) );
  vo = pmalloc( Dimensions_size_state( dims, NU ) );

  /*---Initialize input state array---*/

  initialize_state( vi, dims, NU );

  /*---Initialize output state array---*/
  /*---This is not strictly required for the output vector but might
       have a performance effect from pre-touching pages.
  ---*/

  initialize_state_zero( vo, dims, NU );

  /*---Initialize sweeper---*/

  Sweeper_ctor( &sweeper, dims );

  /*---Call sweeper---*/

  t1 = Env_get_time();

  for( iteration=0; iteration<numiterations; ++iteration )
  {
    Sweeper_sweep( &sweeper,
                   iteration%2==0 ? vo : vi,
                   iteration%2==0 ? vi : vo,
                   quan,
                   dims );
  }

  t2 = Env_get_time();
  time = t2 - t1;

  /*---Compute flops used---*/

  flops = ( Dimensions_size_state( dims, NU ) * NOCTANT * 2. * dims.na
          + Dimensions_size_state_angles( dims, NU )
                                           * Quantities_flops_per_solve( dims )
          + Dimensions_size_state( dims, NU ) * NOCTANT * 2. * dims.na )
        * numiterations;

  floprate = time <= 0. ? 0. : flops / time / 1e9;

  /*---Compute, print norm squared of result---*/

  get_state_norms( vi, vo, dims, NU, &normsq, &normsqdiff );

  if( Env_do_output() )
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

  Env_finalize();

} /*---main---*/

/*---------------------------------------------------------------------------*/
