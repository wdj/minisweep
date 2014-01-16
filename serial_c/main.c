/*---------------------------------------------------------------------------*/
/*!
 * \file   main.c
 * \author Wayne Joubert
 * \date   Wed May 22 11:22:14 EDT 2013
 * \brief  Main driver for KBA sweep miniapp.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include <sys/time.h>

#include <stdio.h>
#include <stdlib.h>

#ifdef USE_MPI
#include "mpi.h"
#endif

#include "sweeper.h"
#include "memory.h"

/*---------------------------------------------------------------------------*/
/*---Timer utility---*/

double get_time()
{
    struct timeval tp;
    struct timezone tzp;
    int i = gettimeofday( &tp, &tzp );
    double result = ( (double) tp.tv_sec +
                      (double) tp.tv_usec * 1.e-6 );
    return result;
}

/*---------------------------------------------------------------------------*/
/*---Main---*/

int main( int argc, char** argv )
{
  /*---Declarations---*/

  Dimensions dims;
  Quantities quan;
  Sweeper sweeper;

  P* vi = 0;
  P* vo = 0;
  P* a_from_m = 0;
  P* m_from_a = 0;

  const P zero = (P)0.;
  const P one  = (P)1.;

  P normsq = zero;
  P normsqdiff = zero;

  int ix = 0;
  int iy = 0;
  int iz = 0;
  int ie = 0;
  int im = 0;
  int ia = 0;
  int iu = 0;

  int iteration = 0;
  int numiterations = 0;

  double t1 = 0.;
  double t2 = 0.;
  double time = 0.;
  double flops = 0.;
  double floprate = 0.;

  /*---Initialize MPI if needed for Cray system backend node---*/

#ifdef USE_MPI
  MPI_Init( &argc, &argv );
#endif

  /*---Set problem size---*/

  dims.nx           = ( argc>1 && argv[1]!="" ) ? atoi(argv[1]) : 10;
  dims.ny           = ( argc>2 && argv[2]!="" ) ? atoi(argv[2]) : 10;
  dims.nz           = ( argc>3 && argv[3]!="" ) ? atoi(argv[3]) : 10;
  dims.ne           = ( argc>4 && argv[4]!="" ) ? atoi(argv[4]) : 30;
  dims.nm           = ( argc>5 && argv[5]!="" ) ? atoi(argv[5]) : 16;
  dims.na           = ( argc>6 && argv[6]!="" ) ? atoi(argv[6]) : 64;
  dims.tile_octants = ( argc>7 && argv[7]!="" ) ? atoi(argv[7]) : 0;
  numiterations     = ( argc>8 && argv[8]!="" ) ? atoi(argv[8]) : 1;

  /*---Initialize quantities---*/

  Quantities_ctor( &quan, dims );

  /*---Allocate arrays---*/

  vi = pmalloc( dims.nx * dims.ny * dims.nz * dims.ne * dims.nm * NU );
  vo = pmalloc( dims.nx * dims.ny * dims.nz * dims.ne * dims.nm * NU );

  /*---Initialize state array---*/

  for( iz=0; iz<dims.nz; ++iz )
  for( iy=0; iy<dims.ny; ++iy )
  for( ix=0; ix<dims.nx; ++ix )
  for( ie=0; ie<dims.ne; ++ie )
  for( im=0; im<dims.nm; ++im )
  for( iu=0; iu<NU; ++iu )
  {
    *ref_state( vi, dims, ix, iy, iz, ie, im, iu )
                             = Quantities_init_state( ix, iy, iz, ie, im, iu );
  }

  /*---Initialize state array---*/
  /*---This is not strictly required for the output vector but might
       have a performance effect from pre-touching pages.
  ---*/

  for( iz=0; iz<dims.nz; ++iz )
  for( iy=0; iy<dims.ny; ++iy )
  for( ix=0; ix<dims.nx; ++ix )
  for( ie=0; ie<dims.ne; ++ie )
  for( im=0; im<dims.nm; ++im )
  for( iu=0; iu<NU; ++iu )
  {
    *ref_state( vo, dims, ix, iy, iz, ie, im, iu ) = zero;
  }

  /*---Initialize sweeper---*/

  Sweeper_ctor( &sweeper, dims );

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

  /*---Terminate sweeper---*/

  Sweeper_dtor( &sweeper );

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

  normsq = zero;
  normsqdiff = zero;

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

  printf( "Normsq result: %e  diff: %e  %s  time: %.3f  GF/s %.4f\n",
          (double)normsq, (double)normsqdiff,
          normsqdiff==zero ? "PASS" : "FAIL",
          time, floprate );

  /*---Deallocate arrays---*/

  pfree( vi );
  pfree( vo );
  Quantities_dtor( &quan );

  /*---Finalize MPI if needed for Cray system backend node---*/

#ifdef USE_MPI
  MPI_Finalize();
#endif

} /*---main---*/

/*---------------------------------------------------------------------------*/
