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

#include "arguments.h"
#include "env.h"
#include "definitions.h"
#include "dimensions.h"
#include "memory.h"
#include "quantities.h"
#include "array_operations.h"
#include "sweeper.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*===========================================================================*/
/*---Main---*/

int main( int argc, char** argv )
{
  /*---Declarations---*/

  Dimensions  dims_g;
  Dimensions  dims;         /*---dims for the part on this MPI proc---*/
  Arguments   args;
  Quantities  quan;
  Sweeper     sweeper;
  Env         env;

  P* vi = NULL;
  P* vo = NULL;

  P normsq     = P_zero();
  P normsqdiff = P_zero();

  int iteration   = 0;
  int niterations = 0;

  Timer_t t1      = 0.;
  Timer_t t2      = 0.;
  Timer_t time    = 0.;
  double flops    = 0.;
  double floprate = 0.;

  /*---Initialize for execution---*/

  Env_initialize( &env, argc, argv );

  Arguments_ctor( &args, argc, argv );

  Env_set_values( &env, &args );

  /*---Set problem specifications---*/

  dims_g.nx   = Arguments_consume_int_or_default( &args, "--nx",  5 );
  dims_g.ny   = Arguments_consume_int_or_default( &args, "--ny",  5 );
  dims_g.nz   = Arguments_consume_int_or_default( &args, "--nz",  5 );
  dims_g.ne   = Arguments_consume_int_or_default( &args, "--ne", 30 );
  dims_g.nm   = Arguments_consume_int_or_default( &args, "--nm", 16 );
  dims_g.na   = Arguments_consume_int_or_default( &args, "--na", 33 );
  niterations = Arguments_consume_int_or_default( &args, "--niterations",
                                                                    1 );

  Insist( dims_g.nx > 0 && "Invalid nx supplied." );
  Insist( dims_g.ny > 0 && "Invalid ny supplied." );
  Insist( dims_g.nz > 0 && "Invalid nz supplied." );
  Insist( dims_g.ne > 0 && "Invalid ne supplied." );
  Insist( dims_g.nm > 0 && "Invalid nm supplied." );
  Insist( dims_g.na > 0 && "Invalid na supplied." );
  Insist( niterations >= 0 && "Invalid iteration count supplied." );

  /*---Initialize (local) dimensions---*/

  dims = dims_g;

  dims.nx =
      ( ( Env_proc_x_this( &env ) + 1 ) * dims_g.nx ) / Env_nproc_x( &env )
    - ( ( Env_proc_x_this( &env )     ) * dims_g.nx ) / Env_nproc_x( &env );

  dims.ny =
      ( ( Env_proc_y_this( &env ) + 1 ) * dims_g.ny ) / Env_nproc_y( &env )
    - ( ( Env_proc_y_this( &env )     ) * dims_g.ny ) / Env_nproc_y( &env );

  /*---Initialize quantities---*/

  Quantities_ctor( &quan, dims, &env );

  /*---Allocate arrays---*/

  vi = malloc_P( Dimensions_size_state( dims, NU ) );
  vo = malloc_P( Dimensions_size_state( dims, NU ) );

  /*---Initialize input state array---*/

  initialize_state( vi, dims, NU, &quan );

  /*---Initialize output state array---*/
  /*---This is not strictly required for the output vector but might
       have a performance effect from pre-touching pages.
  ---*/

  initialize_state_zero( vo, dims, NU );

  /*---Initialize sweeper---*/

  Sweeper_ctor( &sweeper, dims, &quan, &env, &args );

  /*---Check that all command line args used---*/

  Insist( Arguments_are_all_consumed( &args ) && "Invalid argument detected." );

  /*---Call sweeper---*/

  t1 = Env_get_synced_time();

  for( iteration=0; iteration<niterations; ++iteration )
  {
    Sweeper_sweep( &sweeper,
                   iteration%2==0 ? vo : vi,
                   iteration%2==0 ? vi : vo,
                   &quan,
                   &env );
  }

  t2 = Env_get_synced_time();
  time = t2 - t1;

  /*---Compute flops used---*/

  flops = Env_sum_d( niterations *
            ( Dimensions_size_state( dims, NU ) * NOCTANT * 2. * dims.na
            + Dimensions_size_state_angles( dims, NU )
                                           * Quantities_flops_per_solve( dims )
            + Dimensions_size_state( dims, NU ) * NOCTANT * 2. * dims.na ) );

  floprate = time <= (Timer_t)0. ? 0. : flops / time / 1e9;

  /*---Compute, print norm squared of result---*/

  get_state_norms( vi, vo, dims, NU, &normsq, &normsqdiff );

  if( Env_do_output( &env ) )
  {
    printf( "Normsq result: %.8e  diff: %.3e  %s  time: %.3f  GF/s: %.3f\n",
            (double)normsq, (double)normsqdiff,
            normsqdiff==P_zero() ? "PASS" : "FAIL",
            (double)time, floprate );
  }

  /*---Deallocations---*/

  free_P( vi );
  free_P( vo );
  Sweeper_dtor( &sweeper );
  Quantities_dtor( &quan );
  Arguments_dtor( &args );

  /*---Finalize execution---*/

  Env_finalize( &env );

} /*---main---*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*---------------------------------------------------------------------------*/
