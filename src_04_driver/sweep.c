/*---------------------------------------------------------------------------*/
/*!
 * \file   sweep.c
 * \author Wayne Joubert
 * \date   Wed May 22 11:22:14 EDT 2013
 * \brief  Main driver for sweep miniapp.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include <stdio.h>

#include "arguments.h"
#include "env.h"
#include "definitions.h"
#include "dimensions.h"
#include "memory.h"
#include "pointer.h"
#include "quantities.h"
#include "array_operations.h"
#include "sweeper.h"

#include "run_tools.h"

/*===========================================================================*/
/*---Main---*/

int main( int argc, char** argv )
{
  /*---Declarations---*/
  Env env;

  /*---Initialize for execution---*/

  Env_initialize( &env, argc, argv );

  Arguments args;
  RunData  rundata;

  Arguments_ctor( &args, argc, argv );
  Env_set_values( &env, &args );

  /*---Perform run---*/

  if( Env_is_proc_active( &env ) )
  {
    run_case( &env, &args, &rundata );
  }

  if( Env_is_proc_master( &env ) )
  {
    printf( "Normsq result: %.8e  diff: %.3e  %s  time: %.3f  GF/s: %.3f\n",
            (double)rundata.normsq, (double)rundata.normsqdiff,
            rundata.normsqdiff==P_zero() ? "PASS" : "FAIL",
            (double)rundata.time, rundata.floprate );
  }

  /*---Deallocations---*/

  Arguments_dtor( &args );

  /*---Finalize execution---*/

  Env_finalize( &env );

} /*---main---*/

/*---------------------------------------------------------------------------*/
