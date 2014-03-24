/*---------------------------------------------------------------------------*/
/*!
 * \file   sweeper_tileoctants.h
 * \author Wayne Joubert
 * \date   Wed Jan 15 16:06:28 EST 2014
 * \brief  Declarations for performing a sweep, tileoctants version.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _serial_c__sweeper_tileoctants_h_
#define _serial_c__sweeper_tileoctants_h_

#include "env.h"
#include "definitions.h"
#include "dimensions.h"
#include "arguments.h"
#include "quantities.h"

/*===========================================================================*/
/*---Struct with pointers etc. used to perform sweep---*/

typedef struct
{
  P* __restrict__  facexy;
  P* __restrict__  facexz;
  P* __restrict__  faceyz;
  P* __restrict__  v_local;

  Dimensions       dims;
} Sweeper;

/*===========================================================================*/
/*---Pseudo-constructor for Sweeper struct---*/

void Sweeper_ctor( Sweeper*          sweeper,
                   Dimensions        dims,
                   const Quantities* quan,
                   Env*              env,
                   Arguments*        args );

/*===========================================================================*/
/*---Specify whether to tile octants---*/

static Bool_t Sweeper_tile_octants()
{
  return 1;
}

/*===========================================================================*/
/*---Number of octants in an octant block---*/

static int Sweeper_noctant_per_block( const Sweeper* sweeper )
{
  return Sweeper_tile_octants() ? NOCTANT : 1;
}

/*===========================================================================*/
/*---Pseudo-destructor for Sweeper struct---*/

void Sweeper_dtor( Sweeper* sweeper );

/*===========================================================================*/
/*---Perform a sweep---*/

void Sweeper_sweep(
  Sweeper*               sweeper,
  P* __restrict__        vo,
  const P* __restrict__  vi,
  const Quantities*      quan,
  Env*                   env );

/*===========================================================================*/

#endif /*---_sweeper_tileoctants_h_---*/

/*---------------------------------------------------------------------------*/
