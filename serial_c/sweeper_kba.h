/*---------------------------------------------------------------------------*/
/*!
 * \file   sweeper_kba.h
 * \author Wayne Joubert
 * \date   Tue Jan 28 16:37:41 EST 2014
 * \brief  Declarations for performing a sweep, kba version.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _serial_c__sweeper_kba_h_
#define _serial_c__sweeper_kba_h_

#include "env.h"
#include "definitions.h"
#include "dimensions.h"
#include "quantities.h"

/*===========================================================================*/
/*---Struct with pointers etc. used to perform sweep---*/

typedef struct
{
  P* __restrict__  facexy;
  P* __restrict__  facexz;
  P* __restrict__  faceyz;
  P* __restrict__  v_local;
  Dimensions       dims_b;
  int              nblock_z;
} Sweeper;

/*===========================================================================*/
/*---Struct with info describing a sweep step---*/

typedef struct
{
  int     block_z;
  int     octant;
  Bool_t  is_active;
} Step_Info;

/*===========================================================================*/
/*---Pseudo-constructor for Sweeper struct---*/

void Sweeper_ctor( Sweeper*    sweeper,
                   Dimensions  dims,
                   Env*        env,
                   int         nblock_z );

/*===========================================================================*/
/*---Pseudo-destructor for Sweeper struct---*/

void Sweeper_dtor( Sweeper* sweeper );

/*===========================================================================*/
/*---Number of octants to store for each face---*/

static int Sweeper_num_face_octants_allocated()
{
  return 1;
}

/*===========================================================================*/
/*---Number of block steps executed for an octant in isolation---*/

int Sweeper_nblock( const Sweeper*  sweeper,
                    const Env*      env );

/*===========================================================================*/
/*---Number of kba parallel steps---*/

int Sweeper_nstep( const Sweeper*  sweeper,
                   const Env*      env );

/*===========================================================================*/
/*---Get information describing a sweep step---*/

Step_Info Sweeper_step_info__( const Sweeper*  sweeper,  
                               int             step,
                               int             proc_x,
                               int             proc_y,
                               const Env*      env );

/*===========================================================================*/
/*---Communicate faces---*/

void Sweeper_communicate_faces__(
  Sweeper*         sweeper,
  int              step,
  Dimensions       dims_b,
  Env*             env );

/*===========================================================================*/
/*---Perform a sweep---*/

void Sweeper_sweep(
  Sweeper*               sweeper,
  P* __restrict__        vo,
  const P* __restrict__  vi,
  const Quantities*      quan,
  Dimensions             dims,
  Env*                   env );

/*===========================================================================*/

#endif /*---_sweeper_kba_h_---*/

/*---------------------------------------------------------------------------*/
