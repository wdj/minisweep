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
  Dimensions       dims_block;
  int              nblock_z;
} Sweeper;

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

static int Sweeper_num_face_octants()
{
  return 1;
}

/*===========================================================================*/
/*---Number of kba parallel steps---*/

int Sweeper_nstep( Sweeper* sweeper,
                   Env env );

/*===========================================================================*/
/*---Wehther this proc active for a given sweep step---*/

Bool_t Sweeper_step_active( Sweeper* sweeper,
                            int      step,
                            int      proc_x,
                            int      proc_y,
                            Env      env );

/*===========================================================================*/
/*---Octant to compute for a given sweep step---*/

int Sweeper_octant( Sweeper* sweeper,
                     int     step,
                     int     proc_x,
                     int     proc_y,
                     Env     env );

/*===========================================================================*/
/*---Z block number to compute for a given sweep step---*/

int Sweeper_block_z( Sweeper* sweeper,
                     int      step,
                     int      proc_x,
                     int      proc_y,
                     Env      env );

/*===========================================================================*/
/*---Perform a sweep---*/

void Sweeper_sweep(
  Sweeper*         sweeper,
  P* __restrict__  vo,
  P* __restrict__  vi,
  Quantities       quan,
  Dimensions       dims,
  Env*             env );

/*===========================================================================*/

#endif /*---_sweeper_kba_h_---*/

/*---------------------------------------------------------------------------*/
