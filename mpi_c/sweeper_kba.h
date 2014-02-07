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

  P* __restrict__  facexz0;
  P* __restrict__  faceyz0;

  P* __restrict__  facexz1;
  P* __restrict__  faceyz1;

  P* __restrict__  facexz2;
  P* __restrict__  faceyz2;

  P* __restrict__  v_local;

  Dimensions       dims_b;
  int              nblock_z;

  Request_t        request_send_xz;
  Request_t        request_send_yz;
  Request_t        request_recv_xz;
  Request_t        request_recv_yz;

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
/*---Is face communication done asynchronously---*/

static int Sweeper_is_face_comm_async()
{
  return Bool_true;
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

/*---------------------------------------------------------------------------*/

void Sweeper_send_faces_start__(
  Sweeper*           sweeper,
  int                step,
  Dimensions         dims_b,
  Env*               env );

/*---------------------------------------------------------------------------*/

void Sweeper_send_faces_end__(
  Sweeper*           sweeper,
  int                step,
  Dimensions         dims_b,
  Env*               env );

/*---------------------------------------------------------------------------*/

void Sweeper_recv_faces_start__(
  Sweeper*           sweeper,
  int                step,
  Dimensions         dims_b,
  Env*               env );

/*---------------------------------------------------------------------------*/

void Sweeper_recv_faces_end__(
  Sweeper*           sweeper,
  int                step,
  Dimensions         dims_b,
  Env*               env );

/*===========================================================================*/
/*---Selectors for faces---*/

/*---The face arrays form a circular buffer of length three.  Three are
     needed because at any step there may be a send, a receive, and a
     block-sweep-compute in-flight.
---*/

static P* __restrict__ Sweeper_facexz_c__( Sweeper* sweeper, int step )
{
  assert( sweeper != NULL );
  assert( step >= 0 );
  P* __restrict__ facesxz[3] = { sweeper->facexz0,
                                 sweeper->facexz1,
                                 sweeper->facexz2 };
  return facesxz[(step+3)%3];
}

/*---------------------------------------------------------------------------*/

static P* __restrict__ Sweeper_faceyz_c__( Sweeper* sweeper, int step )
{
  assert( sweeper != NULL );
  assert( step >= 0 );
  P* __restrict__ facesyz[3] = { sweeper->faceyz0,
                                 sweeper->faceyz1,
                                 sweeper->faceyz2 };
  return facesyz[(step+3)%3];
}

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
