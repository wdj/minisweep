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
#include "arguments.h"
#include "quantities.h"
#include "step_scheduler_kba.h"

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

  Dimensions       dims;
  Dimensions       dims_b;
  Dimensions       dims_g;
  int              nblock_z;
  int              nthread_e;

  Request_t        request_send_xz;
  Request_t        request_send_yz;
  Request_t        request_recv_xz;
  Request_t        request_recv_yz;

  Step_Scheduler   step_scheduler;
} Sweeper;

/*===========================================================================*/
/*---Pseudo-constructor for Sweeper struct---*/

void Sweeper_ctor( Sweeper*          sweeper,
                   Dimensions        dims,
                   const Quantities* quan,
                   Env*              env,
                   Arguments*        args );

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
/*---Determine whether to send a face now---*/

Bool_t Sweeper_must_do_send__(
  Sweeper*           sweeper,
  int                step,
  int                axis,
  int                dir_ind,
  Env*               env );

/*===========================================================================*/
/*---Determine whether to receive a face now---*/

Bool_t Sweeper_must_do_recv__(
  Sweeper*           sweeper,
  int                step,
  int                axis,
  int                dir_ind,
  Env*               env );

/*===========================================================================*/
/*---Communicate faces---*/

void Sweeper_communicate_faces__(
  Sweeper*         sweeper,
  int              step,
  Env*             env );

/*---------------------------------------------------------------------------*/

void Sweeper_send_faces_start__(
  Sweeper*           sweeper,
  int                step,
  Env*               env );

/*---------------------------------------------------------------------------*/

void Sweeper_send_faces_end__(
  Sweeper*           sweeper,
  int                step,
  Env*               env );

/*---------------------------------------------------------------------------*/

void Sweeper_recv_faces_start__(
  Sweeper*           sweeper,
  int                step,
  Env*               env );

/*---------------------------------------------------------------------------*/

void Sweeper_recv_faces_end__(
  Sweeper*           sweeper,
  int                step,
  Env*               env );

/*===========================================================================*/
/*---Apply boundary condition: xy face---*/

static void Sweeper_set_boundary_xy(
  const Sweeper*        sweeper,
  const Quantities*     quan,
  P* const __restrict__ facexy_c,
  int                   octant );

/*===========================================================================*/
/*---Apply boundary condition: xz face---*/

static void Sweeper_set_boundary_xz(
  const Sweeper*        sweeper,
  const Quantities*     quan,
  P* const __restrict__ facexz_c,
  int                   octant,
  int                   block_z );

/*===========================================================================*/
/*---Apply boundary condition: yz face---*/

static void Sweeper_set_boundary_yz(
  const Sweeper*        sweeper,
  const Quantities*     quan,
  P* const __restrict__ faceyz_c,
  int                   octant,
  int                   block_z );

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
/*---Select which part of v_local to use for current thread---*/

static inline P* __restrict__ Sweeper_v_local_this__( Sweeper* sweeper,
                                                      int      thread_num )
{
  return sweeper->v_local + sweeper->dims_b.na * NU * thread_num;
}

/*===========================================================================*/
/*---Perform a sweep for a block---*/

void Sweeper_sweep_block(
  Sweeper*               sweeper,
  P* __restrict__        vo,
  const P* __restrict__  vi,
  const Quantities*      quan,
  Env*                   env,
  const Step_Info        step_info,
  const int              thread_num,  
  const int              num_threads, 
  const int              octant_ind,
  P* __restrict__        facexy_c,
  P* __restrict__        facexz_c,
  P* __restrict__        faceyz_c );

/*===========================================================================*/
/*---Perform a sweep---*/

void Sweeper_sweep(
  Sweeper*               sweeper,
  P* __restrict__        vo,
  const P* __restrict__  vi,
  const Quantities*      quan,
  Env*                   env );

/*===========================================================================*/

#endif /*---_sweeper_kba_h_---*/

/*---------------------------------------------------------------------------*/
