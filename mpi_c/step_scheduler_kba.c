/*---------------------------------------------------------------------------*/
/*!
 * \file   step_scheduler_kba_c.h
 * \author Wayne Joubert
 * \date   Tue Jan 28 16:37:41 EST 2014
 * \brief  Definitions for managing sweep step schedule.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _serial_c__step_scheduler_kba_c_h_
#define _serial_c__step_scheduler_kba_c_h_

#include "env.h"
#include "definitions.h"
#include "step_scheduler_kba.h"

/*===========================================================================*/
/*---Pseudo-constructor for Step_Scheduler struct---*/

void Step_Scheduler_ctor( Step_Scheduler* step_scheduler,
                          int             nblock_z,
                          Env*            env )
{
  Insist( nblock_z > 0 && "Invalid z blocking factor supplied." );
  step_scheduler->nblock_z__ = nblock_z;
  step_scheduler->nproc_x__  = env->nproc_x;
  step_scheduler->nproc_y__  = env->nproc_y;
}

/*===========================================================================*/
/*---Pseudo-destructor for Step_Scheduler struct---*/

void Step_Scheduler_dtor( Step_Scheduler* step_scheduler )
{
}

/*===========================================================================*/
/*---Number of block steps executed for a single octant in isolation---*/

int Step_Scheduler_nblock( const Step_Scheduler* step_scheduler )
{
  return step_scheduler->nblock_z__;
}

/*===========================================================================*/
/*---Number of kba parallel steps---*/

int Step_Scheduler_nstep( const Step_Scheduler* step_scheduler )
{
  return NOCTANT * Step_Scheduler_nblock( step_scheduler )
                                       + 3 * ( step_scheduler->nproc_x__ - 1 )
                                       + 2 * ( step_scheduler->nproc_y__ - 1 );
}

/*===========================================================================*/
/*---Get information describing a sweep step---*/

Step_Info Step_Scheduler_step_info( const Step_Scheduler* step_scheduler,  
                                    const int             step,
                                    const int             proc_x,
                                    const int             proc_y )
{
/*
  assert( step >= 0 && step < Step_Scheduler_nstep( step_scheduler ) );
  assert( proc_x >= 0 && proc_x < step_scheduler->nproc_x__ );
  assert( proc_y >= 0 && proc_y < step_scheduler->nproc_y__ );
*/
  const int nproc_x = step_scheduler->nproc_x__;
  const int nproc_y = step_scheduler->nproc_y__;
  const int nblock  = Step_Scheduler_nblock( step_scheduler );
  const int nstep   = Step_Scheduler_nstep( step_scheduler );

  const int octants_visited[NOCTANT] = { 0, 4, 1, 5, 3, 7, 2, 6 };

  Step_Info step_info;

  int octant_index = 0;
  int wave         = 0;
  int step_base    = 0;
  int block        = 0;
  int octant       = 0;
  int dir_x        = 0;
  int dir_y        = 0;
  int dir_z        = 0;
  int start_x      = 0;
  int start_y      = 0;
  int start_z      = 0;

  /*---First compute the octant number, in the order they are visited,
       and the wavefront number for that octant, starting from the
       beginning corner.
       Check every octant/wavefront in sequence to determine which
       one might be active for the proc in question.
  ---*/

  if ( Bool_true )
  {
    wave = step - ( step_base );
    octant_index = 0;
  }
  step_base += nblock;
  if ( step >= ( step_base + proc_x + proc_y ) )
  {
    wave = step - ( step_base );
    octant_index = 1;
  }
  step_base += nblock;
  if ( step >= ( step_base + proc_x + proc_y ) )
  {
    wave = step - ( step_base + (nproc_x-1) );
    octant_index = 2;
  }
  step_base += nblock + (nproc_x-1);
  if ( step >= ( step_base + (nproc_x-1-proc_x) + proc_y ) )
  {
    wave = step - ( step_base );
    octant_index = 3;
  }
  step_base += nblock;
  if ( step >= ( step_base + (nproc_x-1-proc_x) + proc_y ) )
  {
    wave = step - ( step_base + (nproc_y-1) );
    octant_index = 4;
  }
  step_base += nblock + (nproc_y-1);
  if ( step >= ( step_base + (nproc_x-1-proc_x)
                           + (nproc_y-1-proc_y) ) )
  {
    wave = step - ( step_base );
    octant_index = 5;
  }
  step_base += nblock;
  if ( step >= ( step_base + (nproc_x-1-proc_x)
                           + (nproc_y-1-proc_y) ) )
  {
    wave = step - ( step_base + (nproc_x-1) );
    octant_index = 6;
  }
  step_base += nblock + (nproc_x-1);
  if ( step >= ( step_base + proc_x + (nproc_y-1-proc_y) ) )
  {
    wave = step - ( step_base );
    octant_index = 7;
  }

  octant = octants_visited[octant_index];

  /*---Next convert the wavefront number to a block number based on
       location in the domain.  Use the equation that defines the plane.
  ---*/

  dir_x  = Dir_x( octant );
  dir_y  = Dir_y( octant );
  dir_z  = Dir_z( octant );

  /*---Get coordinates of the starting corner block of the wavefront---*/
  start_x = dir_x==Dir_up() ? 0 : ( nproc_x - 1 );
  start_y = dir_y==Dir_up() ? 0 : ( nproc_y - 1 );
  start_z = dir_z==Dir_up() ? 0 : ( nblock  - 1 );

  /*---Get coordinate of block on this processor to be processed---*/
  block = ( wave - ( start_x + proc_x * dir_x)
                 - ( start_y + proc_y * dir_y)
                 - ( start_z ) ) / dir_z;

  /*---Now determine whether the block calculation is active based on whether
       the block in question falls within the physical domain.
  ---*/

  step_info.is_active = block  >= 0 && block  < nblock &&
                        step   >= 0 && step   < nstep &&
                        proc_x >= 0 && proc_x < nproc_x &&
                        proc_y >= 0 && proc_y < nproc_y;

  /*---Set remaining values---*/

  step_info.block_z = step_info.is_active ? block  : -1;
  step_info.octant  = step_info.is_active ? octant : -1;

  return step_info;
}

/*===========================================================================*/

#endif /*---_serial_c__step_scheduler_kba_c_h_---*/

/*---------------------------------------------------------------------------*/
