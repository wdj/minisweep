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

#ifdef __cplusplus
extern "C"
{
#endif

/*===========================================================================*/
/*---Pseudo-constructor for Step_Scheduler struct---*/

void Step_Scheduler_ctor( Step_Scheduler* step_scheduler,
                          int             nblock_z,
                          int             nblock_octant,
                          Env*            env )
{
  Insist( nblock_z > 0 && "Invalid z blocking factor supplied." );
  step_scheduler->nblock_z__      = nblock_z;
  step_scheduler->nproc_x__       = Env_nproc_x( env );
  step_scheduler->nproc_y__       = Env_nproc_y( env );
  step_scheduler->nblock_octant__ = nblock_octant;
}

/*===========================================================================*/
/*---Pseudo-destructor for Step_Scheduler struct---*/

void Step_Scheduler_dtor( Step_Scheduler* step_scheduler )
{
}

/*===========================================================================*/
/*---Accessor: blocks along z axis---*/

int Step_Scheduler_nblock_z( const Step_Scheduler* step_scheduler )
{
  return step_scheduler->nblock_z__;
}

/*===========================================================================*/
/*---Number of block steps executed for a single octant in isolation---*/

int Step_Scheduler_nblock( const Step_Scheduler* step_scheduler )
{
  return step_scheduler->nblock_z__;
}

/*===========================================================================*/
/*---Number of octants per octant block---*/

int Step_Scheduler_noctant_per_block( const Step_Scheduler* step_scheduler )
{
  return NOCTANT / step_scheduler->nblock_octant__;
}

/*===========================================================================*/
/*---Number of kba parallel steps---*/

int Step_Scheduler_nstep( const Step_Scheduler* step_scheduler )
{
  int result;

  switch( step_scheduler->nblock_octant__ )
  {
    case 8:
      result = 8 * Step_Scheduler_nblock( step_scheduler )
                                       + 2 * ( step_scheduler->nproc_x__ - 1 )
                                       + 3 * ( step_scheduler->nproc_y__ - 1 );
      break;

    case 4:
      result = 4 * Step_Scheduler_nblock( step_scheduler )
                                       + 1 * ( step_scheduler->nproc_x__ - 1 )
                                       + 2 * ( step_scheduler->nproc_y__ - 1 );
      break;

    case 2:
      result = 2 * Step_Scheduler_nblock( step_scheduler )
                                       + 1 * ( step_scheduler->nproc_x__ - 1 )
                                       + 1 * ( step_scheduler->nproc_y__ - 1 );
      break;

    case 1:
      result = 1 * Step_Scheduler_nblock( step_scheduler )
                                       + 1 * ( step_scheduler->nproc_x__ - 1 )
                                       + 1 * ( step_scheduler->nproc_y__ - 1 );
      break;

    default:
      assert( Bool_false );
      break;
  }

  return result;
}

/*===========================================================================*/
/*---Get information describing a sweep step---*/

Step_Info Step_Scheduler_step_info( const Step_Scheduler* step_scheduler,  
                                    const int             step,
                                    const int             octant_in_block,
                                    const int             proc_x,
                                    const int             proc_y )
{
  assert( octant_in_block>=0 &&
          octant_in_block * step_scheduler->nblock_octant__ < NOCTANT );

  const int nblock_octant = step_scheduler->nblock_octant__;
  const int nproc_x       = step_scheduler->nproc_x__;
  const int nproc_y       = step_scheduler->nproc_y__;
  const int nblock        = Step_Scheduler_nblock( step_scheduler );
  const int nstep         = Step_Scheduler_nstep( step_scheduler );

  int octant_block = 0;
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

  Step_Info step_info;

  const int octant_selector[NOCTANT] = { 0, 4, 2, 6, 3, 7, 1, 5 };

  /*===========================================================================
    For a given step and octant_in_block, the following computes the
    octant block (i.e., octant step), from which the octant can be
    computed, and the wavefront number, starting from the relevant begin
    corner of the selected octant.
    For the nblock_octant==8 case, the 8 octants are processed in sequence,
    in the order xyz = +++, ++-, -++, -+-, --+, ---, +-+, +--.
    This order is chosen to "pack" the wavefronts to minimize
    the KBA wavefront startup latency.
    For nblock_octant=k for some smaller k, this sequence is divided into
    subsequences of length k, and each subsequence defines the schedule
    for a given octant_in_block.
    The code below is essentially a search into the first subsequence
    to determine where the requested step is located.  Locations in
    the other subsequences can be derived from this.
    NOTE: the following does not address possibility that for a single
    step, two or more octants could update the same block.
  ===========================================================================*/

  if ( Bool_true && nblock_octant >= 1 )
  {
    wave = step - ( step_base );
    octant_block = 0;
  }
  step_base += nblock;
  if ( step >= ( step_base + proc_x + proc_y ) && nblock_octant >= 2 )
  {
    wave = step - ( step_base );
    octant_block = 1;
  }
  step_base += nblock;
  if ( step >= ( step_base + proc_x + proc_y ) && nblock_octant >= 4 )
  {
    wave = step - ( step_base + (nproc_y-1) );
    octant_block = 2;
  }
  step_base += nblock + (nproc_y-1);
  if ( step >= ( step_base + (nproc_y-1-proc_y)
                           +            proc_x ) && nblock_octant >=4 )
  {
    wave = step - ( step_base );
    octant_block = 3;
  }
  step_base += nblock;
  if ( step >= ( step_base + (nproc_y-1-proc_y)
                           +            proc_x ) && nblock_octant >= 8 )
  {
    wave = step - ( step_base + (nproc_x-1) );
    octant_block = 4;
  }
  step_base += nblock + (nproc_x-1);
  if ( step >= ( step_base + (nproc_y-1-proc_y)
                           + (nproc_x-1-proc_x) ) && nblock_octant >= 8 )
  {
    wave = step - ( step_base );
    octant_block = 5;
  }
  step_base += nblock;
  if ( step >= ( step_base + (nproc_y-1-proc_y)
                           + (nproc_x-1-proc_x) ) && nblock_octant >= 8 )
  {
    wave = step - ( step_base + (nproc_y-1) );
    octant_block = 6;
  }
  step_base += nblock + (nproc_y-1);
  if ( step >= ( step_base +            proc_y
                           + (nproc_x-1-proc_x) ) && nblock_octant >= 8 )
  {
    wave = step - ( step_base );
    octant_block = 7;
  }

  octant = octant_selector[ octant_block + nblock_octant * octant_in_block ];

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

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

#endif /*---_serial_c__step_scheduler_kba_c_h_---*/

/*---------------------------------------------------------------------------*/
