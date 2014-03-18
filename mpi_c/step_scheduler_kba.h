/*---------------------------------------------------------------------------*/
/*!
 * \file   step_scheduler_kba.h
 * \author Wayne Joubert
 * \date   Tue Jan 28 16:37:41 EST 2014
 * \brief  Declarations for managing sweep step schedule.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _serial_c__step_scheduler_kba_h_
#define _serial_c__step_scheduler_kba_h_

#include "env.h"

/*===========================================================================*/
/*---Struct with info describing a sweep step---*/

typedef struct
{
  int     block_z;
  int     octant;
  Bool_t  is_active;
} Step_Info;

/*===========================================================================*/
/*---Struct with info to define the sweep step schedule---*/

typedef struct
{
  int nblock_z__;
  int nproc_x__;
  int nproc_y__;
  int nblock_octant__;
} Step_Scheduler;

/*===========================================================================*/
/*---Pseudo-constructor for Step_Scheduler struct---*/

void Step_Scheduler_ctor( Step_Scheduler* step_scheduler,
                          int             nblock_z,
                          int             nblock_octant,
                          Env*            env );

/*===========================================================================*/
/*---Pseudo-destructor for Step_Scheduler struct---*/

void Step_Scheduler_dtor( Step_Scheduler* step_scheduler );

/*===========================================================================*/
/*---Accessor---*/

int Step_Scheduler_nblock_z( const Step_Scheduler* step_scheduler );

/*===========================================================================*/
/*---Number of block steps executed for a single octant in isolation---*/

int Step_Scheduler_nblock( const Step_Scheduler* step_scheduler );

/*===========================================================================*/
/*---Number of kba parallel steps---*/

int Step_Scheduler_nstep( const Step_Scheduler* step_scheduler );

/*===========================================================================*/
/*---Get information describing a sweep step---*/

Step_Info Step_Scheduler_step_info( const Step_Scheduler* step_scheduler,  
                                    const int             step,
                                    const int             octant_per_block,
                                    const int             proc_x,
                                    const int             proc_y );

/*===========================================================================*/

#endif /*---_step_scheduler_kba_h_---*/

/*---------------------------------------------------------------------------*/
