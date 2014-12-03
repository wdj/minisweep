/*---------------------------------------------------------------------------*/
/*!
 * \file   step_scheduler_kba_kernels.h
 * \author Wayne Joubert
 * \date   Tue Jan 28 16:37:41 EST 2014
 * \brief  step_scheduler_kba, code for device kernels.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _mpi_c__step_scheduler_kba_kernels_h_
#define _mpi_c__step_scheduler_kba_kernels_h_

#include "types_kernels.h"
#include "definitions_kernels.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*===========================================================================*/
/*---Struct with info describing a sweep step---*/

typedef struct
{
  int     block_z;
  int     octant;
  Bool_t  is_active;
} Step_Info;

/*===========================================================================*/
/*---8 copies of the same---*/

typedef struct
{
  Step_Info step_info[NOCTANT];
} Step_Info_Values;

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

#endif /*---_mpi_c__step_scheduler_kba_kernels_h_---*/

/*---------------------------------------------------------------------------*/
