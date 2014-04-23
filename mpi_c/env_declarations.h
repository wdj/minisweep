/*---------------------------------------------------------------------------*/
/*!
 * \file   env_declarations.h
 * \author Wayne Joubert
 * \date   Wed Jan 15 16:06:28 EST 2014
 * \brief  Declarations relevant to programming API being used.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _mpi_c__env_declarations_h_
#define _mpi_c__env_declarations_h_

#ifdef __cplusplus
extern "C"
{
#endif

/*===========================================================================*/
/*---Struct with environment information---*/

typedef struct
{
#ifdef USE_MPI
  int nproc_x__;    /*---Number of procs along x axis---*/
  int nproc_y__;    /*---Number of procs along y axis---*/
  int tag__;        /*---Next free message tag---*/
#endif
#ifdef __CUDACC__
  Bool_t is_using_device__;
#endif
} Env;

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

#endif /*---_mpi_c__env_declarations_h_---*/

/*---------------------------------------------------------------------------*/
