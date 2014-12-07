/*---------------------------------------------------------------------------*/
/*!
 * \file   env_declarations.h
 * \author Wayne Joubert
 * \date   Wed Jan 15 16:06:28 EST 2014
 * \brief  Declarations relevant to programming API being used.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _env_declarations_h_
#define _env_declarations_h_

#ifdef USE_MPI
#include "mpi.h"
#endif

#ifdef USE_CUDA
#include "cuda.h"
#endif

#ifdef __cplusplus
extern "C"
{
#endif

/*===========================================================================*/
/*---Types---*/

#ifdef USE_MPI
typedef MPI_Comm    Comm_t;
typedef MPI_Request Request_t;
#else
typedef int Comm_t;
typedef int Request_t;
#endif

#ifdef __CUDACC__
typedef cudaStream_t Stream_t;
#else
typedef int Stream_t;
#endif

/*===========================================================================*/
/*---Struct with environment information---*/

typedef struct
{
#ifdef USE_MPI
  int    nproc_x__;    /*---Number of procs along x axis---*/
  int    nproc_y__;    /*---Number of procs along y axis---*/
  int    tag__;        /*---Next free message tag---*/
  Comm_t active_comm__;
  Bool_t is_proc_active__;
#endif
#ifdef __CUDACC__
  Bool_t   is_using_device__;
  Stream_t stream_send_block__;
  Stream_t stream_recv_block__;
  Stream_t stream_kernel_faces__;
#endif
} Env;

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

#endif /*---_env_declarations_h_---*/

/*---------------------------------------------------------------------------*/
