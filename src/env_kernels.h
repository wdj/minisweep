/*---------------------------------------------------------------------------*/
/*!
 * \file   env_kernels.h
 * \author Wayne Joubert
 * \date   Wed Jan 15 16:06:28 EST 2014
 * \brief  env, code for device kernels.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _env_kernels_h_
#define _env_kernels_h_

/*===========================================================================*/
/*---OpenMP wrapper function definitions---*/

#include "env_openmp_kernels.h"

/*===========================================================================*/
/*---CUDA wrapper function definitions---*/

#include "env_cuda_kernels.h"

/*===========================================================================*/
/*---Intel MIC definitions---*/

#include "env_mic_kernels.h"

#endif /*---_env_kernels_h_---*/

/*---------------------------------------------------------------------------*/
