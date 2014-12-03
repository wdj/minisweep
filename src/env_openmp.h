/*---------------------------------------------------------------------------*/
/*!
 * \file   env_openmp.h
 * \author Wayne Joubert
 * \date   Wed Jan 15 16:06:28 EST 2014
 * \brief  Environment settings for openmp.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _env_openmp_h_
#define _env_openmp_h_

#ifdef USE_OPENMP
#include "omp.h"
#endif

#include "types.h"
#include "env_assert.h"
#include "arguments.h"

#include "env_openmp_kernels.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*===========================================================================*/
/*---Initialize OpenMP---*/

static void Env_omp_initialize__( Env *env, int argc, char** argv )
{
#ifdef USE_OPENMP
#endif
}

/*===========================================================================*/
/*---Finalize OpenMP---*/

static void Env_omp_finalize__( Env* env )
{
#ifdef USE_OPENMP
#endif
}

/*===========================================================================*/
/*---Set values from args---*/

static void Env_omp_set_values__( Env *env, Arguments* args )
{
#ifdef USE_OPENMP
#endif
}

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

#endif /*---_env_openmp_h_---*/

/*---------------------------------------------------------------------------*/
