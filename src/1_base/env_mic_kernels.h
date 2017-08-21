/*---------------------------------------------------------------------------*/
/*!
 * \file   env_mic_kernels.h
 * \author Wayne Joubert
 * \date   Wed Jun 11 09:33:15 EDT 2014
 * \brief  Environment settings for Intel MIC, code for comp. kernel.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _env_mic_kernels_h_
#define _env_mic_kernels_h_

#include "types_kernels.h"
#include "env_assert_kernels.h"

/*===========================================================================*/
/*---Enums---*/

#ifdef IS_MIC
enum{ IS_USING_MIC = Bool_true };
#else
enum{ IS_USING_MIC = Bool_false };
#endif

#ifdef IS_MIC
enum{ VEC_LEN = P_IS_DOUBLE ? 8 : 16 };
#endif

/*===========================================================================*/
/*---Define dummy function/macros for non-MIC case---*/

#ifndef IS_MIC
#define __assume( a )
#define __assume_aligned( a, b )
#define __declspec( a )
#endif

/*===========================================================================*/

#endif /*---_env_mic_kernels_h_---*/

/*---------------------------------------------------------------------------*/
