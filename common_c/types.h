/*---------------------------------------------------------------------------*/
/*!
 * \file   types.h
 * \author Wayne Joubert
 * \date   Wed Jan 15 16:06:28 EST 2014
 * \brief  Declarations for basic types.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _common_c__types_h_
#define _common_c__types_h_

#ifdef __cplusplus
extern "C"
{
#endif

/*===========================================================================*/
/*---Basic types---*/

/*---Boolean type---*/

typedef int Bool_t;

enum{ Bool_true = 1, Bool_false = 0 };

/*---Default floating point type---*/

typedef double P;
enum{ P_IS_DOUBLE = Bool_true };

static inline P P_zero() { return (P)0; }
static inline P P_one()  { return (P)1; }

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

#endif /*---_common_c__types_h_---*/

/*---------------------------------------------------------------------------*/
