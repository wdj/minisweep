/*---------------------------------------------------------------------------*/
/*!
 * \file   sweeper.h
 * \author Wayne Joubert
 * \date   Wed Jan 15 16:06:28 EST 2014
 * \brief  Definitions for performing a sweep.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _sweep_h_
#define _sweep_h_

/* #ifndef SWEEPER_SIMPLE */
/* #ifndef SWEEPER_TILEOCTANTS */
/* #ifndef SWEEPER_KBA */
/* #define SWEEPER_ACC */
/* #endif */
/* #endif */
/* #endif */

#ifdef USE_ACC
#define SWEEPER_ACC
#else
#define SWEEPER_KBA
#endif

#ifdef SWEEPER_SIMPLE
#include "sweeper_simple.h"
#endif

#ifdef SWEEPER_TILEOCTANTS
#include "sweeper_tileoctants.h"
#endif

#ifdef SWEEPER_KBA
#include "sweeper_kba.h"
#endif

#ifdef SWEEPER_ACC
#include "sweeper_acc.h"
#endif

#endif /*---_sweep_h_---*/

/*---------------------------------------------------------------------------*/
