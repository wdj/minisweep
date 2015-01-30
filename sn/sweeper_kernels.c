/*---------------------------------------------------------------------------*/
/*!
 * \file   sweeper_kernels.c
 * \author Wayne Joubert
 * \date   Wed Jan 15 16:06:28 EST 2014
 * \brief  Definitions for performing a sweep. Code for device kernels.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include "sweeper.h"

#ifdef SWEEPER_SIMPLE
#endif

#ifdef SWEEPER_TILEOCTANTS
#endif

#ifdef SWEEPER_KBA
#include "sweeper_kba_c_kernels.h"
#endif

/*---------------------------------------------------------------------------*/
