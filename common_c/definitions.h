/*---------------------------------------------------------------------------*/
/*!
 * \file   defs.h
 * \author Wayne Joubert
 * \date   Wed Jan 15 16:06:28 EST 2014
 * \brief  Basic definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _common_c__defs_h_
#define _common_c__defs_h_

/*---------------------------------------------------------------------------*/

/*---Floating point type for sweep---*/

typedef double P;

/*---Initializations---*/

/*---Number of physical dimensions---*/
enum{ NDIM = 3 };

/*---Number of octant directions---*/
enum{ NOCTANT = 8 };

/*---Number of unknowns per gridcell---*/
/*---Currently hardwired, this could be changed in the future---*/
enum{ NU = 1 };

/*---------------------------------------------------------------------------*/

#endif /*---_common_c__defs_h_---*/

/*---------------------------------------------------------------------------*/
