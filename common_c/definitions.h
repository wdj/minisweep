/*---------------------------------------------------------------------------*/
/*!
 * \file   definitions.h
 * \author Wayne Joubert
 * \date   Wed Jan 15 16:06:28 EST 2014
 * \brief  Basic definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _common_c__definitions_h_
#define _common_c__definitions_h_

/*===========================================================================*/

/*---Initializations---*/

/*---Number of physical dimensions---*/
enum{ NDIM = 3 };

/*---Number of octant directions---*/
enum{ NOCTANT = 8 };

/*===========================================================================*/
/*---Functions to manipulate sweep directions---*/

static inline int Dir_up() { return +1; }
static inline int Dir_dn() { return -1; }

static inline int Dir_hi() { return +1; }
static inline int Dir_lo() { return -1; }

static inline int Dir_x( int octant ) { return octant & (1<<0) ? Dir_dn()
                                                               : Dir_up(); }
static inline int Dir_y( int octant ) { return octant & (1<<1) ? Dir_dn()
                                                               : Dir_up(); }
static inline int Dir_z( int octant ) { return octant & (1<<2) ? Dir_dn()
                                                               : Dir_up(); }

static inline int Dir_inc( int dir ) { return dir; }

/*===========================================================================*/

#endif /*---_common_c__definitions_h_---*/

/*---------------------------------------------------------------------------*/
