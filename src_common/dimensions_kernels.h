/*---------------------------------------------------------------------------*/
/*!
 * \file   dimensions_kernels.h
 * \author Wayne Joubert
 * \date   Wed Jan 15 16:06:28 EST 2014
 * \brief  Declarations for Dimensions struct.  Code for device kernels.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _dimensions_kernels_h_
#define _dimensions_kernels_h_

#include "types_kernels.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*===========================================================================*/
/*---Enums for compile-time sizes---*/

/*---Number of unknowns per gridcell---*/
#ifdef NU_VALUE
enum{ NU = NU_VALUE };
#else
enum{ NU = 4 };
#endif

/*---Number of moments---*/
#ifdef NM_VALUE
enum{ NM = NM_VALUE };
#else
enum{ NM = 16 };
#endif

/*===========================================================================*/
/*---Struct to hold problem dimensions---*/

typedef struct
{
  /*---Grid spatial dimensions---*/
  int nx;
  int ny;
  int nz;

  /*---Number of energy groups---*/
  int ne;

  /*----Number of moments---*/
  int nm;

  /*---Number of angles---*/
  int na;
} Dimensions;

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

#endif /*---_dimensions_kernels_h_---*/

/*---------------------------------------------------------------------------*/
