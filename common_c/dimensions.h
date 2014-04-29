/*---------------------------------------------------------------------------*/
/*!
 * \file   dimensions.h
 * \author Wayne Joubert
 * \date   Wed Jan 15 16:06:28 EST 2014
 * \brief  Declarations for Dimensions struct.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _common_c__dimensions_h_
#define _common_c__dimensions_h_

#include <stdlib.h>

#include "env.h"
#include "types.h"

#ifdef __cplusplus
extern "C"
{
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
/*---Size of state vector---*/

size_t Dimensions_size_state( const Dimensions dims, int nu );

/*===========================================================================*/
/*---Size of state vector in angles space---*/

size_t Dimensions_size_state_angles( const Dimensions dims, int nu );

/*===========================================================================*/
/*---Size of face vectors---*/

size_t Dimensions_size_facexy( const Dimensions dims,
                               int nu,
                               int num_face_octants_allocated );

/*---------------------------------------------------------------------------*/

size_t Dimensions_size_facexz( const Dimensions dims,
                               int nu,
                               int num_face_octants_allocated );

/*---------------------------------------------------------------------------*/

size_t Dimensions_size_faceyz( const Dimensions dims,
                               int nu,
                               int num_face_octants_allocated );

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

#endif /*---_common_c__dimensions_h_---*/

/*---------------------------------------------------------------------------*/
