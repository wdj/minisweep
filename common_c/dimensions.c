/*---------------------------------------------------------------------------*/
/*!
 * \file   dimensions.c
 * \author Wayne Joubert
 * \date   Fri Jan 17 12:21:18 EST 2014
 * \brief  Definitions for Dimensions struct.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include <stdlib.h>

#include "definitions.h"
#include "dimensions.h"

/*===========================================================================*/
/*---Size of state vector---*/

size_t Dimensions_size_state( const Dimensions dims, int nu )
{
  return ( (size_t)dims.nx )
       * ( (size_t)dims.ny )
       * ( (size_t)dims.nz )
       * ( (size_t)dims.ne )
       * ( (size_t)dims.nm )
       * ( (size_t)nu );
}

/*===========================================================================*/
/*---Size of state vector in angles space---*/

size_t Dimensions_size_state_angles( const Dimensions dims, int nu )
{
  return ( (size_t)dims.nx )
       * ( (size_t)dims.ny )
       * ( (size_t)dims.nz )
       * ( (size_t)dims.ne )
       * ( (size_t)dims.na )
       * ( (size_t)nu )
       * ( (size_t)NOCTANT );
}

/*===========================================================================*/
/*---Size of face vectors---*/

size_t Dimensions_size_facexy( const Dimensions dims, 
                               int nu, 
                               int num_face_octants_allocated )
{
  return ( (size_t)dims.nx )
       * ( (size_t)dims.ny )
       * ( (size_t)dims.ne )
       * ( (size_t)dims.na )
       * ( (size_t)nu )
       * ( (size_t)num_face_octants_allocated );
}

/*---------------------------------------------------------------------------*/

size_t Dimensions_size_facexz( const Dimensions dims, 
                               int nu, 
                               int num_face_octants_allocated )
{
  return ( (size_t)dims.nx )
       * ( (size_t)dims.nz )
       * ( (size_t)dims.ne )
       * ( (size_t)dims.na )
       * ( (size_t)nu )
       * ( (size_t)num_face_octants_allocated );
}

/*---------------------------------------------------------------------------*/

size_t Dimensions_size_faceyz( const Dimensions dims, 
                               int nu, 
                               int num_face_octants_allocated )
{
  return ( (size_t)dims.ny )
       * ( (size_t)dims.nz )
       * ( (size_t)dims.ne )
       * ( (size_t)dims.na )
       * ( (size_t)nu )
       * ( (size_t)num_face_octants_allocated );
}

/*===========================================================================*/

/*---------------------------------------------------------------------------*/
