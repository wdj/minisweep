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

/*---------------------------------------------------------------------------*/
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

/*---------------------------------------------------------------------------*/

#endif /*---_common_c__dimensions_h_---*/

/*---------------------------------------------------------------------------*/
