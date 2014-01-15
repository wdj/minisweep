/*---------------------------------------------------------------------------*/
/*!
 * \file   array_accessors.h
 * \author Wayne Joubert
 * \date   Wed Jan 15 16:06:28 EST 2014
 * \brief  Functions for referencing special-purpose multidimensional arrays.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _serial_c__array_accessors_h_
#define _serial_c__array_accessors_h_

#include "definitions.h"
#include "dimensions.h"

/*---------------------------------------------------------------------------*/
/*---Multidimensional array accessor function---*/

static inline P* ref_state(
    P* __restrict__ v,
    Dimensions dims,
    int ix,
    int iy,
    int iz,
    int ie,
    int im,
    int iu )
{
  return & v[
              im + dims.nm * (
              iu + NU      * (
              ix + dims.nx * (
              iy + dims.ny * (
              iz + dims.nz * (
              ie + dims.ne * (
              0 )))))) ];
}

/*---------------------------------------------------------------------------*/
/*---Multidimensional array accessor function---*/

static inline P* ref_v_local(
    P* __restrict__ v,
    Dimensions dims,
    int ia,
    int iu,
    int ie,
    int ioctant,
    int ix,
    int iy,
    int iz)
{
  return & v[
              ia + dims.na * (
              iu + NU      * (
              0 )) ];
}

/*---------------------------------------------------------------------------*/
/*---Multidimensional array accessor function---*/

static inline P* ref_a_from_m(
    P* __restrict__ v,
    Dimensions dims,
    int im,
    int ia )
{
  return & v[
              im + dims.nm * (
              ia + dims.na * (
              0 )) ];
}

/*---------------------------------------------------------------------------*/
/*---Multidimensional array accessor function---*/

static inline P* ref_m_from_a(
    P* __restrict__ v,
    Dimensions dims,
    int im,
    int ia )
{
  return & v[
              ia + dims.na * (
              im + dims.nm * (
              0 )) ];
}

/*---------------------------------------------------------------------------*/
/*---Multidimensional array accessor function---*/

static inline P* ref_facexy(
    P* __restrict__ v,
    Dimensions dims,
    int ix,
    int iy,
    int ie,
    int ia,
    int iu,
    int ioctant )
{
  return & v[
              ia + dims.na * (
              iu + NU      * (
              ix + dims.nx * (
              iy + dims.ny * (
              ie + dims.ne * (
              0 ))))) ];
}

/*---------------------------------------------------------------------------*/
/*---Multidimensional array accessor function---*/

static inline P* ref_facexz(
    P* __restrict__ v,
    Dimensions dims,
    int ix,
    int iz,
    int ie,
    int ia,
    int iu,
    int ioctant )
{
  return & v[
              ia + dims.na * (
              iu + NU      * (
              ix + dims.nx * (
              iz + dims.nz * (
              ie + dims.ne * (
              0 ))))) ];
}

/*---------------------------------------------------------------------------*/
/*---Multidimensional array accessor function---*/

static inline P* ref_faceyz(
    P* __restrict__ v,
    Dimensions dims,
    int iy,
    int iz,
    int ie,
    int ia,
    int iu,
    int ioctant )
{
  return & v[
              ia + dims.na * (
              iu + NU      * (
              iy + dims.ny * (
              iz + dims.nz * (
              ie + dims.ne * (
              0 ))))) ];
}

/*---------------------------------------------------------------------------*/

#endif /*---_serial_c__array_accessors_h_---*/

/*---------------------------------------------------------------------------*/
