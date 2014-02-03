/*---------------------------------------------------------------------------*/
/*!
 * \file   sweeper_simple_c.h
 * \author Wayne Joubert
 * \date   Wed Jan 15 16:06:28 EST 2014
 * \brief  Definitions for performing a sweep, simple version.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _serial_c__sweeper_simple_c_h_
#define _serial_c__sweeper_simple_c_h_

#include "env.h"
#include "definitions.h"
#include "quantities.h"
#include "array_accessors.h"
#include "array_operations.h"
#include "memory.h"
#include "sweeper_simple.h"

/*===========================================================================*/
/*---Pseudo-constructor for Sweeper struct---*/

void Sweeper_ctor( Sweeper*    sweeper,
                   Dimensions  dims,
                   Env*        env,
                   int         nblock_z )
{
  /*---Allocate arrays---*/

  sweeper->v_local = malloc_P( dims.na * NU );
  sweeper->facexy  = malloc_P( dims.nx * dims.ny * dims.ne * dims.na *
                                  NU * Sweeper_num_face_octants_allocated() );
  sweeper->facexz  = malloc_P( dims.nx * dims.nz * dims.ne * dims.na *
                                  NU * Sweeper_num_face_octants_allocated() );
  sweeper->faceyz  = malloc_P( dims.ny * dims.nz * dims.ne * dims.na *
                                  NU * Sweeper_num_face_octants_allocated() );
}

/*===========================================================================*/
/*---Pseudo-destructor for Sweeper struct---*/

void Sweeper_dtor( Sweeper* sweeper )
{
  /*---Deallocate arrays---*/

  free_P( sweeper->v_local );
  free_P( sweeper->facexy );
  free_P( sweeper->facexz );
  free_P( sweeper->faceyz );

  sweeper->v_local = NULL;
  sweeper->facexy  = NULL;
  sweeper->facexz  = NULL;
  sweeper->faceyz  = NULL;
}

/*===========================================================================*/
/*---Perform a sweep---*/

void Sweeper_sweep(
  Sweeper*               sweeper,
  P* __restrict__        vo,
  const P* __restrict__  vi,
  const Quantities*      quan,
  Dimensions             dims,
  Env*                   env )
{
  assert( sweeper );
  assert( vi );
  assert( vo );

  /*---Declarations---*/
  int ix = 0;
  int iy = 0;
  int iz = 0;
  int ie = 0;
  int im = 0;
  int ia = 0;
  int iu = 0;
  int octant = 0;

  /*---Initialize result array to zero---*/

  initialize_state_zero( vo, dims, NU );

  /*---Loop over octants---*/

  for( octant=0; octant<NOCTANT; ++octant )
  {
    const int octant_ind = 0;
    assert( octant_ind >= 0 &&
            octant_ind < Sweeper_num_face_octants_allocated() );

    /*---Decode octant directions from octant number---*/

    const int dir_x = Dir_x( octant );
    const int dir_y = Dir_y( octant );
    const int dir_z = Dir_z( octant );

    /*---Initialize faces---*/

    /*---The semantics of the face arrays are as follows.
         On entering a cell for a solve at the gridcell level,
         the face array is assumed to have a value corresponding to
         "one cell lower" in the relevant direction.
         On leaving the gridcell solve, the face has been updated
         to have the flux at that gridcell.
         Thus, the face is initialized at first to have a value
         "one cell" outside of the domain, e.g., for the XY face,
         either -1 or dims.nx.
         Note also that the face initializer functions now take
         coordinates for all three spatial dimensions --
         the third dimension is used to denote whether it is the
         "lower" or "upper" face and also its exact location
         in that dimension.
    ---*/

    {
      iz = dir_z == Dir_up() ? -1 : dims.nz;
      for( iu=0; iu<NU; ++iu )
      for( iy=0; iy<dims.ny; ++iy )
      for( ix=0; ix<dims.nx; ++ix )
      for( ie=0; ie<dims.ne; ++ie )
      for( ia=0; ia<dims.na; ++ia )
      {
        *ref_facexy( sweeper->facexy, dims, NU,
                     Sweeper_num_face_octants_allocated(),
                     ix, iy, ie, ia, iu, octant_ind )
             = Quantities_init_facexy(
                                  quan, ix, iy, iz, ie, ia, iu, octant, dims );
      }
    }

    {
      iy = dir_y == Dir_up() ? -1 : dims.ny;
      for( iu=0; iu<NU; ++iu )
      for( iz=0; iz<dims.nz; ++iz )
      for( ix=0; ix<dims.nx; ++ix )
      for( ie=0; ie<dims.ne; ++ie )
      for( ia=0; ia<dims.na; ++ia )
      {
        *ref_facexz( sweeper->facexz, dims, NU,
                     Sweeper_num_face_octants_allocated(),
                     ix, iz, ie, ia, iu, octant_ind )
             = Quantities_init_facexz(
                                  quan, ix, iy, iz, ie, ia, iu, octant, dims );
      }
    }

    {
      ix = dir_x == Dir_up() ? -1 : dims.nx;
      for( iu=0; iu<NU; ++iu )
      for( iz=0; iz<dims.nz; ++iz )
      for( iy=0; iy<dims.ny; ++iy )
      for( ie=0; ie<dims.ne; ++ie )
      for( ia=0; ia<dims.na; ++ia )
      {
        *ref_faceyz( sweeper->faceyz, dims, NU,
                     Sweeper_num_face_octants_allocated(),
                     iy, iz, ie, ia, iu, octant_ind )
             = Quantities_init_faceyz(
                                  quan, ix, iy, iz, ie, ia, iu, octant, dims );
      }
    }

    /*---Loop over energy groups---*/

    for( ie=0; ie<dims.ne; ++ie )
    {
      /*---Calculate spatial loop extents---*/

      int ixbeg = dir_x==Dir_up() ? 0 : dims.nx-1;
      int iybeg = dir_y==Dir_up() ? 0 : dims.ny-1;
      int izbeg = dir_z==Dir_up() ? 0 : dims.nz-1;

      int ixend = dir_x==Dir_dn() ? 0 : dims.nx-1;
      int iyend = dir_y==Dir_dn() ? 0 : dims.ny-1;
      int izend = dir_z==Dir_dn() ? 0 : dims.nz-1;

      /*---Loop over gridcells, in proper direction---*/

    for( iz=izbeg; iz!=izend+Dir_inc(dir_z); iz+=Dir_inc(dir_z) )
    for( iy=iybeg; iy!=iyend+Dir_inc(dir_y); iy+=Dir_inc(dir_y) )
    for( ix=ixbeg; ix!=ixend+Dir_inc(dir_x); ix+=Dir_inc(dir_x) )
    {

      /*--------------------*/
      /*---Transform state vector from moments to angles---*/
      /*--------------------*/

      /*---This loads values from the input state vector,
           does the small dense matrix-vector product,
           and stores the result in a relatively small local
           array that is hopefully small enough to fit into
           processor cache.
      ---*/

      for( iu=0; iu<NU; ++iu )
      for( ia=0; ia<dims.na; ++ia )
      {
        P result = P_zero();
        for( im=0; im<dims.nm; ++im )
        {
          result += *const_ref_a_from_m( quan->a_from_m, dims, im, ia ) *
                    *const_ref_state( vi, dims, NU, ix, iy, iz, ie, im, iu );
        }
        *ref_v_local( sweeper->v_local, dims, NU, ia, iu ) = result;
      }

      /*--------------------*/
      /*---Perform solve---*/
      /*--------------------*/

      Quantities_solve( quan, sweeper->v_local,
                        sweeper->facexy, sweeper->facexz, sweeper->faceyz,
                        ix, iy, iz, ie, ix, iy, iz,
                        octant, octant_ind,
                        Sweeper_num_face_octants_allocated(), dims, dims );

      /*--------------------*/
      /*---Transform state vector from angles to moments---*/
      /*--------------------*/

      /*---Perform small dense matrix-vector products and store
           the result in the output state vector.
      ---*/

      for( iu=0; iu<NU; ++iu )
      for( im=0; im<dims.nm; ++im )
      {
        P result = P_zero();
        for( ia=0; ia<dims.na; ++ia )
        {
          result += *const_ref_m_from_a( quan->m_from_a, dims, im, ia ) *
                    *const_ref_v_local( sweeper->v_local, dims, NU, ia, iu );
        }
        *ref_state( vo, dims, NU, ix, iy, iz, ie, im, iu ) += result;
      }

    } /*---ix/iy/iz---*/

    } /*---ie---*/

  } /*---octant---*/

} /*---sweep---*/

/*===========================================================================*/

#endif /*---_serial_c__sweeper_simple_c_h_---*/

/*---------------------------------------------------------------------------*/
