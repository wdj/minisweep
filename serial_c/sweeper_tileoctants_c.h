/*---------------------------------------------------------------------------*/
/*!
 * \file   sweeper_tileoctants_c.h
 * \author Wayne Joubert
 * \date   Wed Jan 15 16:06:28 EST 2014
 * \brief  Definitions for performing a sweep, tileoctants version.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _serial_c__sweeper_tileoctants_c_h_
#define _serial_c__sweeper_tileoctants_c_h_

#include "env.h"
#include "definitions.h"
#include "quantities.h"
#include "array_accessors.h"
#include "array_operations.h"
#include "memory.h"
#include "sweeper_tileoctants.h"

/*===========================================================================*/
/*---Pseudo-constructor for Sweeper struct---*/

void Sweeper_ctor( Sweeper*    sweeper,
                   Dimensions  dims,
                   Env*        env,
                   int         nblock_z )
{
  const Bool_t tile_octants = Sweeper_tile_octants();

  /*---Allocate arrays---*/

  sweeper->v_local = malloc_P( dims.na * NU );
  sweeper->facexy  = malloc_P( dims.nx * dims.ny * dims.ne * dims.na * 
                                             NU * Sweeper_num_face_octants() );
  sweeper->facexz  = malloc_P( dims.nx * dims.nz * dims.ne * dims.na * 
                                             NU * Sweeper_num_face_octants() );
  sweeper->faceyz  = malloc_P( dims.ny * dims.nz * dims.ne * dims.na * 
                                             NU * Sweeper_num_face_octants() );
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
  Sweeper*         sweeper,
  P* __restrict__  vo,
  P* __restrict__  vi,
  Quantities       quan,
  Dimensions       dims,
  Env*             env )
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
  const Bool_t tile_octants = Sweeper_tile_octants();

  int tile_step = 0;

  /*---Initialize result array to zero---*/

  initialize_state_zero( vo, dims, NU );

  /*---Loop over octant tiles---*/

  /*---If tiling is requested, at each of the 8 tile steps,
       each octant direction computes its result on a
       1/8-sized tile of the domain.
       This is scheduled such that: 1) the proper sweep order for each
       octant direction is adhered to, and 2) for each tile step, the
       8 octant directions are working on independent disjoint
       tiles of the domain.
  ---*/

  for( tile_step=0; tile_step<(tile_octants?NOCTANT:1); ++tile_step )
  {

  /*---Loop over octants---*/

  for( octant=0; octant<NOCTANT; ++octant )
  {
    /*---If tiling, then each octant direction needs its own face:
         each octant direction is not processed all-at-once but
         intermittently, thus a need to remember its state.
    ---*/

    const int octant_ind = tile_octants ? octant : 0;
    assert( octant_ind >= 0 && octant_ind < Sweeper_num_face_octants() );

    /*---Decode octant directions from octant number---*/

    const int dir_x = Dir_x( octant );
    const int dir_y = Dir_y( octant );
    const int dir_z = Dir_z( octant );

    /*---Determine tile to be computed---*/

    /*---The scheduling works right if each octant direction visits tiles
         in lexicographical ordering x-varies-fastest, starting at
         the origin corner of the octant direction.
         Note the tiles are defined by lower or upper half of domain
         in each of x, y and z.

         NOTE: this strategy is tested for the single thread case
         but not yet tested under OpenMP.
    ---*/

    const int tile_x = (!tile_octants) ? 0 :
           ( ( tile_step & (1<<0) ) == 0 ) == ( dir_x == Dir_up() ) ? Dir_lo()
                                                                    : Dir_hi();
    const int tile_y = (!tile_octants) ? 0 :
           ( ( tile_step & (1<<1) ) == 0 ) == ( dir_y == Dir_up() ) ? Dir_lo()
                                                                    : Dir_hi();
    const int tile_z = (!tile_octants) ? 0 :
           ( ( tile_step & (1<<2) ) == 0 ) == ( dir_z == Dir_up() ) ? Dir_lo()
                                                                    : Dir_hi();

    /*---Compute tile boundaries---*/

    /*---If no tiling, then whole domain, otherwise 1/2 of
         domain in each direction
    ---*/

    const int tile_xmin = (!tile_octants)  ? 0         :
                          tile_x==Dir_lo() ? 0         : dims.nx/2;
    const int tile_ymin = (!tile_octants)  ? 0         :
                          tile_y==Dir_lo() ? 0         : dims.ny/2;
    const int tile_zmin = (!tile_octants)  ? 0         :
                          tile_z==Dir_lo() ? 0         : dims.nz/2;

    const int tile_xmax = (!tile_octants)  ? dims.nx   :
                          tile_x==Dir_lo() ? dims.nx/2 : dims.nx;
    const int tile_ymax = (!tile_octants)  ? dims.ny   :
                          tile_y==Dir_lo() ? dims.ny/2 : dims.ny;
    const int tile_zmax = (!tile_octants)  ? dims.nz   :
                          tile_z==Dir_lo() ? dims.nz/2 : dims.nz;

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
         Note if no tiling then we initialize all faces here,
         otherwise only the appropriate part on the appropriate
         tiling step.
    ---*/

    if( tile_z != dir_z || !tile_octants )
    {
      iz = dir_z==Dir_up() ? -1 : dims.nz;
      for( iu=0; iu<NU; ++iu )
      for( iy=tile_ymin; iy<tile_ymax; ++iy )
      for( ix=tile_xmin; ix<tile_xmax; ++ix )
      for( ie=0; ie<dims.ne; ++ie )
      for( ia=0; ia<dims.na; ++ia )
      {
        *ref_facexy( sweeper->facexy, dims, NU, ix, iy, ie, ia, iu, octant_ind )
                    = Quantities_init_facexy( ix, iy, iz, ie, ia, iu, dims );
      }
    }

    if( tile_y != dir_y || !tile_octants )
    {
      iy = dir_y==Dir_up() ? -1 : dims.ny;
      for( iu=0; iu<NU; ++iu )
      for( iz=tile_zmin; iz<tile_zmax; ++iz )
      for( ix=tile_xmin; ix<tile_xmax; ++ix )
      for( ie=0; ie<dims.ne; ++ie )
      for( ia=0; ia<dims.na; ++ia )
      {
        *ref_facexz( sweeper->facexz, dims, NU, ix, iz, ie, ia, iu, octant_ind )
                    = Quantities_init_facexz( ix, iy, iz, ie, ia, iu, dims );
      }
    }

    if( tile_x != dir_x || !tile_octants )
    {
      ix = dir_x==Dir_up() ? -1 : dims.nx;
      for( iu=0; iu<NU; ++iu )
      for( iz=tile_zmin; iz<tile_zmax; ++iz )
      for( iy=tile_ymin; iy<tile_ymax; ++iy )
      for( ie=0; ie<dims.ne; ++ie )
      for( ia=0; ia<dims.na; ++ia )
      {
        *ref_faceyz( sweeper->faceyz, dims, NU, iy, iz, ie, ia, iu, octant_ind )
                    = Quantities_init_faceyz( ix, iy, iz, ie, ia, iu, dims );
      }
    }

    /*---Loop over energy groups---*/

    for( ie=0; ie<dims.ne; ++ie )
    {
      /*---Calculate spatial loop extents, possibly based on tiling---*/

      const int ixbeg = dir_x==Dir_up() ? tile_xmin : tile_xmax-1;
      const int iybeg = dir_y==Dir_up() ? tile_ymin : tile_ymax-1;
      const int izbeg = dir_z==Dir_up() ? tile_zmin : tile_zmax-1;

      const int ixend = dir_x==Dir_dn() ? tile_xmin : tile_xmax-1;
      const int iyend = dir_y==Dir_dn() ? tile_ymin : tile_ymax-1;
      const int izend = dir_z==Dir_dn() ? tile_zmin : tile_zmax-1;

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
          result += *ref_a_from_m( quan.a_from_m, dims, im, ia ) *
                    *ref_state( vi, dims, NU, ix, iy, iz, ie, im, iu );
        }
        *ref_v_local( sweeper->v_local, dims, NU, ia, iu ) = result;
      }

      /*--------------------*/
      /*---Perform solve---*/
      /*--------------------*/

      Quantities_solve( sweeper->v_local,
                        sweeper->facexy, sweeper->facexz, sweeper->faceyz,
                        ix, iy, iz, ie, octant, octant_ind, quan, dims );

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
          result += *ref_m_from_a( quan.m_from_a, dims, im, ia ) *
                    *ref_v_local( sweeper->v_local, dims, NU, ia, iu );
        }
        *ref_state( vo, dims, NU, ix, iy, iz, ie, im, iu ) += result;
      }

    } /*---ix/iy/iz---*/

    } /*---ie---*/

  } /*---octant---*/

  } /*---octant_tile---*/

} /*---sweep---*/

/*===========================================================================*/

#endif /*---_serial_c__sweeper_tileoctants_c_h_---*/

/*---------------------------------------------------------------------------*/
