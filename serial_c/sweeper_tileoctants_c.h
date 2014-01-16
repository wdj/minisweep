/*---------------------------------------------------------------------------*/
/*!
 * \file   sweeper_tileoctnats_c.h
 * \author Wayne Joubert
 * \date   Wed Jan 15 16:06:28 EST 2014
 * \brief  Definitions for performing a sweep, tileoctants version.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _serial_c__sweeper_tileoctants_c_h_
#define _serial_c__sweeper_tileoctants_c_h_

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
                   int         tile_octants )
{
  sweeper->tile_octants = tile_octants;

  /*---Allocate arrays---*/

  sweeper->v_local = pmalloc( dims.na * NU );
  sweeper->facexy  = pmalloc( dims.nx * dims.ny * dims.ne * dims.na * NU
                                   * ( sweeper->tile_octants ? NOCTANT : 1 ) );
  sweeper->facexz  = pmalloc( dims.nx * dims.nz * dims.ne * dims.na * NU
                                   * ( sweeper->tile_octants ? NOCTANT : 1 ) );
  sweeper->faceyz  = pmalloc( dims.ny * dims.nz * dims.ne * dims.na * NU
                                   * ( sweeper->tile_octants ? NOCTANT : 1 ) );
}

/*===========================================================================*/
/*---Pseudo-destructor for Sweeper struct---*/

void Sweeper_dtor( Sweeper* sweeper )
{
  /*---Deallocate arrays---*/

  pfree( sweeper->v_local );
  pfree( sweeper->facexy );
  pfree( sweeper->facexz );
  pfree( sweeper->faceyz );
}

/*===========================================================================*/
/*---Perform a sweep---*/

void Sweeper_sweep(
  Sweeper*         sweeper,
  P* __restrict__  vo,
  P* __restrict__  vi,
  Quantities       quan,
  Dimensions       dims )
{
  assert( sweeper );
  assert( vi );
  assert( vo );

  /*---Declarations---*/
  int tile_step = 0;

  /*---Initialize result array to zero---*/

  initialize_state_zero( vo, dims );

  /*---Loop over octant tiles---*/

  /*---If tiling is requested, at each of the 8 tile steps,
       each octant direction computes its result on a
       1/8-sized tile of the domain.
       This is scheduled such that: 1) the proper sweep order for each
       octant direction is adhered to, and 2) for each tile step, the
       8 octant directions are working on independent disjoint
       tiles of the domain.
  ---*/

  for( tile_step=0; tile_step<(sweeper->tile_octants?NOCTANT:1); ++tile_step )
  {
    int ix = 0;
    int iy = 0;
    int iz = 0;
    int ie = 0;
    int im = 0;
    int ia = 0;
    int iu = 0;
    int ioctant = 0;

    /*---Loop over octants---*/

    for( ioctant=0; ioctant<NOCTANT; ++ioctant )
    {
      /*---If tiling, then each octant direction needs its own face:
           each octant direction is not processed all-at-once but
           intermittently, thus a need to remember its state.
      ---*/

      const int octant_index = sweeper->tile_octants ? ioctant : 0;

      /*---Decode octant directions from octant number---*/
      /*--- -1 for downward direction, +1 for upward direction---*/

      const int idirx = ioctant & (1<<0) ? -1 : 1;
      const int idiry = ioctant & (1<<1) ? -1 : 1;
      const int idirz = ioctant & (1<<2) ? -1 : 1;

      /*---Determine tile to be computed---*/

      /*---The scheduling works right if each octant direction visits tiles
           in lexicographical ordering x-varies-fastest, starting at
           the origin corner of the octant direction.
           Note the tiles are defined by lower or upper half of domain
           in each of x, y and z.

           NOTE: this strategy is tested for the single thread case
           but not yet tested under OpenMP.
      ---*/

      const int tile_x = (!sweeper->tile_octants) ? 0 :
                 ( ( tile_step & (1<<0) ) == 0 ) == ( idirx == +1 ) ? -1 : 1;
      const int tile_y = (!sweeper->tile_octants) ? 0 :
                 ( ( tile_step & (1<<1) ) == 0 ) == ( idiry == +1 ) ? -1 : 1;
      const int tile_z = (!sweeper->tile_octants) ? 0 :
                 ( ( tile_step & (1<<2) ) == 0 ) == ( idirz == +1 ) ? -1 : 1;

      /*---Compute tile boundaries---*/

      /*---If no tiling, then whole domain, otherwise 1/2 of
           domain in each direction
      ---*/

      const int tile_xmin = (!sweeper->tile_octants) ? 0 :
                            tile_x==-1           ? 0 : dims.nx/2;
      const int tile_ymin = (!sweeper->tile_octants) ? 0 :
                            tile_y==-1           ? 0 : dims.ny/2;
      const int tile_zmin = (!sweeper->tile_octants) ? 0 :
                            tile_z==-1           ? 0 : dims.nz/2;

      const int tile_xmax = (!sweeper->tile_octants) ? dims.nx   :
                            tile_x==-1           ? dims.nx/2 : dims.nx;
      const int tile_ymax = (!sweeper->tile_octants) ? dims.ny   :
                            tile_y==-1           ? dims.ny/2 : dims.ny;
      const int tile_zmax = (!sweeper->tile_octants) ? dims.nz   :
                            tile_z==-1           ? dims.nz/2 : dims.nz;

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

      if( tile_z != idirz || !sweeper->tile_octants )
      {
        iz = idirz==+1 ? -1 : dims.nz;
        for( iu=0; iu<NU; ++iu )
        for( iy=tile_ymin; iy<tile_ymax; ++iy )
        for( ix=tile_xmin; ix<tile_xmax; ++ix )
        for( ie=0; ie<dims.ne; ++ie )
        for( ia=0; ia<dims.na; ++ia )
        {
          *ref_facexy( sweeper->facexy, dims, ix, iy, ie, ia, iu, octant_index )
                            = Quantities_init_facexy( ix, iy, iz, ie, ia, iu );
        }
      }

      if( tile_y != idiry || !sweeper->tile_octants )
      {
        iy = idiry==+1 ? -1 : dims.ny;
        for( iu=0; iu<NU; ++iu )
        for( iz=tile_zmin; iz<tile_zmax; ++iz )
        for( ix=tile_xmin; ix<tile_xmax; ++ix )
        for( ie=0; ie<dims.ne; ++ie )
        for( ia=0; ia<dims.na; ++ia )
        {
          *ref_facexz( sweeper->facexz, dims, ix, iz, ie, ia, iu, octant_index )
                            = Quantities_init_facexz( ix, iy, iz, ie, ia, iu );
        }
      }

      if( tile_x != idirx || !sweeper->tile_octants )
      {
        ix = idirx==+1 ? -1 : dims.nx;
        for( iu=0; iu<NU; ++iu )
        for( iz=tile_zmin; iz<tile_zmax; ++iz )
        for( iy=tile_ymin; iy<tile_ymax; ++iy )
        for( ie=0; ie<dims.ne; ++ie )
        for( ia=0; ia<dims.na; ++ia )
        {
          *ref_faceyz( sweeper->faceyz, dims, iy, iz, ie, ia, iu, octant_index )
                            = Quantities_init_faceyz( ix, iy, iz, ie, ia, iu );
        }
      }

      /*---Loop over energy groups---*/

      for( ie=0; ie<dims.ne; ++ie )
      {
        /*---Calculate spatial loop extents, possibly based on tiling---*/

        const int ixbeg = idirx==+1 ? tile_xmin : tile_xmax-1;
        const int iybeg = idiry==+1 ? tile_ymin : tile_ymax-1;
        const int izbeg = idirz==+1 ? tile_zmin : tile_zmax-1;

        const int ixend = idirx==-1 ? tile_xmin : tile_xmax-1;
        const int iyend = idiry==-1 ? tile_ymin : tile_ymax-1;
        const int izend = idirz==-1 ? tile_zmin : tile_zmax-1;

        /*---Loop over gridcells, in proper direction---*/

        for( iz=izbeg; iz!=izend+idirz; iz+=idirz )
        for( iy=iybeg; iy!=iyend+idiry; iy+=idiry )
        for( ix=ixbeg; ix!=ixend+idirx; ix+=idirx )
        {

          /*---Transform state vector from moments to angles---*/

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
                        *ref_state( vi, dims, ix, iy, iz, ie, im, iu );
            }
            *ref_v_local( sweeper->v_local, dims, ia, iu ) = result;
          }

          /*---Perform solve---*/

          Quantities_solve( sweeper->v_local,
                            sweeper->facexy, sweeper->facexz, sweeper->faceyz,
                            ix, iy, iz, ie, octant_index, quan, dims );

          /*---Transform state vector from angles to moments---*/

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
                        *ref_v_local( sweeper->v_local, dims, ia, iu );
            }
            *ref_state( vo, dims, ix, iy, iz, ie, im, iu ) += result;
          }

        } /*---ix/iy/iz---*/

      } /*---ie---*/

    } /*---ioctant---*/

  } /*---octant_tile---*/

} /*---sweep---*/

/*===========================================================================*/

#endif /*---_serial_c__sweeper_tileoctants_c_h_---*/

/*===========================================================================*/
