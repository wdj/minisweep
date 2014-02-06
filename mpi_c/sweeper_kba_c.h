/*---------------------------------------------------------------------------*/
/*!
 * \file   sweeper_kba_c.h
 * \author Wayne Joubert
 * \date   Tue Jan 28 16:37:41 EST 2014
 * \brief  Definitions for performing a sweep, kba version.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _serial_c__sweeper_kba_c_h_
#define _serial_c__sweeper_kba_c_h_

#include "env.h"
#include "definitions.h"
#include "quantities.h"
#include "array_accessors.h"
#include "array_operations.h"
#include "memory.h"
#include "sweeper_kba.h"

/*===========================================================================*/
/*---Pseudo-constructor for Sweeper struct---*/

void Sweeper_ctor( Sweeper*    sweeper,
                   Dimensions  dims,
                   Env*        env,
                   int         nblock_z )
{
  Insist( dims.nz % nblock_z == 0 &&
          "KBA sweeper currently requires all blocks have same z dimension" );
  Insist( dims.nx > 0 && "KBA sweeper currently requires all blocks nonempty" );
  Insist( dims.ny > 0 && "KBA sweeper currently requires all blocks nonempty" );
  Insist( dims.nz > 0 && "KBA sweeper currently requires all blocks nonempty" );

  /*---Set up dimensions of kba block---*/
  sweeper->nblock_z = nblock_z;
  sweeper->dims_b = dims;
  sweeper->dims_b.nz = dims.nz / nblock_z;

  /*---Allocate arrays---*/

  sweeper->v_local = malloc_P( sweeper->dims_b.na * NU );
  sweeper->facexy  = malloc_P( Dimensions_size_facexy( sweeper->dims_b, NU,
                                      Sweeper_num_face_octants_allocated() ) );
  sweeper->facexz  = malloc_P( Dimensions_size_facexz( sweeper->dims_b, NU,
                                      Sweeper_num_face_octants_allocated() ) );
  sweeper->faceyz  = malloc_P( Dimensions_size_faceyz( sweeper->dims_b, NU,
                                      Sweeper_num_face_octants_allocated() ) );
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
/*---Number of block steps executed for an octant in isolation---*/

int Sweeper_nblock( const Sweeper*  sweeper,
                    const Env*      env )
{
  return sweeper->nblock_z;
}

/*===========================================================================*/
/*---Number of kba parallel steps---*/

int Sweeper_nstep( const Sweeper*  sweeper,
                   const Env*      env )
{
  return NOCTANT * Sweeper_nblock( sweeper, env )
                                              + 3 * ( Env_nproc_x( env ) - 1 )
                                              + 2 * ( Env_nproc_y( env ) - 1 );
}

/*===========================================================================*/
/*---Wehther this proc active for a given sweep step---*/

Step_Info Sweeper_step_info__( const Sweeper*  sweeper,
                               int             step,
                               int             proc_x,
                               int             proc_y,
                               const Env*      env )
{
/*
  assert( step >= 0 && step < Sweeper_nstep( sweeper, env ) );
  assert( proc_x >= 0 && proc_x < Env_nproc_x( env ) );
  assert( proc_y >= 0 && proc_y < Env_nproc_y( env ) );
*/

  const int nproc_x = Env_nproc_x( env );
  const int nproc_y = Env_nproc_y( env );
  const int nblock  = Sweeper_nblock( sweeper, env );

  const int octants_visited[NOCTANT] = { 0, 4, 1, 5, 3, 7, 2, 6 };

  Step_Info step_info;
  int octant_index = 0;
  int wave = 0;
  int step_base = 0;
  int block = 0;
  int octant = 0;
  int dir_x = 0;
  int dir_y = 0;
  int dir_z = 0;
  int start_x = 0;
  int start_y = 0;
  int start_z = 0;

  /*---First compute the octant number, in the order they are visited,
       and the wavefront number for that octant, starting from the
       beginning corner.
       Check every octant/wavefront in sequence to determine which
       one might be active for the proc in question.
  ---*/

  if ( Bool_true )
  {
    wave = step - ( step_base );
    octant_index = 0;
  }
  step_base += nblock;
  if ( step >= ( step_base + proc_x + proc_y ) )
  {
    wave = step - ( step_base );
    octant_index = 1;
  }
  step_base += nblock;
  if ( step >= ( step_base + proc_x + proc_y ) )
  {
    wave = step - ( step_base + (nproc_x-1) );
    octant_index = 2;
  }
  step_base += nblock + (nproc_x-1);
  if ( step >= ( step_base + (nproc_x-1-proc_x) + proc_y ) )
  {
    wave = step - ( step_base );
    octant_index = 3;
  }
  step_base += nblock;
  if ( step >= ( step_base + (nproc_x-1-proc_x) + proc_y ) )
  {
    wave = step - ( step_base + (nproc_y-1) );
    octant_index = 4;
  }
  step_base += nblock + (nproc_y-1);
  if ( step >= ( step_base + (nproc_x-1-proc_x)
                           + (nproc_y-1-proc_y) ) )
  {
    wave = step - ( step_base );
    octant_index = 5;
  }
  step_base += nblock;
  if ( step >= ( step_base + (nproc_x-1-proc_x)
                           + (nproc_y-1-proc_y) ) )
  {
    wave = step - ( step_base + (nproc_x-1) );
    octant_index = 6;
  }
  step_base += nblock + (nproc_x-1);
  if ( step >= ( step_base + proc_x + (nproc_y-1-proc_y) ) )
  {
    wave = step - ( step_base );
    octant_index = 7;
  }

  octant = octants_visited[octant_index];

  /*---Next convert the wavefront number to a block number based on
       location in the domain.  Use the equation that defines the plane.
  ---*/

  dir_x  = Dir_x( octant );
  dir_y  = Dir_y( octant );
  dir_z  = Dir_z( octant );

  /*---Get coordinates of the starting corner block of the wavefront---*/
  start_x = dir_x==Dir_up() ? 0 : ( nproc_x - 1 );
  start_y = dir_y==Dir_up() ? 0 : ( nproc_y - 1 );
  start_z = dir_z==Dir_up() ? 0 : ( nblock  - 1 );

  /*---Get coordinate of block on this processor to be processed---*/
  block = ( wave - ( start_x + proc_x * dir_x)
                 - ( start_y + proc_y * dir_y)
                 - ( start_z ) ) / dir_z;

  /*---Now determine whether the block calculation is active based on whether
       the block in question falls within the physical domain.
  ---*/

  step_info.is_active = block  >= 0 && block  < nblock &&
                        step   >= 0 && step   < Sweeper_nstep( sweeper, env ) &&
                        proc_x >= 0 && proc_x < Env_nproc_x( env ) &&
                        proc_y >= 0 && proc_y < Env_nproc_y( env );

  /*---Set remaining values---*/

  step_info.block_z = step_info.is_active ? block  : -1;
  step_info.octant  = step_info.is_active ? octant : -1;

  return step_info;
}

/*===========================================================================*/
/*---Communicate faces---*/

void Sweeper_communicate_faces__(
  Sweeper*           sweeper,
  int                step,
  Dimensions         dims_b,
  Env*               env )
{
  const int proc_x = Env_proc_x_this( env );
  const int proc_y = Env_proc_y_this( env );

  const size_t size_facexz = Dimensions_size_facexz( dims_b, NU,
                                        Sweeper_num_face_octants_allocated() );
  const size_t size_faceyz = Dimensions_size_faceyz( dims_b, NU,
                                        Sweeper_num_face_octants_allocated() );

  /*---Allocate temporary face buffers---*/

  P* __restrict__ buf_facexz  = malloc_P( size_facexz );
  P* __restrict__ buf_faceyz  = malloc_P( size_faceyz );

  /*---Communicate +/-X, +/-Y---*/

  int axis = 0;

  for( axis=0; axis<2; ++axis )
  {
    const Bool_t axis_x = axis==0;
    const Bool_t axis_y = axis==1;

    const int proc_axis = axis_x ? proc_x : proc_y;

    const size_t    size_face    = axis_x ? size_faceyz     : size_facexz;
    P* __restrict__ buf_face     = axis_x ? buf_faceyz      : buf_facexz;
    P* __restrict__ sweeper_face = axis_x ? sweeper->faceyz : sweeper->facexz;
    int dir_ind = 0;

    for( dir_ind=0; dir_ind<2; ++dir_ind )
    {
      const int dir = dir_ind==0 ? Dir_up() : Dir_dn();
      const int inc_x = axis_x ? Dir_inc( dir ) : 0;
      const int inc_y = axis_y ? Dir_inc( dir ) : 0;

      /*---Get step info for processors involved in communication---*/

      const Step_Info step_info_send_source =
        Sweeper_step_info__( sweeper, step,   proc_x,       proc_y,       env );

      const Step_Info step_info_send_target =
        Sweeper_step_info__( sweeper, step+1, proc_x+inc_x, proc_y+inc_y, env );

      const Step_Info step_info_recv_source =
        Sweeper_step_info__( sweeper, step,   proc_x-inc_x, proc_y-inc_y, env );

      const Step_Info step_info_recv_target =
        Sweeper_step_info__( sweeper, step+1, proc_x,       proc_y,       env );

      /*---Determine whether to communicate---*/

      Bool_t const do_send = step_info_send_source.is_active
                          && step_info_send_target.is_active
                          && step_info_send_source.octant ==
                             step_info_send_target.octant
                          && step_info_send_source.block_z ==
                             step_info_send_target.block_z
                          && ( axis_x ?
                               Dir_x( step_info_send_target.octant ) :
                               Dir_y( step_info_send_target.octant ) ) == dir;

      Bool_t const do_recv = step_info_recv_source.is_active
                          && step_info_recv_target.is_active
                          && step_info_recv_source.octant ==
                             step_info_recv_target.octant
                          && step_info_recv_source.block_z ==
                             step_info_recv_target.block_z
                          && ( axis_x ?
                               Dir_x( step_info_recv_target.octant ) :
                               Dir_y( step_info_recv_target.octant ) ) == dir;

      /*---Communicate as needed - use red/black coloring to avoid deadlock---*/

      int color = 0;

      Bool_t use_buf = Bool_false;

      for( color=0; color<2; ++color )
      {
        if( color == 0 )
        {
           if( proc_axis % 2 == 0 )
           {
             if( do_send )
             {
               const int proc_other
                                = Env_proc( env, proc_x+inc_x, proc_y+inc_y );
               Env_send_P( sweeper_face, size_face, proc_other, env->tag );
             }
           }
           else
           {
             if( do_recv )
             {
               const int proc_other
                                = Env_proc( env, proc_x-inc_x, proc_y-inc_y );
               /*---save copy else color 0 recv will destroy color 1 send---*/
               copy_vector( buf_face, sweeper_face, size_face );
               use_buf = Bool_true;
               Env_recv_P( sweeper_face, size_face, proc_other, env->tag );
             }
           }
        }
        else /*---if color---*/
        {
           if( proc_axis % 2 == 0 )
           {
             if( do_recv )
             {
               const int proc_other
                                = Env_proc( env, proc_x-inc_x, proc_y-inc_y );
               Env_recv_P( sweeper_face, size_face, proc_other, env->tag );
             }
           }
           else
           {
             if( do_send )
             {
               const int proc_other
                                = Env_proc( env, proc_x+inc_x, proc_y+inc_y );
               Env_send_P( use_buf ? buf_face : sweeper_face,
                                             size_face, proc_other, env->tag );
             }
           }
        } /*---if color---*/
      } /*---color---*/
    } /*---dir_ind---*/
  } /*---axis---*/

  /*---Deallocations---*/

  free_P( buf_facexz );
  free_P( buf_faceyz );

  /*---Increment message tag---*/

  env->tag++;
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

  int ix_b = 0;
  int iy_b = 0;
  int iz_b = 0;
  int ix = 0;
  int iy = 0;
  int iz = 0;
  int ie = 0;
  int im = 0;
  int ia = 0;
  int iu = 0;
  int step = -1;

  const int octant_ind = 0;
  const Dimensions dims_b = sweeper->dims_b;
  const int proc_x = Env_proc_x_this( env );
  const int proc_y = Env_proc_y_this( env );
  const int ix_base = quan->ix_base;
  const int iy_base = quan->iy_base;

  Dimensions dims_g = dims;
  dims_g.nx = quan->nx_g;
  dims_g.ny = quan->ny_g;

  /*---Initialize result array to zero---*/

  initialize_state_zero( vo, dims, NU );

  /*--------------------*/
  /*---Loop over kba parallel steps---*/
  /*--------------------*/

  for( step=0; step<Sweeper_nstep( sweeper, env ); ++step )
  {
    /*---Get step info for this proc---*/

    const Step_Info step_info = Sweeper_step_info__( sweeper, step,
                                                   proc_x, proc_y, env );

    if( step_info.is_active )
    {
      const int octant  = step_info.octant;
      const int block_z = step_info.block_z;

      const int iz_base = block_z * dims_b.nz;
      const int nz_b    = ( dims.nz - iz_base ) < dims_b.nz ?
                          ( dims.nz - iz_base ) : dims_b.nz;

      const int dir_x = Dir_x( octant );
      const int dir_y = Dir_y( octant );
      const int dir_z = Dir_z( octant );

      /*--------------------*/
      /*---Set XY face boundary conditions---*/
      /*--------------------*/

      if( ( dir_z == Dir_up() && block_z == 0 ) ||
          ( dir_z == Dir_dn() && block_z == sweeper->nblock_z-1 ) )
      {
        const int iz_g = dir_z == Dir_up() ? -1 : dims_g.nz;
        for( iu=0; iu<NU; ++iu )
        {
        for( iy_b=0; iy_b<dims_b.ny; ++iy_b )
        {
          const int iy_g = iy_b + iy_base;
        for( ix_b=0; ix_b<dims_b.nx; ++ix_b )
        {
          const int ix_g = ix_b + ix_base;
        for( ie=0; ie<dims.ne; ++ie )
        {
        for( ia=0; ia<dims.na; ++ia )
        {
          *ref_facexy( sweeper->facexy, dims_b, NU,
                       Sweeper_num_face_octants_allocated(),
                       ix_b, iy_b, ie, ia, iu, octant_ind )
              = Quantities_init_facexy(
                          quan, ix_g, iy_g, iz_g, ie, ia, iu, octant, dims_g );
        }
        }
        }
        }
        }
      } /*---dir_z---*/

      /*--------------------*/
      /*---Set XZ face boundary conditions---*/
      /*--------------------*/

      if( ( dir_y == Dir_up() && proc_y == 0 ) ||
          ( dir_y == Dir_dn() && proc_y == Env_nproc_y( env )-1 ) )
      {
        const int iy_g = dir_y == Dir_up() ? -1 : dims_g.ny;
        for( iu=0; iu<NU; ++iu )
        {
        for( iz_b=0; iz_b<nz_b; ++iz_b )
        {
          const int iz_g = iz_b + iz_base;
        for( ix_b=0; ix_b<dims_b.nx; ++ix_b )
        {
          const int ix_g = ix_b + ix_base;
        for( ie=0; ie<dims.ne; ++ie )
        {
        for( ia=0; ia<dims.na; ++ia )
        {
          *ref_facexz( sweeper->facexz, dims_b, NU,
                       Sweeper_num_face_octants_allocated(),
                       ix_b, iz_b, ie, ia, iu, octant_ind )
              = Quantities_init_facexz(
                           quan, ix_g, iy_g, iz_g, ie, ia, iu, octant, dims_g );
        }
        }
        }
        }
        }
      } /*---dir_y---*/

      /*--------------------*/
      /*---Set YZ face boundary conditions---*/
      /*--------------------*/

      if( ( dir_x == Dir_up() && proc_x == 0 ) ||
          ( dir_x == Dir_dn() && proc_x == Env_nproc_x( env )-1 ) )
      {
        const int ix_g = dir_x == Dir_up() ? -1 : dims_g.nx;
        for( iu=0; iu<NU; ++iu )
        {
        for( iz_b=0; iz_b<nz_b; ++iz_b )
        {
          const int iz_g = iz_b + iz_base;
        for( iy_b=0; iy_b<dims_b.ny; ++iy_b )
        {
          const int iy_g = iy_b + iy_base;
        for( ie=0; ie<dims.ne; ++ie )
        {
        for( ia=0; ia<dims.na; ++ia )
        {
          *ref_faceyz( sweeper->faceyz, dims_b, NU,
                       Sweeper_num_face_octants_allocated(),
                       iy_b, iz_b, ie, ia, iu, octant_ind )
              = Quantities_init_faceyz(
                          quan, ix_g, iy_g, iz_g, ie, ia, iu, octant, dims_g );
        }
        }
        }
        }
        }
      } /*---dir_x---*/

      /*--------------------*/
      /*---Loop over energy groups---*/
      /*--------------------*/

      for( ie=0; ie<dims.ne; ++ie )
      {

        /*---Calculate spatial loop extents---*/

        const int ixbeg = dir_x==Dir_up() ? 0       : dims_b.nx-1;
        const int iybeg = dir_y==Dir_up() ? 0       : dims_b.ny-1;
        const int izbeg = dir_z==Dir_up() ? iz_base : iz_base+nz_b-1;

        const int ixend = dir_x==Dir_dn() ? 0       : dims_b.nx-1;
        const int iyend = dir_y==Dir_dn() ? 0       : dims_b.ny-1;
        const int izend = dir_z==Dir_dn() ? iz_base : iz_base+nz_b-1;

        /*--------------------*/
        /*---Loop over gridcells, in proper direction---*/
        /*--------------------*/

        for( iz=izbeg; iz!=izend+Dir_inc(dir_z); iz+=Dir_inc(dir_z) )
        for( iy=iybeg; iy!=iyend+Dir_inc(dir_y); iy+=Dir_inc(dir_y) )
        for( ix=ixbeg; ix!=ixend+Dir_inc(dir_x); ix+=Dir_inc(dir_x) )
        {


          /*--------------------*/
          /*---Transform state vector from moments to angles---*/
          /*--------------------*/

          for( iu=0; iu<NU; ++iu )
          for( ia=0; ia<dims.na; ++ia )
          {
            P result = P_zero();
            for( im=0; im<dims.nm; ++im )
            {
              result +=
                  *const_ref_a_from_m( quan->a_from_m, dims, im, ia, octant ) *
                  *const_ref_state( vi, dims, NU, ix, iy, iz, ie, im, iu );
            }
            *ref_v_local( sweeper->v_local, dims, NU, ia, iu ) = result;
          }

          /*--------------------*/
          /*---Perform solve---*/
          /*--------------------*/

          Quantities_solve( quan, sweeper->v_local,
                            sweeper->facexy, sweeper->facexz, sweeper->faceyz,
                            ix, iy, iz-iz_base, ie, ix+ix_base, iy+iy_base, iz,
                            octant, octant_ind,
                            Sweeper_num_face_octants_allocated(),
                            dims_b, dims_g );

          /*--------------------*/
          /*---Transform state vector from angles to moments---*/
          /*--------------------*/

          for( iu=0; iu<NU; ++iu )
          for( im=0; im<dims.nm; ++im )
          {
            P result = P_zero();
            for( ia=0; ia<dims.na; ++ia )
            {
              result +=
                  *const_ref_m_from_a( quan->m_from_a, dims, im, ia, octant ) *
                  *const_ref_v_local( sweeper->v_local, dims, NU, ia, iu );
            }
            *ref_state( vo, dims, NU, ix, iy, iz, ie, im, iu ) += result;
          }

        } /*---ix/iy/iz---*/

      } /*---ie---*/

    }  /*---is_active---*/

    /*--------------------*/
    /*---Communicate faces---*/
    /*--------------------*/

    Sweeper_communicate_faces__( sweeper, step, dims_b, env );

  } /*---step---*/

} /*---sweep---*/

/*===========================================================================*/

#endif /*---_serial_c__sweeper_kba_c_h_---*/

/*---------------------------------------------------------------------------*/
