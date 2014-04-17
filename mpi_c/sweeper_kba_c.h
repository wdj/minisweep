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
#include "step_scheduler_kba.h"
#include "sweeper_kba.h"


/*===========================================================================*/
/*---Pseudo-constructor for Sweeper struct---*/

void Sweeper_ctor( Sweeper*          sweeper,
                   Dimensions        dims,
                   const Quantities* quan,
                   Env*              env,
                   Arguments*        args )
{
  Insist( dims.nx > 0 && "KBA sweeper currently requires all blocks nonempty" );
  Insist( dims.ny > 0 && "KBA sweeper currently requires all blocks nonempty" );
  Insist( dims.nz > 0 && "KBA sweeper currently requires all blocks nonempty" );

  /*---Set up number of kba blocks---*/
  sweeper->nblock_z = Arguments_consume_int_or_default( args, "--nblock_z", 1);
  Insist( sweeper->nblock_z > 0 && "Invalid z blocking factor supplied" );
  Insist( dims.nz % sweeper->nblock_z == 0 &&
          "KBA sweeper currently requires all blocks have same z dimension" );

  /*---Set up number of octant threads---*/
  sweeper->nthread_octant
              = Arguments_consume_int_or_default( args, "--nthread_octant", 1);
  /*---Require a power of 2 between 1 and 8 inclusive---*/
  Insist( sweeper->nthread_octant>0 && sweeper->nthread_octant<=NOCTANT
          && ((sweeper->nthread_octant&(sweeper->nthread_octant-1))==0)
          && "Invalid octant thread count supplied" );
  sweeper->noctant_per_block = sweeper->nthread_octant;
  sweeper->nblock_octant = NOCTANT / sweeper->noctant_per_block;

  /*---Set up number of semiblock steps---*/
  sweeper->nsemiblock = Arguments_consume_int_or_default(
                               args, "--nsemiblock", sweeper->nthread_octant );
  Insist( sweeper->nsemiblock>0 && sweeper->nsemiblock<=NOCTANT
          && ((sweeper->nsemiblock&(sweeper->nsemiblock-1))==0)
          && "Invalid semiblock count supplied" );
  Insist( ( sweeper->nsemiblock >= sweeper->nthread_octant ||
            IS_USING_OPENMP_VO_ATOMIC ) &&
          "Incomplete set of semiblock steps requires atomic vo update" );

  /*---Set up number of energy threads---*/
  sweeper->nthread_e
                   = Arguments_consume_int_or_default( args, "--nthread_e", 1);
  Insist( sweeper->nthread_e > 0 && "Invalid e thread count supplied." );

  /*---Set up step scheduler---*/
  Step_Scheduler_ctor( &(sweeper->step_scheduler),
                              sweeper->nblock_z, sweeper->nblock_octant, env );

  /*---Set up dims structs---*/
  sweeper->dims = dims;

  sweeper->dims_b = sweeper->dims;
  sweeper->dims_b.nz = sweeper->dims.nz / sweeper->nblock_z;

  sweeper->dims_g = sweeper->dims;
  sweeper->dims_g.nx = quan->nx_g;
  sweeper->dims_g.ny = quan->ny_g;

  /*---Allocate arrays---*/

  sweeper->v_local = malloc_P( sweeper->dims_b.na * NU *
                               sweeper->nthread_e * sweeper->nthread_octant );

  sweeper->facexy0 = malloc_P( Dimensions_size_facexy( sweeper->dims_b, NU,
                                                sweeper->noctant_per_block ) );
  sweeper->facexz0 = malloc_P( Dimensions_size_facexz( sweeper->dims_b, NU,
                                                sweeper->noctant_per_block ) );
  sweeper->faceyz0 = malloc_P( Dimensions_size_faceyz( sweeper->dims_b, NU,
                                                sweeper->noctant_per_block ) );

  sweeper->facexz1 = NULL;
  sweeper->facexz2 = NULL;
  sweeper->faceyz1 = NULL;
  sweeper->faceyz2 = NULL;

  if( Sweeper_is_face_comm_async() )
  {
    sweeper->facexz1 = malloc_P( Dimensions_size_facexz( sweeper->dims_b, NU,
                                                sweeper->noctant_per_block ) );
    sweeper->facexz2 = malloc_P( Dimensions_size_facexz( sweeper->dims_b, NU,
                                                sweeper->noctant_per_block ) );
    sweeper->faceyz1 = malloc_P( Dimensions_size_faceyz( sweeper->dims_b, NU,
                                                sweeper->noctant_per_block ) );
    sweeper->faceyz2 = malloc_P( Dimensions_size_faceyz( sweeper->dims_b, NU,
                                                sweeper->noctant_per_block ) );
  }
}

/*===========================================================================*/
/*---Pseudo-destructor for Sweeper struct---*/

void Sweeper_dtor( Sweeper* sweeper )
{
  /*---Deallocate arrays---*/

  free_P( sweeper->v_local );
  sweeper->v_local = NULL;

  free_P( sweeper->facexy0 );
  free_P( sweeper->facexz0 );
  free_P( sweeper->faceyz0 );
  sweeper->facexy0 = NULL;
  sweeper->facexz0 = NULL;
  sweeper->faceyz0 = NULL;

  if( Sweeper_is_face_comm_async() )
  {
    free_P( sweeper->facexz1 );
    free_P( sweeper->facexz2 );
    free_P( sweeper->faceyz1 );
    free_P( sweeper->faceyz2 );
    sweeper->facexz1 = NULL;
    sweeper->facexz2 = NULL;
    sweeper->faceyz1 = NULL;
    sweeper->faceyz2 = NULL;
  }

  Step_Scheduler_dtor( &( sweeper->step_scheduler ) );
}

/*===========================================================================*/
/*---Determine whether to send a face computed at step, used at step+1---*/

Bool_t Sweeper_must_do_send__(
  Sweeper*           sweeper,
  int                step,
  int                axis,
  int                dir_ind,
  int                octant_in_block,
  Env*               env )
{
  const int proc_x = Env_proc_x_this( env );
  const int proc_y = Env_proc_y_this( env );

  const Bool_t axis_x = axis==0;
  const Bool_t axis_y = axis==1;

  const int dir = dir_ind==0 ? Dir_up() : Dir_dn();
  const int inc_x = axis_x ? Dir_inc( dir ) : 0;
  const int inc_y = axis_y ? Dir_inc( dir ) : 0;

  /*---Get step info for processors involved in communication---*/

  const Step_Info step_info_send_source_step = Step_Scheduler_step_info(
    &(sweeper->step_scheduler), step,   octant_in_block,
                                                  proc_x,       proc_y);

  const Step_Info step_info_send_target_step = Step_Scheduler_step_info(
    &(sweeper->step_scheduler), step+1, octant_in_block,
                                                  proc_x+inc_x, proc_y+inc_y );

  /*---Determine whether to communicate---*/

  Bool_t const do_send = step_info_send_source_step.is_active
                      && step_info_send_target_step.is_active
                      && step_info_send_source_step.octant ==
                         step_info_send_target_step.octant
                      && step_info_send_source_step.block_z ==
                         step_info_send_target_step.block_z
                      && ( axis_x ?
                           Dir_x( step_info_send_target_step.octant ) :
                           Dir_y( step_info_send_target_step.octant ) ) == dir;

  return do_send;
}

/*===========================================================================*/
/*---Determine whether to recv a face computed at step, used at step+1---*/

Bool_t Sweeper_must_do_recv__(
  Sweeper*           sweeper,
  int                step,
  int                axis,
  int                dir_ind,
  int                octant_in_block,
  Env*               env )
{
  const int proc_x = Env_proc_x_this( env );
  const int proc_y = Env_proc_y_this( env );

  const Bool_t axis_x = axis==0;
  const Bool_t axis_y = axis==1;

  const int dir = dir_ind==0 ? Dir_up() : Dir_dn();
  const int inc_x = axis_x ? Dir_inc( dir ) : 0;
  const int inc_y = axis_y ? Dir_inc( dir ) : 0;

  /*---Get step info for processors involved in communication---*/

  const Step_Info step_info_recv_source_step = Step_Scheduler_step_info(
    &(sweeper->step_scheduler), step,   octant_in_block,
                                                  proc_x-inc_x, proc_y-inc_y );

  const Step_Info step_info_recv_target_step = Step_Scheduler_step_info(
    &(sweeper->step_scheduler), step+1, octant_in_block,
                                                  proc_x,       proc_y );

  /*---Determine whether to communicate---*/

  Bool_t const do_recv = step_info_recv_source_step.is_active
                      && step_info_recv_target_step.is_active
                      && step_info_recv_source_step.octant ==
                         step_info_recv_target_step.octant
                      && step_info_recv_source_step.block_z ==
                         step_info_recv_target_step.block_z
                      && ( axis_x ?
                           Dir_x( step_info_recv_target_step.octant ) :
                           Dir_y( step_info_recv_target_step.octant ) ) == dir;

  return do_recv;
}

/*===========================================================================*/
/*---Communicate faces computed at step, used at step+1---*/

void Sweeper_communicate_faces__(
  Sweeper*           sweeper,
  int                step,
  Env*               env )
{
  assert( ! Sweeper_is_face_comm_async() );

  const int proc_x = Env_proc_x_this( env );
  const int proc_y = Env_proc_y_this( env );

  const size_t size_facexz_per_octant = Dimensions_size_facexz( sweeper->dims_b,
                 NU, sweeper->noctant_per_block ) / sweeper->noctant_per_block;
  const size_t size_faceyz_per_octant = Dimensions_size_faceyz( sweeper->dims_b,
                 NU, sweeper->noctant_per_block ) / sweeper->noctant_per_block;

  /*---Allocate temporary face buffers---*/

  P* __restrict__ buf_xz  = malloc_P( size_facexz_per_octant );
  P* __restrict__ buf_yz  = malloc_P( size_faceyz_per_octant );

  /*---Loop over octants---*/

  int octant_in_block = 0;

  for( octant_in_block=0; octant_in_block<sweeper->noctant_per_block;
                                                            ++octant_in_block )
  {
    /*---Communicate +/-X, +/-Y---*/

    int axis = 0;

    for( axis=0; axis<2; ++axis )  /*---Loop: X, Y---*/
    {
      const Bool_t axis_x = axis==0;
      const Bool_t axis_y = axis==1;

      const int proc_axis = axis_x ? proc_x : proc_y;

      const size_t    size_face_per_octant    = axis_x ? size_faceyz_per_octant
                                                       : size_facexz_per_octant;
      P* __restrict__ buf                     = axis_x ? buf_yz
                                                       : buf_xz;
      P* __restrict__ face_per_octant = axis_x ?
        ref_faceyz( Sweeper_faceyz__( sweeper, step ), sweeper->dims_b, NU,
                  sweeper->noctant_per_block, 0, 0, 0, 0, 0, octant_in_block ) :
        ref_facexz( Sweeper_facexz__( sweeper, step ), sweeper->dims_b, NU,
                  sweeper->noctant_per_block, 0, 0, 0, 0, 0, octant_in_block );

      int dir_ind = 0;

      for( dir_ind=0; dir_ind<2; ++dir_ind ) /*---Loop: up, down---*/
      {
        const int dir = dir_ind==0 ? Dir_up() : Dir_dn();
        const int inc_x = axis_x ? Dir_inc( dir ) : 0;
        const int inc_y = axis_y ? Dir_inc( dir ) : 0;

        /*---Determine whether to communicate---*/

        Bool_t const do_send = Sweeper_must_do_send__(
                          sweeper, step, axis, dir_ind, octant_in_block, env );

        Bool_t const do_recv = Sweeper_must_do_recv__(
                          sweeper, step, axis, dir_ind, octant_in_block, env );

        /*---Communicate as needed - red/black coloring to avoid deadlock---*/

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
                Env_send_P( face_per_octant, size_face_per_octant,
                            proc_other, Env_tag( env )+octant_in_block );
              }
            }
            else
            {
              if( do_recv )
              {
                const int proc_other
                                 = Env_proc( env, proc_x-inc_x, proc_y-inc_y );
                /*---save copy else color 0 recv will destroy color 1 send---*/
                copy_vector( buf, face_per_octant, size_face_per_octant );
                use_buf = Bool_true;
                Env_recv_P( face_per_octant, size_face_per_octant,
                            proc_other, Env_tag( env )+octant_in_block );
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
                Env_recv_P( face_per_octant, size_face_per_octant,
                            proc_other, Env_tag( env )+octant_in_block );
              }
            }
            else
            {
              if( do_send )
              {
                const int proc_other
                                 = Env_proc( env, proc_x+inc_x, proc_y+inc_y );
                Env_send_P( use_buf ? buf : face_per_octant,
                  size_face_per_octant, proc_other, Env_tag( env )+octant_in_block );
              }
            }
          } /*---if color---*/
        } /*---color---*/
      } /*---dir_ind---*/
    } /*---axis---*/
  } /*---octant_in_block---*/

  /*---Deallocations---*/

  free_P( buf_xz );
  free_P( buf_yz );
}

/*===========================================================================*/
/*---Asynchronously send faces computed at step, used at step+1: start---*/

void Sweeper_send_faces_start__(
  Sweeper*           sweeper,
  int                step,
  Env*               env )
{
  assert( Sweeper_is_face_comm_async() );

  const int proc_x = Env_proc_x_this( env );
  const int proc_y = Env_proc_y_this( env );

  const size_t size_facexz_per_octant = Dimensions_size_facexz( sweeper->dims_b,
                 NU, sweeper->noctant_per_block ) / sweeper->noctant_per_block;
  const size_t size_faceyz_per_octant = Dimensions_size_faceyz( sweeper->dims_b,
                 NU, sweeper->noctant_per_block ) / sweeper->noctant_per_block;

  /*---Loop over octants---*/

  int octant_in_block = 0;

  for( octant_in_block=0; octant_in_block<sweeper->noctant_per_block;
                                                            ++octant_in_block )
  {
    /*---Communicate +/-X, +/-Y---*/

    int axis = 0;

    for( axis=0; axis<2; ++axis )
    {
      const Bool_t axis_x = axis==0;
      const Bool_t axis_y = axis==1;

      const int proc_axis = axis_x ? proc_x : proc_y;

      /*---Send values computed on this step---*/

      const size_t    size_face_per_octant    = axis_x ? size_faceyz_per_octant
                                                       : size_facexz_per_octant;
      P* __restrict__ face_per_octant = axis_x ?
        ref_faceyz( Sweeper_faceyz__( sweeper, step ), sweeper->dims_b, NU,
                  sweeper->noctant_per_block, 0, 0, 0, 0, 0, octant_in_block ) :
        ref_facexz( Sweeper_facexz__( sweeper, step ), sweeper->dims_b, NU,
                  sweeper->noctant_per_block, 0, 0, 0, 0, 0, octant_in_block );

      int dir_ind = 0;

      for( dir_ind=0; dir_ind<2; ++dir_ind )
      {
        const int dir = dir_ind==0 ? Dir_up() : Dir_dn();
        const int inc_x = axis_x ? Dir_inc( dir ) : 0;
        const int inc_y = axis_y ? Dir_inc( dir ) : 0;

        /*---Determine whether to communicate---*/

        Bool_t const do_send = Sweeper_must_do_send__(
                          sweeper, step, axis, dir_ind, octant_in_block, env );

        if( do_send )
        {
          const int proc_other = Env_proc( env, proc_x+inc_x, proc_y+inc_y );
          Request_t* request = axis_x ?
                                   & sweeper->request_send_xz[octant_in_block]
                                 : & sweeper->request_send_yz[octant_in_block];
          Env_asend_P( face_per_octant, size_face_per_octant,
                       proc_other, Env_tag( env )+octant_in_block, request );
        }
      } /*---dir_ind---*/
    } /*---axis---*/
  } /*---octant_in_block---*/
}

/*===========================================================================*/
/*---Asynchronously send faces computed at step, used at step+1: end---*/

void Sweeper_send_faces_end__(
  Sweeper*           sweeper,
  int                step,
  Env*               env )
{
  assert( Sweeper_is_face_comm_async() );

  const int proc_x = Env_proc_x_this( env );
  const int proc_y = Env_proc_y_this( env );

  /*---Loop over octants---*/

  int octant_in_block = 0;

  for( octant_in_block=0; octant_in_block<sweeper->noctant_per_block;
                                                            ++octant_in_block )
  {
    /*---Communicate +/-X, +/-Y---*/

    int axis = 0;

    for( axis=0; axis<2; ++axis )
    {
      const Bool_t axis_x = axis==0;
      const Bool_t axis_y = axis==1;

      int dir_ind = 0;

      for( dir_ind=0; dir_ind<2; ++dir_ind )
      {

        /*---Determine whether to communicate---*/

        Bool_t const do_send = Sweeper_must_do_send__(
                          sweeper, step, axis, dir_ind, octant_in_block, env );

        if( do_send )
        {
          Request_t* request = axis_x ?
                                   & sweeper->request_send_xz[octant_in_block]
                                 : & sweeper->request_send_yz[octant_in_block];
          Env_wait( request );
        }
      } /*---dir_ind---*/
    } /*---axis---*/
  } /*---octant_in_block---*/
}

/*===========================================================================*/
/*---Asynchronously recv faces computed at step, used at step+1: start---*/

void Sweeper_recv_faces_start__(
  Sweeper*           sweeper,
  int                step,
  Env*               env )
{
  assert( Sweeper_is_face_comm_async() );

  const int proc_x = Env_proc_x_this( env );
  const int proc_y = Env_proc_y_this( env );

  const size_t size_facexz_per_octant = Dimensions_size_facexz( sweeper->dims_b,
                 NU, sweeper->noctant_per_block ) / sweeper->noctant_per_block;
  const size_t size_faceyz_per_octant = Dimensions_size_faceyz( sweeper->dims_b,
                 NU, sweeper->noctant_per_block ) / sweeper->noctant_per_block;

  /*---Loop over octants---*/

  int octant_in_block = 0;

  for( octant_in_block=0; octant_in_block<sweeper->noctant_per_block;
                                                            ++octant_in_block )
  {
    /*---Communicate +/-X, +/-Y---*/

    int axis = 0;

    for( axis=0; axis<2; ++axis )
    {
      const Bool_t axis_x = axis==0;
      const Bool_t axis_y = axis==1;

      const int proc_axis = axis_x ? proc_x : proc_y;

      /*---Receive values computed on the next step---*/

      const size_t    size_face_per_octant    = axis_x ? size_faceyz_per_octant
                                                       : size_facexz_per_octant;
      P* __restrict__ face_per_octant = axis_x ?
        ref_faceyz( Sweeper_faceyz__( sweeper, step+1 ), sweeper->dims_b, NU,
                  sweeper->noctant_per_block, 0, 0, 0, 0, 0, octant_in_block ) :
        ref_facexz( Sweeper_facexz__( sweeper, step+1 ), sweeper->dims_b, NU,
                  sweeper->noctant_per_block, 0, 0, 0, 0, 0, octant_in_block );

      int dir_ind = 0;

      for( dir_ind=0; dir_ind<2; ++dir_ind )
      {
        const int dir = dir_ind==0 ? Dir_up() : Dir_dn();
        const int inc_x = axis_x ? Dir_inc( dir ) : 0;
        const int inc_y = axis_y ? Dir_inc( dir ) : 0;

        /*---Determine whether to communicate---*/

        Bool_t const do_recv = Sweeper_must_do_recv__(
                          sweeper, step, axis, dir_ind, octant_in_block, env );

        if( do_recv )
        {
          const int proc_other = Env_proc( env, proc_x-inc_x, proc_y-inc_y );
          Request_t* request = axis_x ?
                                   & sweeper->request_recv_xz[octant_in_block]
                                 : & sweeper->request_recv_yz[octant_in_block];
          Env_arecv_P( face_per_octant, size_face_per_octant,
                       proc_other, Env_tag( env )+octant_in_block, request );
        }
      } /*---dir_ind---*/
    } /*---axis---*/
  } /*---octant_in_block---*/
}

/*===========================================================================*/
/*---Asynchronously recv faces computed at step, used at step+1: end---*/

void Sweeper_recv_faces_end__(
  Sweeper*           sweeper,
  int                step,
  Env*               env )
{
  assert( Sweeper_is_face_comm_async() );

  const int proc_x = Env_proc_x_this( env );
  const int proc_y = Env_proc_y_this( env );

  const size_t size_facexz_per_octant = Dimensions_size_facexz( sweeper->dims_b,
                 NU, sweeper->noctant_per_block ) / sweeper->noctant_per_block;
  const size_t size_faceyz_per_octant = Dimensions_size_faceyz( sweeper->dims_b,
                 NU, sweeper->noctant_per_block ) / sweeper->noctant_per_block;

  /*---Loop over octants---*/

  int octant_in_block = 0;

  for( octant_in_block=0; octant_in_block<sweeper->noctant_per_block;
                                                            ++octant_in_block )
  {
    /*---Communicate +/-X, +/-Y---*/

    int axis = 0;

    for( axis=0; axis<2; ++axis )
    {
      const Bool_t axis_x = axis==0;
      const Bool_t axis_y = axis==1;

      int dir_ind = 0;

      for( dir_ind=0; dir_ind<2; ++dir_ind )
      {
        /*---Determine whether to communicate---*/

        Bool_t const do_recv = Sweeper_must_do_recv__(
                          sweeper, step, axis, dir_ind, octant_in_block, env );

        if( do_recv )
        {
          Request_t* request = axis_x ?
                                   & sweeper->request_recv_xz[octant_in_block]
                                 : & sweeper->request_recv_yz[octant_in_block];
          Env_wait( request );
        }
      } /*---dir_ind---*/
    } /*---axis---*/
  } /*---octant_in_block---*/
}

/*===========================================================================*/
/*---Apply boundary condition: xy face---*/

static void Sweeper_set_boundary_xy(
  const Sweeper*        sweeper,
  P* const __restrict__ facexy,
  const Quantities*     quan,
  int                   octant,
  int                   octant_in_block, 
  const int             ixmin_b,
  const int             ixmax_b,
  const int             iymin_b,
  const int             iymax_b )
{
  const int ix_base = quan->ix_base;
  const int iy_base = quan->iy_base;
  const int dir_z = Dir_z( octant );
  const int iz_g = dir_z == Dir_up() ? -1 : sweeper->dims_g.nz;

  int ie = 0;

#ifdef USE_OPENMP_E
#pragma omp parallel for num_threads( sweeper->nthread_e )
#endif
  for( ie=0; ie<sweeper->dims_b.ne; ++ie )
  {
    int ix_b = 0;
    int iy_b = 0;
    int ia   = 0;
    int iu   = 0;
  for( iu=0; iu<NU; ++iu )
  {
  for( iy_b=iymin_b; iy_b<=iymax_b; ++iy_b )
  {
    const int iy_g = iy_b + iy_base;
  for( ix_b=ixmin_b; ix_b<=ixmax_b; ++ix_b )
  {
    const int ix_g = ix_b + ix_base;
  for( ia=0; ia<sweeper->dims_b.na; ++ia )
  {
    *ref_facexy( facexy, sweeper->dims_b, NU,
                 sweeper->noctant_per_block,
                 ix_b, iy_b, ie, ia, iu, octant_in_block )
        = Quantities_init_facexy(
                 quan, ix_g, iy_g, iz_g, ie, ia, iu, octant, sweeper->dims_g );
  }
  }
  }
  }
  }
}

/*===========================================================================*/
/*---Apply boundary condition: xz face---*/

static void Sweeper_set_boundary_xz(
  const Sweeper*        sweeper,
  P* const __restrict__ facexz,
  const Quantities*     quan,
  int                   block_z, 
  int                   octant,
  int                   octant_in_block, 
  const int             ixmin_b,
  const int             ixmax_b,
  const int             izmin_b,
  const int             izmax_b )
{
  const int ix_base = quan->ix_base;
  const int iz_base = block_z * sweeper->dims_b.nz;
  const int dir_y = Dir_y( octant );
  const int iy_g = dir_y == Dir_up() ? -1 : sweeper->dims_g.ny;

  int ie = 0;

#ifdef USE_OPENMP_E
#pragma omp parallel for num_threads( sweeper->nthread_e )
#endif
  for( ie=0; ie<sweeper->dims_b.ne; ++ie )
  {
    int ix_b = 0;
    int iz_b = 0;
    int ia   = 0;
    int iu   = 0;
  for( iu=0; iu<NU; ++iu )
  {
  for( iz_b=izmin_b; iz_b<=izmax_b; ++iz_b )
  {
    const int iz_g = iz_b + iz_base;
  for( ix_b=ixmin_b; ix_b<=ixmax_b; ++ix_b )
  {
    const int ix_g = ix_b + ix_base;
  for( ia=0; ia<sweeper->dims_b.na; ++ia )
  {
    *ref_facexz( facexz, sweeper->dims_b, NU,
                 sweeper->noctant_per_block,
                 ix_b, iz_b, ie, ia, iu, octant_in_block )
        = Quantities_init_facexz(
                 quan, ix_g, iy_g, iz_g, ie, ia, iu, octant, sweeper->dims_g );
  }
  }
  }
  }
  }
}

/*===========================================================================*/
/*---Apply boundary condition: yz face---*/

static void Sweeper_set_boundary_yz(
  const Sweeper*        sweeper,
  P* const __restrict__ faceyz,
  const Quantities*     quan,
  int                   block_z, 
  int                   octant,
  int                   octant_in_block, 
  const int             iymin_b,
  const int             iymax_b,
  const int             izmin_b,
  const int             izmax_b )
{
  const int iy_base = quan->iy_base;
  const int iz_base = block_z * sweeper->dims_b.nz;
  const int dir_x = Dir_x( octant );
  const int ix_g = dir_x == Dir_up() ? -1 : sweeper->dims_g.nx;

  int ie = 0;

#ifdef USE_OPENMP_E
#pragma omp parallel for num_threads( sweeper->nthread_e )
#endif
  for( ie=0; ie<sweeper->dims_b.ne; ++ie )
  {
    int iy_b = 0;
    int iz_b = 0;
    int ia   = 0;
    int iu   = 0;
  for( iu=0; iu<NU; ++iu )
  {
  for( iz_b=izmin_b; iz_b<=izmax_b; ++iz_b )
  {
    const int iz_g = iz_b + iz_base;
  for( iy_b=iymin_b; iy_b<=iymax_b; ++iy_b )
  {
    const int iy_g = iy_b + iy_base;
  for( ia=0; ia<sweeper->dims_b.na; ++ia )
  {
    *ref_faceyz( faceyz, sweeper->dims_b, NU,
                 sweeper->noctant_per_block,
                 iy_b, iz_b, ie, ia, iu, octant_in_block )
        = Quantities_init_faceyz(
                 quan, ix_g, iy_g, iz_g, ie, ia, iu, octant, sweeper->dims_g );
  }
  }
  }
  }
  }
}

/*===========================================================================*/
/*---Perform a sweep for a block---*/

void Sweeper_sweep_block(
  Sweeper*               sweeper,
  P* __restrict__        vo,
  const P* __restrict__  vi,
  P* __restrict__        facexy,
  P* __restrict__        facexz,
  P* __restrict__        faceyz,
  const Quantities*      quan,
  Env*                   env, 
  const Step_Info        step_info,
  const int              thread_num,
  const int              num_threads,
  const int              octant_in_block,
  const int              ixmin,
  const int              ixmax,
  const int              iymin,
  const int              iymax,
  const int              izmin,
  const int              izmax )
{
  /*---Declarations---*/

  const int octant  = step_info.octant;
  const int block_z = step_info.block_z;
  const int iz_base = block_z * sweeper->dims_b.nz;

  const int dir_x = Dir_x( octant );
  const int dir_y = Dir_y( octant );
  const int dir_z = Dir_z( octant );

  /*---Calculate spatial loop extents---*/

  const int ixbeg = dir_x==Dir_up() ? ixmin : ixmax;
  const int iybeg = dir_y==Dir_up() ? iymin : iymax;
  const int izbeg = dir_z==Dir_up() ? izmin : izmax;

  const int ixend = dir_x==Dir_dn() ? ixmin : ixmax;
  const int iyend = dir_y==Dir_dn() ? iymin : iymax;
  const int izend = dir_z==Dir_dn() ? izmin : izmax;

  int ie = 0;

#ifdef USE_OPENMP_E
#pragma omp parallel for num_threads( sweeper->nthread_e )
#endif
  for( ie=0; ie<sweeper->dims.ne; ++ie )
  {
    /*---Compute thread information---*/

    const int thread_num_outer  = thread_num;
    const int num_threads_outer = num_threads;
    const int thread_num  = thread_num_outer + num_threads_outer *
                  ( IS_USING_OPENMP_E ? Env_thread_this( env ) : 0 );
    const int num_threads =  num_threads_outer *
                  ( IS_USING_OPENMP_E ? Env_num_threads( env ) : 1 );

    /*---Get v_local part to be used by this thread---*/

    P* __restrict__ v_local = Sweeper_v_local_this__( sweeper, thread_num );

    int ix = 0;
    int iy = 0;
    int iz = 0;

    /*--------------------*/
    /*---Loop over gridcells, in proper direction---*/
    /*--------------------*/

    for( iz=izbeg; iz!=izend+Dir_inc(dir_z); iz+=Dir_inc(dir_z) )
    for( iy=iybeg; iy!=iyend+Dir_inc(dir_y); iy+=Dir_inc(dir_y) )
    for( ix=ixbeg; ix!=ixend+Dir_inc(dir_x); ix+=Dir_inc(dir_x) )
    {
      int im = 0;
      int ia = 0;
      int iu = 0;

      /*--------------------*/
      /*---Transform state vector from moments to angles---*/
      /*--------------------*/

      for( iu=0; iu<NU; ++iu )
      {
      for( ia=0; ia<sweeper->dims.na; ++ia )
      {
        P result = P_zero();
        for( im=0; im<sweeper->dims.nm; ++im )
        {
          result +=
            *const_ref_a_from_m( quan->a_from_m, sweeper->dims, im, ia, octant )
            * *const_ref_state( vi, sweeper->dims, NU, ix, iy, iz, ie, im, iu );
        }
        *ref_v_local( v_local, sweeper->dims, NU, ia, iu ) = result;
      }
      }

      /*--------------------*/
      /*---Perform solve---*/
      /*--------------------*/

      Quantities_solve( quan, v_local,
                        facexy, facexz, faceyz,
                        ix, iy, iz-iz_base, ie,
                        ix+quan->ix_base, iy+quan->iy_base, iz,
                        octant, octant_in_block,
                        sweeper->noctant_per_block,
                        sweeper->dims_b, sweeper->dims_g );

      /*--------------------*/
      /*---Transform state vector from angles to moments---*/
      /*--------------------*/

      for( iu=0; iu<NU; ++iu )
      {
      for( im=0; im<sweeper->dims.nm; ++im )
      {
        P result = P_zero();
        for( ia=0; ia<sweeper->dims.na; ++ia )
        {
          result +=
            *const_ref_m_from_a( quan->m_from_a, sweeper->dims, im, ia, octant )
            * *const_ref_v_local( v_local, sweeper->dims, NU, ia, iu );
        }
#ifdef USE_OPENMP_VO_ATOMIC
#pragma omp atomic
#endif
        *ref_state( vo, sweeper->dims, NU, ix, iy, iz, ie, im, iu ) += result;
      }
      }

    } /*---ix/iy/iz---*/

  } /*---ie---*/
}

/*===========================================================================*/
/*---Perform a sweep---*/

void Sweeper_sweep(
  Sweeper*               sweeper,
  P* __restrict__        vo,
  const P* __restrict__  vi,
  const Quantities*      quan,
  Env*                   env )
{
  assert( sweeper );
  assert( vi );
  assert( vo );

  /*---Declarations---*/

  const int proc_x = Env_proc_x_this( env );
  const int proc_y = Env_proc_y_this( env );

  int step = -1;

  /*---Initialize result array to zero---*/

  initialize_state_zero( vo, sweeper->dims, NU );

  /*--------------------*/
  /*---Loop over kba parallel steps---*/
  /*--------------------*/

  for( step=0; step<Step_Scheduler_nstep( &(sweeper->step_scheduler) ); ++step )
  {
    int semiblock = 0;

    /*---Pick up needed face pointers---*/

    /*=========================================================================
    =    Order is important here.
    =    The _r face for a step must match the _c face for the next step.
    =    The _s face for a step must match the _c face for the prev step.
    =========================================================================*/

    P* const __restrict__ facexy = Sweeper_facexy__( sweeper, step );
    P* const __restrict__ facexz = Sweeper_facexz__( sweeper, step );
    P* const __restrict__ faceyz = Sweeper_faceyz__( sweeper, step );

    /*--------------------*/
    /*---Communicate faces---*/
    /*--------------------*/

    /*=========================================================================
    =    Faces are triple buffered via a circular buffer of face arrays.
    =    The following shows the pattern of face usage over a step:
    =
    =                         step:     ...    i    i+1   i+2   i+3   ...
    =    ------------------------------------------------------------------
    =    Recv face for this step wait   ...  face0 face1 face2 face0  ...
    =    Recv face for next step start  ...  face1 face2 face0 face1  ...
    =    Compute this step using face   ...  face0 face1 face2 face0  ...
    =    Send face from last step wait  ...  face2 face0 face1 face2  ...
    =    Send face from this step start ...  face0 face1 face2 face0  ...
    =========================================================================*/

    if( Sweeper_is_face_comm_async() )
    {
      Sweeper_recv_faces_end__  (  sweeper, step-1, env );
      Sweeper_recv_faces_start__(  sweeper, step,   env );
    }

    /*=========================================================================
    =    OpenMP-parallelizing octants leads to the problem that for the same
    =    step, two octants may be updating the same location in a state vector.
    =
    =    One solution is to make the state vector update atomic, which is
    =    likely to be inefficient depending on the system.
    =
    =    The alternative used here is to break the step into sub steps
    =    and break the block into subregions such that during a sub-step,
    =    different octants in different threads update disjoint subregions.
    =
    =    First, note that octants are assigned to threads as follows:
    =      nthread_octant==1: one thread for all octants.
    =      nthread_octant==2: -x and +x octants assigned to different threads.
    =      nthread_octant==4: -y and +y octants also have different threads.
    =      nthread_octant==8: -z and +z octants also have different threads.
    =
    =    Along each coordinate axis for which two threads are assigned,
    =    the block is divided into two halves.  This gives a set of semiblocks.
    =
    =    The semiblocks are visited by the semiblock loop in an ordering
    =    which is lexicographical, either forward or reverse direction
    =    depending on the direction specified by that octant along the axis.
    =    This is set up so that (1) the disjointness condition described
    =    above holds, and (2) the cells are visited in an order that
    =    satisfies the sweep recursion.
    =
    =    NOTES:
    =    - For the unthreaded case, nsemiblock and noctant_per_block
    =      can be set to any of the allowed values and the algorithm will
    =      work properly.
    =    - If nsemiblock==noctant_per_block, then any value of nthread_octant
    =      applied to the OpenMP loop will work ok.
    =    - If nsemiblock<noctant_per_block==nthread_octant, then
    =      a potential race condition will occur.  This can be fixed by
    =      making the update of vo at the end of Sweeper_sweep_block atomic.
    =      What is in question here is the overhead of the semiblock loop.
    =      One might want to reduce the number of semiblocks while keeping
    =      noctant_per_block==nthread_octant high to get more thread
    =      parallelism but possibly not too high so as to control the
    =      wavefront latency.
    =========================================================================*/

    /*--------------------*/
    /*---Loop over semiblocks---*/
    /*--------------------*/

    for( semiblock=0; semiblock<sweeper->nsemiblock; ++semiblock )
    {
      /*---Initialize OpenMP thread number and thread count---*/

      const int thread_num  = 0;
      const int num_threads = 1;

      int octant_in_block = 0;

      /*--------------------*/
      /*---Loop over octants in octant block---*/
      /*--------------------*/

#ifdef USE_OPENMP_OCTANT
#pragma omp parallel for num_threads( sweeper->nthread_octant )
#endif
      for( octant_in_block=0; octant_in_block<sweeper->noctant_per_block;
                                                            ++octant_in_block )
      {
        /*---Compute thread information---*/

        const int thread_num_outer  = thread_num;
        const int num_threads_outer = num_threads;
        const int thread_num  = thread_num_outer + num_threads_outer *
                      ( IS_USING_OPENMP_OCTANT ? Env_thread_this( env ) : 0 );
        const int num_threads =  num_threads_outer *
                      ( IS_USING_OPENMP_OCTANT ? Env_num_threads( env ) : 1 );

        /*---Get step info---*/

        const Step_Info step_info = Step_Scheduler_step_info(
           &(sweeper->step_scheduler), step, octant_in_block, proc_x, proc_y );

        /*--------------------*/
        /*---Begin compute section---*/
        /*--------------------*/

        if( step_info.is_active )
        {
          const int iz_base = step_info.block_z * sweeper->dims_b.nz;

          const int dir_x = Dir_x( step_info.octant );
          const int dir_y = Dir_y( step_info.octant );
          const int dir_z = Dir_z( step_info.octant );

          /*--------------------*/
          /*---Compute semiblock bounds---*/
          /*--------------------*/

          /*===================================================================
          =    is_x_semiblocked: indicate whether the block is broken into
          =      semiblocks along the x axis.
          =    is_semiblock_x_lo: on this semiblock step for this thread
          =      do we process the lower or the higher semiblock along
          =      the x axis.  Only meaningful if is_x_semiblocked.
          =    has_x_lo: does this semiblock contain the lowest cell of the
          =      block along the x axis.
          =    has_x_hi: similarly.
          =    ixmin_b: the lowest cell boundary of the semiblock within the
          =      block along the x axis, inclusive of endpoints.
          =    ixmax_b: similarly.
          ===================================================================*/


          const Bool_t is_x_semiblocked = sweeper->nsemiblock > (1<<0);
          const Bool_t is_semiblock_x_lo = ( ( semiblock & (1<<0) ) == 0 ) ==
                                           ( dir_x == Dir_up() );

          const Bool_t has_x_lo =     is_semiblock_x_lo   || ! is_x_semiblocked;
          const Bool_t has_x_hi = ( ! is_semiblock_x_lo ) || ! is_x_semiblocked;

          const int ixmin_b =     has_x_lo ? 0
                                             : sweeper->dims_b.nx/2;
          const int ixmax_b = ( ! has_x_hi ) ? sweeper->dims_b.nx/2 - 1
                                             : sweeper->dims_b.nx   - 1;

          /*--------------------*/

          const Bool_t is_y_semiblocked = sweeper->nsemiblock > (1<<1);
          const Bool_t is_semiblock_y_lo = ( ( semiblock & (1<<1) ) == 0 ) ==
                                           ( dir_y == Dir_up() );

          const Bool_t has_y_lo =     is_semiblock_y_lo   || ! is_y_semiblocked;
          const Bool_t has_y_hi = ( ! is_semiblock_y_lo ) || ! is_y_semiblocked;

          const int iymin_b =     has_y_lo ? 0
                                             : sweeper->dims_b.ny/2;
          const int iymax_b = ( ! has_y_hi ) ? sweeper->dims_b.ny/2 - 1
                                             : sweeper->dims_b.ny   - 1;

          /*--------------------*/

          const Bool_t is_z_semiblocked = sweeper->nsemiblock > (1<<2);
          const Bool_t is_semiblock_z_lo = ( ( semiblock & (1<<2) ) == 0 ) ==
                                           ( dir_z == Dir_up() );

          const Bool_t has_z_lo =     is_semiblock_z_lo   || ! is_z_semiblocked;
          const Bool_t has_z_hi = ( ! is_semiblock_z_lo ) || ! is_z_semiblocked;

          const int izmin_b =     has_z_lo ? 0
                                             : sweeper->dims_b.nz/2;
          const int izmax_b = ( ! has_z_hi ) ? sweeper->dims_b.nz/2 - 1
                                             : sweeper->dims_b.nz   - 1;

          /*--------------------*/
          /*---Set physical boundary conditions if part of semiblock---*/
          /*--------------------*/

          if( ( dir_z == Dir_up() && step_info.block_z == 0 && has_z_lo ) ||
              ( dir_z == Dir_dn() && step_info.block_z ==
                                      sweeper->nblock_z - 1 && has_z_hi ) )
          {
            Sweeper_set_boundary_xy( sweeper, facexy, quan,
                                     step_info.octant, octant_in_block,
                                     ixmin_b, ixmax_b, iymin_b, iymax_b );
          }

          /*--------------------*/

          if( ( dir_y == Dir_up() && proc_y == 0 && has_y_lo ) ||
              ( dir_y == Dir_dn() && proc_y ==
                            Env_nproc_y( env )-1 && has_y_hi ) )
          {
            Sweeper_set_boundary_xz( sweeper, facexz, quan, step_info.block_z,
                                     step_info.octant, octant_in_block,
                                     ixmin_b, ixmax_b, izmin_b, izmax_b );
          }

          /*--------------------*/

          if( ( dir_x == Dir_up() && proc_x == 0 && has_x_lo ) ||
              ( dir_x == Dir_dn() && proc_x ==
                            Env_nproc_x( env )-1 && has_x_hi ) )
          {
            Sweeper_set_boundary_yz( sweeper, faceyz, quan, step_info.block_z,
                                     step_info.octant, octant_in_block,
                                     iymin_b, iymax_b, izmin_b, izmax_b );
          }

          /*--------------------*/
          /*---Perform sweep on relevant block---*/
          /*--------------------*/

          Sweeper_sweep_block( sweeper, vo, vi, facexy, facexz, faceyz,
                               quan, env, step_info,
                               thread_num, num_threads, octant_in_block,
                               ixmin_b,         ixmax_b,
                               iymin_b,         iymax_b,
                               izmin_b+iz_base, izmax_b+iz_base );

        }  /*---is_active---*/
      } /*---octant_in_block---*/
    } /*---semiblock---*/

    /*--------------------*/
    /*---Communicate faces---*/
    /*--------------------*/

    if( Sweeper_is_face_comm_async() )
    {
      Sweeper_send_faces_end__  (  sweeper, step-1, env );
      Sweeper_send_faces_start__(  sweeper, step, env );
    }
    else
    {
      Sweeper_communicate_faces__( sweeper, step, env );
    }

  } /*---step---*/

  /*---Increment message tag---*/

  Env_increment_tag( env, sweeper->noctant_per_block );

} /*---sweep---*/

/*===========================================================================*/

#endif /*---_serial_c__sweeper_kba_c_h_---*/

/*---------------------------------------------------------------------------*/
