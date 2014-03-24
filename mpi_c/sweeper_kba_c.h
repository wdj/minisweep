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
  Insist( sweeper->nblock_z > 0 && "Invalid z blocking factor supplied." );
  Insist( dims.nz % sweeper->nblock_z == 0 &&
          "KBA sweeper currently requires all blocks have same z dimension" );

  /*---Set up number of octant threads---*/
  sweeper->nthread_octant
              = Arguments_consume_int_or_default( args, "--nthread_octant", 1);
  /*---Require a power of 2 between 1 and 8 inclusive---*/
  Insist( sweeper->nthread_octant>0 && sweeper->nthread_octant<=NOCTANT
          && ((sweeper->nthread_octant&(sweeper->nthread_octant-1))==0)
          && "Invalid octant thread count supplied." );

  /*---Set up number of energy threads---*/
  sweeper->nthread_e
                   = Arguments_consume_int_or_default( args, "--nthread_e", 1);
  Insist( sweeper->nthread_e > 0 && "Invalid e thread count supplied." );

  sweeper->nblock_octant = NOCTANT / sweeper->nthread_octant;
  Step_Scheduler_ctor( &(sweeper->step_scheduler),
                              sweeper->nblock_z, sweeper->nblock_octant, env );

  sweeper->dims = dims;

  sweeper->dims_b = sweeper->dims;
  sweeper->dims_b.nz = sweeper->dims.nz / sweeper->nblock_z;

  sweeper->dims_g = sweeper->dims;
  sweeper->dims_g.nx = quan->nx_g;
  sweeper->dims_g.ny = quan->ny_g;

  /*---Allocate arrays---*/

  sweeper->v_local = malloc_P( sweeper->dims_b.na * NU * sweeper->nthread_e );

  sweeper->facexy0 = malloc_P( Dimensions_size_facexy( sweeper->dims_b, NU,
                                      Sweeper_noctant_per_block( sweeper ) ) );
  sweeper->facexz0 = malloc_P( Dimensions_size_facexz( sweeper->dims_b, NU,
                                      Sweeper_noctant_per_block( sweeper ) ) );
  sweeper->faceyz0 = malloc_P( Dimensions_size_faceyz( sweeper->dims_b, NU,
                                      Sweeper_noctant_per_block( sweeper ) ) );

  sweeper->facexz1 = NULL;
  sweeper->facexz2 = NULL;
  sweeper->faceyz1 = NULL;
  sweeper->faceyz2 = NULL;

  if( Sweeper_is_face_comm_async() )
  {
    sweeper->facexz1 = malloc_P( Dimensions_size_facexz( sweeper->dims_b, NU,
                                      Sweeper_noctant_per_block( sweeper ) ) );
    sweeper->facexz2 = malloc_P( Dimensions_size_facexz( sweeper->dims_b, NU,
                                      Sweeper_noctant_per_block( sweeper ) ) );
    sweeper->faceyz1 = malloc_P( Dimensions_size_faceyz( sweeper->dims_b, NU,
                                      Sweeper_noctant_per_block( sweeper ) ) );
    sweeper->faceyz2 = malloc_P( Dimensions_size_faceyz( sweeper->dims_b, NU,
                                      Sweeper_noctant_per_block( sweeper ) ) );
  }
}

/*===========================================================================*/
/*---Pseudo-destructor for Sweeper struct---*/

void Sweeper_dtor( Sweeper* sweeper )
{
  /*---Deallocate arrays---*/

  free_P( sweeper->v_local );

  free_P( sweeper->facexy0 );
  free_P( sweeper->facexz0 );
  free_P( sweeper->faceyz0 );

  if( Sweeper_is_face_comm_async() )
  {
    free_P( sweeper->facexz1 );
    free_P( sweeper->facexz2 );
    free_P( sweeper->faceyz1 );
    free_P( sweeper->faceyz2 );
  }

  sweeper->v_local = NULL;

  sweeper->facexy0 = NULL;
  sweeper->facexz0 = NULL;
  sweeper->faceyz0 = NULL;

  sweeper->facexz1 = NULL;
  sweeper->facexz2 = NULL;
  sweeper->faceyz1 = NULL;
  sweeper->faceyz2 = NULL;

  Step_Scheduler_dtor( &( sweeper->step_scheduler ) );
}

/*===========================================================================*/
/*---Determine whether to send a face computed at step, used at step+1---*/

Bool_t Sweeper_must_do_send__(
  Sweeper*           sweeper,
  int                step,
  int                axis,
  int                dir_ind,
  Env*               env )
{
  const int octant_in_block = 0;
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
  Env*               env )
{
  const int octant_in_block = 0;
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

  const size_t size_facexz = Dimensions_size_facexz( sweeper->dims_b, NU,
                                        Sweeper_noctant_per_block( sweeper ) );
  const size_t size_faceyz = Dimensions_size_faceyz( sweeper->dims_b, NU,
                                        Sweeper_noctant_per_block( sweeper ) );

  /*---Allocate temporary face buffers---*/

  P* __restrict__ buf_facexz  = malloc_P( size_facexz );
  P* __restrict__ buf_faceyz  = malloc_P( size_faceyz );

  /*---Communicate +/-X, +/-Y---*/

  int axis = 0;

  for( axis=0; axis<2; ++axis )  /*---Loop: X, Y---*/
  {
    const Bool_t axis_x = axis==0;
    const Bool_t axis_y = axis==1;

    const int proc_axis = axis_x ? proc_x : proc_y;

    const size_t    size_face    = axis_x ? size_faceyz      : size_facexz;
    P* __restrict__ buf_face     = axis_x ? buf_faceyz       : buf_facexz;
    P* __restrict__ sweeper_face = axis_x ? sweeper->faceyz0 : sweeper->facexz0;
    int dir_ind = 0;

    for( dir_ind=0; dir_ind<2; ++dir_ind ) /*---Loop: up, down---*/
    {
      const int dir = dir_ind==0 ? Dir_up() : Dir_dn();
      const int inc_x = axis_x ? Dir_inc( dir ) : 0;
      const int inc_y = axis_y ? Dir_inc( dir ) : 0;

      /*---Determine whether to communicate---*/

      Bool_t const do_send = Sweeper_must_do_send__(
                                           sweeper, step, axis, dir_ind, env );

      Bool_t const do_recv = Sweeper_must_do_recv__(
                                           sweeper, step, axis, dir_ind, env );

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

  const size_t size_facexz = Dimensions_size_facexz( sweeper->dims_b, NU,
                                        Sweeper_noctant_per_block( sweeper ) );
  const size_t size_faceyz = Dimensions_size_faceyz( sweeper->dims_b, NU,
                                        Sweeper_noctant_per_block( sweeper ) );

  /*---Communicate +/-X, +/-Y---*/

  int axis = 0;

  for( axis=0; axis<2; ++axis )
  {
    const Bool_t axis_x = axis==0;
    const Bool_t axis_y = axis==1;

    const int proc_axis = axis_x ? proc_x : proc_y;

    /*---Send values computed on this step---*/

    const size_t    size_face         = axis_x ? size_faceyz : size_facexz;
    P* __restrict__ sweeper_face_send = axis_x ?
                                             Sweeper_faceyz__( sweeper, step )
                                           : Sweeper_facexz__( sweeper, step );

    int dir_ind = 0;

    for( dir_ind=0; dir_ind<2; ++dir_ind )
    {
      const int dir = dir_ind==0 ? Dir_up() : Dir_dn();
      const int inc_x = axis_x ? Dir_inc( dir ) : 0;
      const int inc_y = axis_y ? Dir_inc( dir ) : 0;

      /*---Determine whether to communicate---*/

      Bool_t const do_send = Sweeper_must_do_send__(
                                           sweeper, step, axis, dir_ind, env );

      if( do_send )
      {
        const int proc_other = Env_proc( env, proc_x+inc_x, proc_y+inc_y );
        Request_t* request = axis_x ? & sweeper->request_send_xz
                                    : & sweeper->request_send_yz;
        Env_asend_P(
                 sweeper_face_send, size_face, proc_other, env->tag, request );
      }

    } /*---dir_ind---*/
  } /*---axis---*/
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
                                           sweeper, step, axis, dir_ind, env );

      if( do_send )
      {
        Request_t* request = axis_x ? & sweeper->request_send_xz
                                    : & sweeper->request_send_yz;
        Env_wait( request );
      }

    } /*---dir_ind---*/
  } /*---axis---*/
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

  const size_t size_facexz = Dimensions_size_facexz( sweeper->dims_b, NU,
                                        Sweeper_noctant_per_block( sweeper ) );
  const size_t size_faceyz = Dimensions_size_faceyz( sweeper->dims_b, NU,
                                        Sweeper_noctant_per_block( sweeper ) );

  /*---Communicate +/-X, +/-Y---*/

  int axis = 0;

  for( axis=0; axis<2; ++axis )
  {
    const Bool_t axis_x = axis==0;
    const Bool_t axis_y = axis==1;

    const int proc_axis = axis_x ? proc_x : proc_y;

    /*---Receive values computed on the next step---*/

    const size_t    size_face         = axis_x ? size_faceyz     : size_facexz;
    P* __restrict__ sweeper_face_recv = axis_x ?
                                           Sweeper_faceyz__( sweeper, step+1 )
                                         : Sweeper_facexz__( sweeper, step+1 );
    int dir_ind = 0;

    for( dir_ind=0; dir_ind<2; ++dir_ind )
    {
      const int dir = dir_ind==0 ? Dir_up() : Dir_dn();
      const int inc_x = axis_x ? Dir_inc( dir ) : 0;
      const int inc_y = axis_y ? Dir_inc( dir ) : 0;

      /*---Determine whether to communicate---*/

      Bool_t const do_recv = Sweeper_must_do_recv__(
                                           sweeper, step, axis, dir_ind, env );

      if( do_recv )
      {
        const int proc_other
                         = Env_proc( env, proc_x-inc_x, proc_y-inc_y );
        Request_t* request = axis_x ? & sweeper->request_recv_xz
                                    : & sweeper->request_recv_yz;
        Env_arecv_P(
                sweeper_face_recv, size_face, proc_other, env->tag, request );
      }

    } /*---dir_ind---*/
  } /*---axis---*/
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
                                           sweeper, step, axis, dir_ind, env );

      if( do_recv )
      {
        Request_t* request = axis_x ? & sweeper->request_recv_xz
                                    : & sweeper->request_recv_yz;
        Env_wait( request );
      }

    } /*---dir_ind---*/
  } /*---axis---*/
}

/*===========================================================================*/
/*---Apply boundary condition: xy face---*/

static void Sweeper_set_boundary_xy(
  const Sweeper*        sweeper,
  const Quantities*     quan,
  P* const __restrict__ facexy,
  int                   octant,
  int                   octant_in_block )
{
  const int ix_base = quan->ix_base;
  const int iy_base = quan->iy_base;
  const int dir_z = Dir_z( octant );

  int ix_b = 0;
  int ix_g = 0;
  int iy_b = 0;
  int iy_g = 0;
  int ia   = 0;
  int ie   = 0;
  int iu   = 0;

  const int iz_g = dir_z == Dir_up() ? -1 : sweeper->dims_g.nz;

#ifdef USE_OPENMP_E
#pragma omp parallel for num_threads( sweeper->nthread_e )
#endif
  for( ie=0; ie<sweeper->dims_b.ne; ++ie )
  {
  for( iu=0; iu<NU; ++iu )
  {
  for( iy_b=0; iy_b<sweeper->dims_b.ny; ++iy_b )
  {
    const int iy_g = iy_b + iy_base;
  for( ix_b=0; ix_b<sweeper->dims_b.nx; ++ix_b )
  {
    const int ix_g = ix_b + ix_base;
  for( ia=0; ia<sweeper->dims_b.na; ++ia )
  {
    *ref_facexy( facexy, sweeper->dims_b, NU,
                 Sweeper_noctant_per_block( sweeper ),
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
  const Quantities*     quan,
  P* const __restrict__ facexz,
  int                   octant,
  int                   block_z, 
  int                   octant_in_block )
{
  const int ix_base = quan->ix_base;
  const int iz_base = block_z * sweeper->dims_b.nz;
  const int dir_y = Dir_y( octant );
  const int nz_b = sweeper->dims_b.nz;

  int ix_b = 0;
  int ix_g = 0;
  int iz_b = 0;
  int iz_g = 0;
  int ia   = 0;
  int ie   = 0;
  int iu   = 0;

  const int iy_g = dir_y == Dir_up() ? -1 : sweeper->dims_g.ny;

#ifdef USE_OPENMP_E
#pragma omp parallel for num_threads( sweeper->nthread_e )
#endif
  for( ie=0; ie<sweeper->dims_b.ne; ++ie )
  {
  for( iu=0; iu<NU; ++iu )
  {
  for( iz_b=0; iz_b<nz_b; ++iz_b )
  {
    const int iz_g = iz_b + iz_base;
  for( ix_b=0; ix_b<sweeper->dims_b.nx; ++ix_b )
  {
    const int ix_g = ix_b + ix_base;
  for( ia=0; ia<sweeper->dims_b.na; ++ia )
  {
    *ref_facexz( facexz, sweeper->dims_b, NU,
                 Sweeper_noctant_per_block( sweeper ),
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
  const Quantities*     quan,
  P* const __restrict__ faceyz,
  int                   octant,
  int                   block_z, 
  int                   octant_in_block )
{
  const int iy_base = quan->iy_base;
  const int iz_base = block_z * sweeper->dims_b.nz;
  const int dir_x = Dir_x( octant );
  const int nz_b = sweeper->dims_b.nz;

  int iy_b = 0;
  int iy_g = 0;
  int iz_b = 0;
  int iz_g = 0;
  int ia   = 0;
  int ie   = 0;
  int iu   = 0;

  const int ix_g = dir_x == Dir_up() ? -1 : sweeper->dims_g.nx;

#ifdef USE_OPENMP_E
#pragma omp parallel for num_threads( sweeper->nthread_e )
#endif
  for( ie=0; ie<sweeper->dims_b.ne; ++ie )
  {
  for( iu=0; iu<NU; ++iu )
  {
  for( iz_b=0; iz_b<nz_b; ++iz_b )
  {
    const int iz_g = iz_b + iz_base;
  for( iy_b=0; iy_b<sweeper->dims_b.ny; ++iy_b )
  {
    const int iy_g = iy_b + iy_base;
  for( ia=0; ia<sweeper->dims_b.na; ++ia )
  {
    *ref_faceyz( faceyz, sweeper->dims_b, NU,
                 Sweeper_noctant_per_block( sweeper ),
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
  const Quantities*      quan,
  Env*                   env, 
  const Step_Info        step_info,
  const int              thread_num,
  const int              num_threads,
  const int              octant_in_block,
  P* __restrict__        facexy,
  P* __restrict__        facexz,
  P* __restrict__        faceyz )
{
  const int octant  = step_info.octant;
  const int block_z = step_info.block_z;
  const int iz_base = block_z * sweeper->dims_b.nz;
  const int dir_x = Dir_x( octant );
  const int dir_y = Dir_y( octant );
  const int dir_z = Dir_z( octant );

  int ie = 0;

#ifdef USE_OPENMP_E
#pragma omp parallel for num_threads( sweeper->nthread_e )
#endif
  for( ie=0; ie<sweeper->dims.ne; ++ie )
  {
    /*---Calculate spatial loop extents---*/

    const int ixbeg = dir_x==Dir_up() ? 0       : sweeper->dims_b.nx-1;
    const int iybeg = dir_y==Dir_up() ? 0       : sweeper->dims_b.ny-1;
    const int izbeg = dir_z==Dir_up() ? iz_base : iz_base+sweeper->dims_b.nz-1;

    const int ixend = dir_x==Dir_dn() ? 0       : sweeper->dims_b.nx-1;
    const int iyend = dir_y==Dir_dn() ? 0       : sweeper->dims_b.ny-1;
    const int izend = dir_z==Dir_dn() ? iz_base : iz_base+sweeper->dims_b.nz-1;

    /*---Compute thread information---*/

    const int thread_num_outer  = thread_num;
    const int num_threads_outer = num_threads;

    const int thread_num  = thread_num_outer + num_threads_outer *
                                                 Env_thread_this( env );
    const int num_threads =  num_threads_outer * Env_num_threads( env );

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
                        Sweeper_noctant_per_block( sweeper ),
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

  const int proc_x     = Env_proc_x_this( env );
  const int proc_y     = Env_proc_y_this( env );
  const int octant_in_block = 0;

  int step = -1;

  /*---Initialize result array to zero---*/

  initialize_state_zero( vo, sweeper->dims, NU );

  /*--------------------*/
  /*---Loop over kba parallel steps---*/
  /*--------------------*/

  for( step=0; step<Step_Scheduler_nstep( &(sweeper->step_scheduler) ); ++step )
  {
    /*---Get step info for this proc---*/

    const Step_Info step_info = Step_Scheduler_step_info(
           &(sweeper->step_scheduler), step, octant_in_block, proc_x, proc_y );

    /*---Pick up needed face pointers---*/

    /*=========================================================================
    =    Order is important here.
    =    The _r face for a step must match the _c face for the next step.
    =    The _s face for a step must match the _c face for the prev step.
    =========================================================================*/

    P* const __restrict__ facexy = Sweeper_facexy__( sweeper, step );
    P* const __restrict__ facexz = Sweeper_facexz__( sweeper, step );
    P* const __restrict__ faceyz = Sweeper_faceyz__( sweeper, step );

    /*---Initialize OpenMP thread number and thread count---*/

    int thread_num  = 0;
    int num_threads = 1;

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
      /*---Set face boundary conditions as needed---*/
      /*--------------------*/

      if( ( dir_z == Dir_up() && step_info.block_z == 0 ) ||
          ( dir_z == Dir_dn() && step_info.block_z == Step_Scheduler_nblock_z(
                                           &(sweeper->step_scheduler) ) - 1 ) )
      {
        Sweeper_set_boundary_xy( sweeper, quan, facexy,
                                           step_info.octant, octant_in_block );
      }

      /*--------------------*/

      if( ( dir_y == Dir_up() && proc_y == 0 ) ||
          ( dir_y == Dir_dn() && proc_y == Env_nproc_y( env )-1 ) )
      {
        Sweeper_set_boundary_xz( sweeper, quan, facexz,
                        step_info.octant, step_info.block_z, octant_in_block );
      }

      /*--------------------*/

      if( ( dir_x == Dir_up() && proc_x == 0 ) ||
          ( dir_x == Dir_dn() && proc_x == Env_nproc_x( env )-1 ) )
      {
        Sweeper_set_boundary_yz( sweeper, quan, faceyz,
                        step_info.octant, step_info.block_z, octant_in_block );
      }

      /*--------------------*/
      /*---Perform sweep on relevant block---*/
      /*--------------------*/

      Sweeper_sweep_block( sweeper, vo, vi, quan, env, step_info,
                           thread_num, num_threads, octant_in_block,
                           facexy, facexz, faceyz );

    }  /*---is_active---*/

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

  env->tag++;

} /*---sweep---*/

/*===========================================================================*/

#endif /*---_serial_c__sweeper_kba_c_h_---*/

/*---------------------------------------------------------------------------*/
