/*---------------------------------------------------------------------------*/
/*!
 * \file   sweeper_kba_c.h
 * \author Wayne Joubert
 * \date   Tue Jan 28 16:37:41 EST 2014
 * \brief  Definitions for performing a sweep, kba version.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _sweeper_kba_c_h_
#define _sweeper_kba_c_h_

#include "function_attributes.h"
#include "env.h"
#include "definitions.h"
#include "quantities.h"
#include "array_accessors.h"
#include "array_operations.h"
#include "memory.h"
#include "step_scheduler_kba.h"
#include "sweeper_kba.h"

#include "sweeper_kba_c_kernels.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*===========================================================================*/
/*---Pseudo-constructor for Sweeper struct---*/

void Sweeper_ctor( Sweeper*          sweeper,
                   Dimensions        dims,
                   const Quantities* quan,
                   Env*              env,
                   Arguments*        args )
{
  /*====================*/
  /*---Declarations---*/
  /*====================*/

  Bool_t is_face_comm_async = Arguments_consume_int_or_default( args,
                                           "--is_face_comm_async", Bool_true );

  Insist( dims.nx > 0 ? "Currently requires all spatial blocks nonempty" : 0 );
  Insist( dims.ny > 0 ? "Currently requires all spatial blocks nonempty" : 0 );
  Insist( dims.nz > 0 ? "Currently requires all spatial blocks nonempty" : 0 );

  /*====================*/
  /*---Set up number of kba blocks---*/
  /*====================*/

  sweeper->nblock_z = Arguments_consume_int_or_default( args, "--nblock_z", 1);

  Insist( sweeper->nblock_z > 0 ? "Invalid z blocking factor supplied" : 0 );
  Insist( dims.nz % sweeper->nblock_z == 0
                  ? "Currently require all blocks have same z dimension" : 0 );

  /*====================*/
  /*---Set up number of octant threads---*/
  /*====================*/

  sweeper->nthread_octant
              = Arguments_consume_int_or_default( args, "--nthread_octant", 1);

  /*---Require a power of 2 between 1 and 8 inclusive---*/
  Insist( sweeper->nthread_octant>0 && sweeper->nthread_octant<=NOCTANT
          && ((sweeper->nthread_octant&(sweeper->nthread_octant-1))==0)
                                ? "Invalid octant thread count supplied" : 0 );
  /*---Don't allow threading in cases where it doesn't make sense---*/
  Insist( sweeper->nthread_octant==1 || IS_USING_OPENMP_THREADS
                                     || Env_cuda_is_using_device( env ) ?
          "Threading not allowed for this case" : 0 );

  sweeper->noctant_per_block = sweeper->nthread_octant;
  sweeper->nblock_octant     = NOCTANT / sweeper->noctant_per_block;

  /*====================*/
  /*---Set up number of semiblock steps---*/
  /*====================*/

  sweeper->nsemiblock = Arguments_consume_int_or_default(
                               args, "--nsemiblock", sweeper->nthread_octant );

  Insist( sweeper->nsemiblock>0 && sweeper->nsemiblock<=NOCTANT
          && ((sweeper->nsemiblock&(sweeper->nsemiblock-1))==0)
                                ? "Invalid semiblock count supplied" : 0 );
  Insist( ( sweeper->nsemiblock >= sweeper->nthread_octant ||
            IS_USING_OPENMP_VO_ATOMIC )
         ? "Incomplete set of semiblock steps requires atomic vo update" : 0 );

  /*====================*/
  /*---Set up number of energy threads---*/
  /*====================*/

  sweeper->nthread_e
                   = Arguments_consume_int_or_default( args, "--nthread_e", 1);

  Insist( sweeper->nthread_e > 0 ? "Invalid e thread count supplied." : 0 );
  /*---Don't allow threading in cases where it doesn't make sense---*/
  Insist( sweeper->nthread_e==1 || IS_USING_OPENMP_THREADS
                                || Env_cuda_is_using_device( env ) ?
          "Threading not allowed for this case" : 0 );

  /*====================*/
  /*---Set up step scheduler---*/
  /*====================*/

  Step_Scheduler_ctor( &(sweeper->step_scheduler),
                              sweeper->nblock_z, sweeper->nblock_octant, env );

  /*====================*/
  /*---Set up amu threads---*/
  /*====================*/

  Insist( NU * 1 > 0 );
  Insist(      Sweeper_nthread_u( sweeper, env ) > 0 );
  Insist( NU % Sweeper_nthread_u( sweeper, env ) == 0 );
  Insist( Sweeper_nthread_a( sweeper, env ) > 0 );
  if( ! IS_USING_MIC )
  {
    Insist( Sweeper_nthread_a( sweeper, env ) %
            Sweeper_nthread_u( sweeper, env ) == 0 );
    Insist( Sweeper_nthread_a( sweeper, env ) ==
            Sweeper_nthread_u( sweeper, env ) *
            Sweeper_nthread_m( sweeper, env ) );
  }
  if( IS_USING_MIC )
  {
    /*---For alignment, make this assumption.  For user case, assume this
         may mean some padding---*/
    Insist( dims.na % VEC_LEN == 0 );
  }

  /*====================*/
  /*---Set up dims structs---*/
  /*====================*/

  sweeper->dims = dims;

  sweeper->dims_b = sweeper->dims;
  sweeper->dims_b.nz = sweeper->dims.nz / sweeper->nblock_z;

  sweeper->dims_g = sweeper->dims;
  sweeper->dims_g.nx = quan->nx_g;
  sweeper->dims_g.ny = quan->ny_g;

  /*====================*/
  /*---Allocate arrays---*/
  /*====================*/

  sweeper->vilocal_host__ = Env_cuda_is_using_device( env ) ?
                            ( (P*) NULL ) :
                            malloc_host_P( Sweeper_nvilocal__( sweeper, env ) );

  sweeper->vslocal_host__ = Env_cuda_is_using_device( env ) ?
                            ( (P*) NULL ) :
                            malloc_host_P( Sweeper_nvslocal__( sweeper, env ) );

  sweeper->volocal_host__ = Env_cuda_is_using_device( env ) ?
                            ( (P*) NULL ) :
                            malloc_host_P( Sweeper_nvolocal__( sweeper, env ) );

  /*====================*/
  /*---Allocate faces---*/
  /*====================*/

  Faces_ctor( &(sweeper->faces), sweeper->dims_b,
                sweeper->noctant_per_block, is_face_comm_async, env );
}

/*===========================================================================*/
/*---Pseudo-destructor for Sweeper struct---*/

void Sweeper_dtor( Sweeper* sweeper,
                   Env*     env )
{
  /*====================*/
  /*---Deallocate arrays---*/
  /*====================*/

  if( ! Env_cuda_is_using_device( env ) )
  {
    if( sweeper->vilocal_host__ )
    {
      free_host_P( sweeper->vilocal_host__ );
    }
    if( sweeper->vslocal_host__ )
    {
      free_host_P( sweeper->vslocal_host__ );
    }
    if( sweeper->volocal_host__ )
    {
      free_host_P( sweeper->volocal_host__ );
    }
    sweeper->vilocal_host__ = NULL;
    sweeper->vslocal_host__ = NULL;
    sweeper->volocal_host__ = NULL;
  }

  /*====================*/
  /*---Deallocate faces---*/
  /*====================*/

  Faces_dtor( &(sweeper->faces) );

  /*====================*/
  /*---Terminate scheduler---*/
  /*====================*/

  Step_Scheduler_dtor( &( sweeper->step_scheduler ) );
}

/*===========================================================================*/
/*---Extract Sweeper_Lite from Sweeper---*/

Sweeper_Lite Sweeper_sweeper_lite( Sweeper sweeper )
{
  Sweeper_Lite sweeper_lite;

  sweeper_lite.vilocal_host__ = sweeper.vilocal_host__;
  sweeper_lite.vslocal_host__ = sweeper.vslocal_host__;
  sweeper_lite.volocal_host__ = sweeper.volocal_host__;

  sweeper_lite.dims   = sweeper.dims;
  sweeper_lite.dims_b = sweeper.dims_b;
  sweeper_lite.dims_g = sweeper.dims_g;

  sweeper_lite.nthread_e      = sweeper.nthread_e;
  sweeper_lite.nthread_octant = sweeper.nthread_octant;

  sweeper_lite.nblock_z          = sweeper.nblock_z;
  sweeper_lite.nblock_octant     = sweeper.nblock_octant;
  sweeper_lite.noctant_per_block = sweeper.noctant_per_block;
  sweeper_lite.nsemiblock        = sweeper.nsemiblock;

  return sweeper_lite;
}

/*===========================================================================*/
/*---Perform a sweep for a block, implementation, global---*/

TARGET_G void Sweeper_sweep_block_adapter(
  Sweeper_Lite           sweeper,
        P* __restrict__  vo,
  const P* __restrict__  vi,
        P* __restrict__  facexy,
        P* __restrict__  facexz,
        P* __restrict__  faceyz,
  const P* __restrict__  a_from_m,
  const P* __restrict__  m_from_a,
  int                    step,
  const Quantities       quan,
  Bool_t                 proc_x_min,
  Bool_t                 proc_x_max,
  Bool_t                 proc_y_min,
  Bool_t                 proc_y_max,
  Step_Info_Values       step_info_values,
  unsigned long int      do_block_init )
{
    Sweeper_sweep_block_impl( &sweeper, vo, vi, facexy, facexz, faceyz,
                              a_from_m, m_from_a, step, &quan,
                              proc_x_min, proc_x_max, proc_y_min, proc_y_max,
                              step_info_values, do_block_init );
}

/*===========================================================================*/
/*---Perform a sweep for a block---*/

void Sweeper_sweep_block(
  Sweeper*               sweeper,
  Pointer*               vo,
  Pointer*               vi,
  int*                   is_block_init,
  Pointer*               facexy,
  Pointer*               facexz,
  Pointer*               faceyz,
  const Pointer*         a_from_m,
  const Pointer*         m_from_a,
  int                    step,
  const Quantities*      quan,
  Env*                   env )
{
  /*---Declarations---*/

  const int proc_x = Env_proc_x_this( env );
  const int proc_y = Env_proc_y_this( env );

  const int noctant_per_block = sweeper->noctant_per_block;

  Step_Info_Values step_info_values;
                                /*---But only use noctant_per_block values---*/

  int octant_in_block = 0;

  int semiblock = 0;

  unsigned long int do_block_init = 0;

  /*---Precalculate step_info for required octants---*/

  for( octant_in_block=0; octant_in_block<sweeper->noctant_per_block;
                                                            ++octant_in_block )
  {
    step_info_values.step_info[octant_in_block] = Step_Scheduler_step_info(
      &(sweeper->step_scheduler), step, octant_in_block, proc_x, proc_y );

  }

  /*---Precalculate initialization schedule---*/

  for( semiblock=0; semiblock<sweeper->nsemiblock; ++semiblock )
  {
#pragma novector
    for( octant_in_block=0; octant_in_block<sweeper->noctant_per_block;
                                                            ++octant_in_block )
    {
      const Step_Info step_info = step_info_values.step_info[octant_in_block];
      if( step_info.is_active )
      {
        const int dir_x = Dir_x( step_info.octant );
        const int dir_y = Dir_y( step_info.octant );
        const int dir_z = Dir_z( step_info.octant );
        const Bool_t is_x_semiblocked = sweeper->nsemiblock > (1<<0);
        const Bool_t is_semiblock_x_lo = ( ( semiblock & (1<<0) ) == 0 ) ==
                                         ( dir_x == DIR_UP );
        const Bool_t has_x_lo =     is_semiblock_x_lo   || ! is_x_semiblocked;
        const Bool_t is_y_semiblocked = sweeper->nsemiblock > (1<<1);
        const Bool_t is_semiblock_y_lo = ( ( semiblock & (1<<1) ) == 0 ) ==
                                         ( dir_y == DIR_UP );
        const Bool_t has_y_lo =     is_semiblock_y_lo   || ! is_y_semiblocked;
        const Bool_t is_z_semiblocked = sweeper->nsemiblock > (1<<2);
        const Bool_t is_semiblock_z_lo = ( ( semiblock & (1<<2) ) == 0 ) ==
                                         ( dir_z == DIR_UP );
        const Bool_t has_z_lo =     is_semiblock_z_lo   || ! is_z_semiblocked;
        const int semiblock_num = ( has_x_lo ? 0 : 1 ) + 2 * (
                                  ( has_y_lo ? 0 : 1 ) + 2 * (
                                  ( has_z_lo ? 0 : 1 ) ));
        if( ! ( is_block_init[ step_info.block_z ] & ( 1 << semiblock_num ) ) )
        {
          do_block_init |= ( ((unsigned long int)1) <<
                             ( octant_in_block + noctant_per_block *
                               semiblock ) );
          is_block_init[ step_info.block_z ] |= ( 1 << semiblock_num );
        }
      }
    } /*---octant_in_block---*/
  } /*---semiblock---*/

  /*---Call sweep block implementation function---*/

  Sweeper_Lite sweeper_lite = Sweeper_sweeper_lite( *sweeper );

  if( Env_cuda_is_using_device( env ) )
  {
    Sweeper_sweep_block_adapter
#ifdef __CUDACC__
                 <<< dim3( Sweeper_nthreadblock( sweeper, 0, env ),
                           Sweeper_nthreadblock( sweeper, 1, env ),
                           Sweeper_nthreadblock( sweeper, 2, env ) ),
                     dim3( Sweeper_nthread_in_threadblock( sweeper, 0, env ),
                           Sweeper_nthread_in_threadblock( sweeper, 1, env ),
                           Sweeper_nthread_in_threadblock( sweeper, 2, env ) ),
                     Sweeper_shared_size__( sweeper, env ),
                     Env_cuda_stream_kernel_faces( env )
                 >>>
#endif
                            ( sweeper_lite,
                              Pointer_d( vo ),
                              Pointer_d( vi ),
                              Pointer_d( facexy ),
                              Pointer_d( facexz ),
                              Pointer_d( faceyz ),
                              Pointer_const_d( a_from_m ),
                              Pointer_const_d( m_from_a ),
                              step, *quan,
                              proc_x==0,
                              proc_x==Env_nproc_x( env )-1,
                              proc_y==0,
                              proc_y==Env_nproc_y( env )-1,
                              step_info_values,
                              do_block_init );
    Assert( Env_cuda_last_call_succeeded() );
  }
  else
#ifdef USE_OPENMP_THREADS
#pragma omp parallel num_threads( sweeper->nthread_e * sweeper->nthread_octant )
#endif
  {
    Sweeper_sweep_block_impl( &sweeper_lite,
                              Pointer_h( vo ),
                              Pointer_h( vi ),
                              Pointer_h( facexy ),
                              Pointer_h( facexz ),
                              Pointer_h( faceyz ),
                              Pointer_const_h( a_from_m ),
                              Pointer_const_h( m_from_a ),
                              step, quan,
                              proc_x==0,
                              proc_x==Env_nproc_x( env )-1,
                              proc_y==0,
                              proc_y==Env_nproc_y( env )-1,
                              step_info_values,
                              do_block_init );

  } /*---OPENMP---*/
}

/*===========================================================================*/
/*---Perform a sweep---*/

void Sweeper_sweep(
  Sweeper*               sweeper,
  Pointer*               vo,
  Pointer*               vi,
  const Quantities*      quan,
  Env*                   env )
{
  Assert( sweeper );
  Assert( vi );
  Assert( vo );

  /*---Declarations---*/

  const int nblock_z = sweeper->nblock_z;

  const int nstep = Step_Scheduler_nstep( &(sweeper->step_scheduler) );
  int step = -1;

  const size_t size_state_block = Dimensions_size_state( sweeper->dims, NU )
                                                                   / nblock_z;

  Bool_t* is_block_init = (Bool_t*) malloc( nblock_z * sizeof( Bool_t ) );

  int i = 0;

  for( i=0; i<nblock_z; ++i )
  {
    is_block_init[i] = 0;
  }

  /*---Initialize result array to zero if needed---*/

  if( sweeper->nsemiblock < sweeper->nthread_octant )
  {
    initialize_state_zero( Pointer_h( vo ), sweeper->dims, NU );
    Pointer_update_d_stream( vo, Env_cuda_stream_kernel_faces( env ) );
  }

  /*--------------------*/
  /*---Loop over kba parallel steps---*/
  /*--------------------*/

  for( step=0-1; step<nstep+1; ++step )
  {
    const Bool_t is_sweep_step = step>=0 && step<nstep;

    Pointer vi_b = Pointer_null();
    Pointer vo_b = Pointer_null();

    int i = 0;

    /*---Pick up needed face pointers---*/

    /*=========================================================================
    =    Order is important here.
    =    The _r face for a step must match the _c face for the next step.
    =    The _s face for a step must match the _c face for the prev step.
    =========================================================================*/

    Pointer* facexy = Faces_facexy_step( &(sweeper->faces), step );
    Pointer* facexz = Faces_facexz_step( &(sweeper->faces), step );
    Pointer* faceyz = Faces_faceyz_step( &(sweeper->faces), step );

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

    /*====================*/
    /*---Recv face via MPI WAIT (i)---*/
    /*====================*/

    if( is_sweep_step &&  Faces_is_face_comm_async( &(sweeper->faces)) )
    {
      Faces_recv_faces_end( &(sweeper->faces), &(sweeper->step_scheduler),
                            sweeper->dims_b, step-1, env );
    }

    /*====================*/
    /*---Send face to device START (i)---*/
    /*---Send face to device WAIT (i)---*/
    /*====================*/

    if( is_sweep_step )
    {
      if( step == 0 )
      {
        Pointer_update_d_stream( facexy, Env_cuda_stream_kernel_faces( env ) );
      }
      Pointer_update_d_stream(   facexz, Env_cuda_stream_kernel_faces( env ) );
      Pointer_update_d_stream(   faceyz, Env_cuda_stream_kernel_faces( env ) );
    }
    Env_cuda_stream_wait( Env_cuda_stream_kernel_faces( env ) );

    /*====================*/
    /*---Recv face via MPI START (i+1)---*/
    /*====================*/

    if( is_sweep_step &&  Faces_is_face_comm_async( &(sweeper->faces)) )
    {
      Faces_recv_faces_start( &(sweeper->faces), &(sweeper->step_scheduler),
                            sweeper->dims_b, step, env );
    }

    /*====================*/
    /*---Perform the sweep on the block START (i)---*/
    /*====================*/

    if( is_sweep_step )
    {
      Sweeper_sweep_block( sweeper, vo, vi, is_block_init,
                           facexy, facexz, faceyz,
                           & quan->a_from_m, & quan->m_from_a,
                           step, quan, env );
    }

    /*====================*/
    /*---Send block to device START (i+1)---*/
    /*====================*/

    for( i=0; i<2; ++i )
    {
      /*---Determine blocks needing transfer, counting from top/bottom z---*/
      /*---NOTE: for case of one octant thread, can speed this up by only
           send/recv of one block per step, not two---*/

      const int stept = step + 1;
      const int    block_to_send[2] = {                                stept,
                                        ( nblock_z-1 ) -               stept };
      const Bool_t do_block_send[2] = { block_to_send[0] <  nblock_z/2,
                                        block_to_send[1] >= nblock_z/2 };
      Assert( nstep >= nblock_z );  /*---Sanity check---*/
      if( do_block_send[i] )
      {
        Pointer_ctor_alias(      &vi_b, vi, size_state_block * block_to_send[i],
                                            size_state_block );
        Pointer_update_d_stream( &vi_b, Env_cuda_stream_send_block( env ) );
        Pointer_dtor(            &vi_b );

        /*---Initialize result array to zero if needed---*/
        /*---NOTE: this is not performance-optimal---*/
        if( sweeper->nsemiblock < sweeper->nthread_octant )
        {
          Pointer_ctor_alias(    &vo_b, vi, size_state_block * block_to_send[i],
                                            size_state_block );
          initialize_state_zero( Pointer_h( &vo_b ), sweeper->dims, NU );
          Pointer_update_d_stream( &vo_b, Env_cuda_stream_send_block( env ) );
          Pointer_dtor(            &vo_b );
        }
      }
    }

    /*====================*/
    /*---Recv block from device START (i-1)---*/
    /*====================*/

    for( i=0; i<2; ++i )
    {
      /*---Determine blocks needing transfer, counting from top/bottom z---*/
      /*---NOTE: for case of one octant thread, can speed this up by only
           send/recv of one block per step, not two---*/

      const int stept = step - 1;
      const int    block_to_recv[2] = { ( nblock_z-1 ) - ( nstep-1 - stept ),
                                                         ( nstep-1 - stept ) };
      const Bool_t do_block_recv[2] = { block_to_recv[0] >= nblock_z/2,
                                        block_to_recv[1] <  nblock_z/2 };
      Assert( nstep >= nblock_z );  /*---Sanity check---*/
      if( do_block_recv[i] )
      {
        Pointer_ctor_alias(      &vo_b, vo, size_state_block * block_to_recv[i],
                                            size_state_block );
        Pointer_update_h_stream( &vo_b, Env_cuda_stream_recv_block( env ) );
        Pointer_dtor(            &vo_b );
      }
    }

    /*====================*/
    /*---Send block to device WAIT (i+1)---*/
    /*---Recv block from device WAIT (i-1)---*/
    /*====================*/

    Env_cuda_stream_wait( Env_cuda_stream_send_block( env ) );
    Env_cuda_stream_wait( Env_cuda_stream_recv_block( env ) );

    /*====================*/
    /*---Send face via MPI WAIT (i-1)---*/
    /*====================*/

    if( is_sweep_step && Faces_is_face_comm_async( &(sweeper->faces)) )
    {
      Faces_send_faces_end( &(sweeper->faces), &(sweeper->step_scheduler),
                            sweeper->dims_b, step-1, env );
    }

    /*====================*/
    /*---Perform the sweep on the block WAIT (i)---*/
    /*====================*/

    Env_cuda_stream_wait( Env_cuda_stream_kernel_faces( env ) );

    /*====================*/
    /*---Recv face from device START (i)---*/
    /*---Recv face from device WAIT (i)---*/
    /*====================*/

    if( is_sweep_step )
    {
      if( step == nstep-1 )
      {
        Pointer_update_h_stream( facexy, Env_cuda_stream_kernel_faces( env ) );
      }
      Pointer_update_h_stream(   facexz, Env_cuda_stream_kernel_faces( env ) );
      Pointer_update_h_stream(   faceyz, Env_cuda_stream_kernel_faces( env ) );
    }
    Env_cuda_stream_wait( Env_cuda_stream_kernel_faces( env ) );

    /*====================*/
    /*---Send face via MPI START (i)---*/
    /*====================*/

    if( is_sweep_step && Faces_is_face_comm_async( &(sweeper->faces)) )
    {
      Faces_send_faces_start( &(sweeper->faces), &(sweeper->step_scheduler),
                            sweeper->dims_b, step, env );
    }

    /*====================*/
    /*---Communicate faces (synchronous)---*/
    /*====================*/

    if( is_sweep_step && ! Faces_is_face_comm_async( &(sweeper->faces)) )
    {
      Faces_communicate_faces( &(sweeper->faces), &(sweeper->step_scheduler),
                            sweeper->dims_b, step, env );
    }

  } /*---step---*/

  /*---Increment message tag---*/

  Env_increment_tag( env, sweeper->noctant_per_block );

  /*---Finish---*/

  free( (void*) is_block_init );

} /*---sweep---*/

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

#endif /*---_sweeper_kba_c_h_---*/

/*---------------------------------------------------------------------------*/
