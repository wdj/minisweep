/*---------------------------------------------------------------------------*/
/*!
 * \file   sweeper_kba.h
 * \author Wayne Joubert
 * \date   Tue Jan 28 16:37:41 EST 2014
 * \brief  Declarations for performing a sweep, kba version.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _serial_c__sweeper_kba_h_
#define _serial_c__sweeper_kba_h_

#include "function_attributes.h"
#include "env.h"
#include "definitions.h"
#include "dimensions.h"
#include "arguments.h"
#include "pointer.h"
#include "quantities.h"
#include "step_scheduler_kba.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*===========================================================================*/
/*---Set up enums---*/

#ifdef USE_OPENMP_THREADS
enum{ IS_USING_OPENMP_THREADS = 1 };
#else
enum{ IS_USING_OPENMP_THREADS = 0 };
#endif

#ifdef USE_OPENMP_VO_ATOMIC
enum{ IS_USING_OPENMP_VO_ATOMIC = 1 };
#else
enum{ IS_USING_OPENMP_VO_ATOMIC = 0 };
#endif

/*===========================================================================*/
/*---Struct with pointers etc. used to perform sweep---*/

typedef struct
{
  Pointer          facexy0;

  Pointer          facexz0;
  Pointer          facexz1;
  Pointer          facexz2;

  Pointer          faceyz0;
  Pointer          faceyz1;
  Pointer          faceyz2;

  P* __restrict__  vslocal;

  Dimensions       dims;
  Dimensions       dims_b;
  Dimensions       dims_g;

  int              nblock_z;
  int              nblock_octant;
  int              nthread_e;
  int              nthread_octant;
  int              noctant_per_block;
  int              nsemiblock;

  Request_t        request_send_xz[NOCTANT];
  Request_t        request_send_yz[NOCTANT];
  Request_t        request_recv_xz[NOCTANT];
  Request_t        request_recv_yz[NOCTANT];

  Step_Scheduler   step_scheduler;
} Sweeper;

/*===========================================================================*/
/*---Pseudo-constructor for Sweeper struct---*/

void Sweeper_ctor( Sweeper*          sweeper,
                   Dimensions        dims,
                   const Quantities* quan,
                   Env*              env,
                   Arguments*        args );

/*===========================================================================*/
/*---Pseudo-destructor for Sweeper struct---*/

void Sweeper_dtor( Sweeper* sweeper,
                   Env*     env );

/*===========================================================================*/
/*---Number of octants in an octant block---*/

static int Sweeper_noctant_per_block( const Sweeper* sweeper )
{
  return sweeper->noctant_per_block;
}

/*===========================================================================*/
/*---Is face communication done asynchronously---*/

static int Sweeper_is_face_comm_async()
{
  return Bool_true;
}

/*===========================================================================*/
/*---Determine whether to send a face now---*/

Bool_t Sweeper_must_do_send__(
  Sweeper*           sweeper,
  int                step,
  int                axis,
  int                dir_ind,
  int                octant_in_block,
  Env*               env );

/*===========================================================================*/
/*---Determine whether to receive a face now---*/

Bool_t Sweeper_must_do_recv__(
  Sweeper*           sweeper,
  int                step,
  int                axis,
  int                dir_ind,
  int                octant_in_block,
  Env*               env );

/*===========================================================================*/
/*---Communicate faces---*/

void Sweeper_communicate_faces__(
  Sweeper*         sweeper,
  int              step,
  Env*             env );

/*---------------------------------------------------------------------------*/

void Sweeper_send_faces_start__(
  Sweeper*           sweeper,
  int                step,
  Env*               env );

/*---------------------------------------------------------------------------*/

void Sweeper_send_faces_end__(
  Sweeper*           sweeper,
  int                step,
  Env*               env );

/*---------------------------------------------------------------------------*/

void Sweeper_recv_faces_start__(
  Sweeper*           sweeper,
  int                step,
  Env*               env );

/*---------------------------------------------------------------------------*/

void Sweeper_recv_faces_end__(
  Sweeper*           sweeper,
  int                step,
  Env*               env );

/*===========================================================================*/
/*---Apply boundary condition: xy face---*/

TARGET_HD static void Sweeper_set_boundary_xy(
  const Sweeper*        sweeper,
  P* const __restrict__ facexy,
  const Quantities*     quan,
  int                   octant,
  int                   octant_in_block, 
  const int             ixmin_b,
  const int             ixmax_b,
  const int             iymin_b,
  const int             iymax_b );

/*===========================================================================*/
/*---Apply boundary condition: xz face---*/

TARGET_HD static void Sweeper_set_boundary_xz(
  const Sweeper*        sweeper,
  P* const __restrict__ facexz,
  const Quantities*     quan,
  int                   block_z,
  int                   octant,
  int                   octant_in_block, 
  const int             ixmin_b,
  const int             ixmax_b,
  const int             izmin_b,
  const int             izmax_b );

/*===========================================================================*/
/*---Apply boundary condition: yz face---*/

TARGET_HD static void Sweeper_set_boundary_yz(
  const Sweeper*        sweeper,
  P* const __restrict__ faceyz,
  const Quantities*     quan,
  int                   block_z,
  int                   octant,
  int                   octant_in_block,
  const int             iymin_b,
  const int             iymax_b,
  const int             izmin_b,
  const int             izmax_b );

/*===========================================================================*/
/*---Selectors for faces---*/

/*---The xz and yz face arrays form a circular buffer of length three.
     Three are needed because at any step there may be a send, a receive, and a
     block-sweep-compute in-flight.
---*/

static Pointer* Sweeper_facexy__( Sweeper* sweeper, int i )
{
  Assert( sweeper != NULL );
  Assert( i >= 0 && i < 1 );
  return & sweeper->facexy0;
}

/*---------------------------------------------------------------------------*/

static Pointer* Sweeper_facexz__( Sweeper* sweeper, int i )
{
  Assert( sweeper != NULL );
  Assert( i >= 0 && i < ( Sweeper_is_face_comm_async() ? NDIM : 1 ) );
  Pointer* facesxz[NDIM] = { & sweeper->facexz0,
                             & sweeper->facexz1,
                             & sweeper->facexz2 };
  return facesxz[i];
}

/*---------------------------------------------------------------------------*/

static Pointer* Sweeper_faceyz__( Sweeper* sweeper, int i )
{
  Assert( sweeper != NULL );
  Assert( i >= 0 && i < ( Sweeper_is_face_comm_async() ? NDIM : 1 ) );
  Pointer* facesyz[NDIM] = { & sweeper->faceyz0,
                             & sweeper->faceyz1,
                             & sweeper->faceyz2 };
  return facesyz[i];
}

/*---------------------------------------------------------------------------*/


static Pointer* Sweeper_facexy_step__( Sweeper* sweeper, int step )
{
  Assert( sweeper != NULL );
  Assert( step >= -1 );
  return Sweeper_facexy__( sweeper, 0 );
}

/*---------------------------------------------------------------------------*/

static Pointer* Sweeper_facexz_step__( Sweeper* sweeper, int step )
{
  Assert( sweeper != NULL );
  Assert( step >= -1 );

  return Sweeper_facexz__( sweeper,
                           Sweeper_is_face_comm_async() ? (step+3)%3 : 0 );
}

/*---------------------------------------------------------------------------*/

static Pointer* Sweeper_faceyz_step__( Sweeper* sweeper, int step )
{
  Assert( sweeper != NULL );
  Assert( step >= -1 );

  return Sweeper_faceyz__( sweeper,
                           Sweeper_is_face_comm_async() ? (step+3)%3 : 0 );
}

/*===========================================================================*/
/*---Thread indexers---*/

TARGET_HD static inline int Sweeper_thread_e( const Sweeper* sweeper )
{
#ifdef __CUDA_ARCH__
  return Env_cuda_threadblock( 0 );
#else
  Assert( sweeper->nthread_e * sweeper->nthread_octant == 1 ||
                                                        Env_omp_in_parallel() );
  return Env_omp_thread() % sweeper->nthread_e;
#endif
}

/*---------------------------------------------------------------------------*/

TARGET_HD static inline int Sweeper_thread_octant( const Sweeper* sweeper )
{
#ifdef __CUDA_ARCH__
  return Env_cuda_thread_in_threadblock( 0 );
#else
  Assert( sweeper->nthread_e * sweeper->nthread_octant == 1 ||
                                                       Env_omp_in_parallel() );
  return Env_omp_thread() / sweeper->nthread_e;
#endif
}

/*---------------------------------------------------------------------------*/

TARGET_HD static inline int Sweeper_thread( const Sweeper* sweeper )
{
  return Sweeper_thread_octant( sweeper ) + sweeper->nthread_octant *
         Sweeper_thread_e(      sweeper );
}

/*===========================================================================*/
/*---Synchronize across octant threads---*/

TARGET_HD static void Sweeper_sync_octant_threads( Sweeper* sweeper )
{
#ifdef __CUDA_ARCH__
  /*---NOTE: this is not needed if octant threads are mapped in-warp---*/
  Env_cuda_sync_threadblock();
#else
#ifdef USE_OPENMP_THREADS
#pragma omp barrier
#endif
#endif
}

/*===========================================================================*/
/*---CUDA thread information---*/

/*---------------------------------------------------------------------------*/

static int Sweeper_nthreadblock( const Sweeper* sweeper, int axis )
{
  Assert( axis >= 0 && axis < 2 );

  return axis==0 ? sweeper->nthread_e : 1;
}

/*---------------------------------------------------------------------------*/

static int Sweeper_nthread_in_threadblock( const Sweeper* sweeper, int axis )
{
  Assert( axis >= 0 && axis < 2 );

  return axis==0 ? sweeper->nthread_octant : 1;
}

/*===========================================================================*/
/*---Number of elements of vslocal---*/

TARGET_HD static inline int Sweeper_nvslocal__( Sweeper* sweeper,
                                                Env*     env )
{
  return Env_cuda_is_using_device( env )
      ?
         sweeper->dims_b.na *
         NU *
         Sweeper_nthread_in_threadblock( sweeper, 0 ) *
         Sweeper_nthread_in_threadblock( sweeper, 1 ) *
         Sweeper_nthread_in_threadblock( sweeper, 2 )
       :
         sweeper->dims_b.na *
         NU *
         sweeper->nthread_e *
         sweeper->nthread_octant
       ;
}

/*===========================================================================*/
/*---Select which part of vslocal to use for current thread---*/

TARGET_HD static inline P* __restrict__ Sweeper_vslocal_this__(
                                                            Sweeper* sweeper )
{
#ifdef __CUDA_ARCH__
  return ( (P*) Env_cuda_shared_memory() )
    + sweeper->dims_b.na * NU *
      Sweeper_thread_octant( sweeper )
  ;
#else
  return sweeper->vslocal
    + sweeper->dims_b.na * NU *
      Sweeper_thread( sweeper )
  ;
#endif
}

/*===========================================================================*/
/*---Full size of vslocal---*/

static int Sweeper_shared_size__( Sweeper* sweeper,
                                Env* env )
{
  return Sweeper_nvslocal__( sweeper, env ) * sizeof( P );
}

/*===========================================================================*/
/*---Perform a sweep for a semiblock---*/

TARGET_HD void Sweeper_sweep_semiblock(
  Sweeper*               sweeper,
  P* __restrict__        vo,
  const P* __restrict__  vi,
  P* __restrict__        facexy,
  P* __restrict__        facexz,
  P* __restrict__        faceyz,
  const P* __restrict__  a_from_m,
  const P* __restrict__  m_from_a,
  const Quantities*      quan,
  const Step_Info        step_info,
  const int              octant_in_block,
  const int              ixmin,
  const int              ixmax,
  const int              iymin,
  const int              iymax,
  const int              izmin,
  const int              izmax,
  const Bool_t           do_block_init_this );

/*===========================================================================*/
/*---Perform a sweep for a block, implementation---*/

TARGET_HD void Sweeper_sweep_block_impl(
  Sweeper*               sweeper,
  P* __restrict__        vo,
  const P* __restrict__  vi,
  P* __restrict__        facexy,
  P* __restrict__        facexz,
  P* __restrict__        faceyz,
  const P* __restrict__  a_from_m,
  const P* __restrict__  m_from_a,
  int                    step,
  const Quantities*      quan,
  Bool_t                 proc_x_min,
  Bool_t                 proc_x_max,
  Bool_t                 proc_y_min,
  Bool_t                 proc_y_max,
  Step_Info_Values       step_info_values,
  unsigned long int      do_block_init );

/*===========================================================================*/
/*---Perform a sweep for a block, implementation, global---*/

TARGET_G void Sweeper_sweep_block_impl_global(
  Sweeper                sweeper,
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
  unsigned long int      do_block_init );

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
  Env*                   env );

/*===========================================================================*/
/*---Perform a sweep---*/

void Sweeper_sweep(
  Sweeper*               sweeper,
  Pointer*               vo,
  Pointer*               vi,
  const Quantities*      quan,
  Env*                   env );

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

#endif /*---_sweeper_kba_h_---*/

/*---------------------------------------------------------------------------*/
