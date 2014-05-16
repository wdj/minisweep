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
#include "faces_kba.h"

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

/*---NOTE: these should NOT be accessed outside of the Sweeper pseudo-class---*/

enum{ NTHREAD_DEVICE_U = ( NU * NM * 1 <= WARP_SIZE * 1 ) ? NU :
                           NM - 0 == 16 && NU - 0 == 4    ?  2 :
                                                            NU };
enum{ NTHREAD_DEVICE_M = WARP_SIZE / NTHREAD_DEVICE_U };
enum{ NTHREAD_DEVICE_A = NTHREAD_DEVICE_U * NTHREAD_DEVICE_M };

#ifdef __CUDA_ARCH__
  enum{ NTHREAD_A = NTHREAD_DEVICE_A };
  enum{ NTHREAD_M = NTHREAD_DEVICE_M };
  enum{ NTHREAD_U = NTHREAD_DEVICE_U };
#else
  enum{ NTHREAD_A = 1 };
  enum{ NTHREAD_M = 1 };
  enum{ NTHREAD_U = 1 };
#endif

enum{ NU_PER_THREAD = NU / NTHREAD_U };

/*===========================================================================*/
/*---Struct with pointers etc. used to perform sweep---*/

typedef struct
{
  P* __restrict__  vilocal_host__;
  P* __restrict__  vslocal_host__;
  P* __restrict__  volocal_host__;

  Dimensions       dims;
  Dimensions       dims_b;
  Dimensions       dims_g;

  int              nthread_e;
  int              nthread_octant;

  int              nblock_z;
  int              nblock_octant;
  int              noctant_per_block;
  int              nsemiblock;

  Step_Scheduler   step_scheduler;

  Faces            faces;
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
  return Env_cuda_thread_in_threadblock( 1 );
#else
  Assert( sweeper->nthread_e * sweeper->nthread_octant == 1 ||
                                                       Env_omp_in_parallel() );
  return Env_omp_thread() / sweeper->nthread_e;
#endif
}

/*---------------------------------------------------------------------------*/

TARGET_HD static inline int Sweeper_thread_a( const Sweeper* sweeper )
{
#ifdef __CUDA_ARCH__
  return Env_cuda_thread_in_threadblock( 0 );
#else
  return 0;
#endif
}

/*---------------------------------------------------------------------------*/

TARGET_HD static inline int Sweeper_thread_m( const Sweeper* sweeper )
{
#ifdef __CUDA_ARCH__
  return Env_cuda_thread_in_threadblock( 0 ) / NTHREAD_U;
#else
  return 0;
#endif
}

/*---------------------------------------------------------------------------*/

TARGET_HD static inline int Sweeper_thread_u( const Sweeper* sweeper )
{
#ifdef __CUDA_ARCH__
  return Env_cuda_thread_in_threadblock( 0 ) % NTHREAD_U;
#else
  return 0;
#endif
}

/*===========================================================================*/
/*---Thread synchronization---*/

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

/*---------------------------------------------------------------------------*/

TARGET_HD static void Sweeper_sync_amu_threads( Sweeper* sweeper )
{
#ifdef __CUDA_ARCH__
  Env_cuda_sync_threadblock();
#endif
}

/*===========================================================================*/
/*---Number of elements to allocate for v*local---*/

TARGET_HD static inline int Sweeper_nvilocal__( Sweeper* sweeper,
                                                Env*     env )
{
  return Env_cuda_is_using_device( env )
      ?
         NTHREAD_DEVICE_M *
         NU *
         sweeper->nthread_octant
       :
         NTHREAD_M *
         NU *
         sweeper->nthread_octant *
         sweeper->nthread_e
       ;
}

/*---------------------------------------------------------------------------*/

TARGET_HD static inline int Sweeper_nvslocal__( Sweeper* sweeper,
                                                Env*     env )
{
  return Env_cuda_is_using_device( env )
      ?
         NTHREAD_DEVICE_A *
         NU *
         sweeper->nthread_octant
       :
         NTHREAD_A *
         NU *
         sweeper->nthread_octant *
         sweeper->nthread_e
       ;
}

/*---------------------------------------------------------------------------*/

TARGET_HD static inline int Sweeper_nvolocal__( Sweeper* sweeper,
                                                Env*     env )
{
  return Env_cuda_is_using_device( env )
      ?
         NTHREAD_DEVICE_M *
         NU *
         sweeper->nthread_octant
       :
         NTHREAD_M *
         NU *
         sweeper->nthread_octant *
         sweeper->nthread_e
       ;
}

/*===========================================================================*/
/*---Select which part of v*local to use for current thread/block---*/

TARGET_HD static inline P* __restrict__ Sweeper_vilocal_this__(
                                                            Sweeper* sweeper )
{
#ifdef __CUDA_ARCH__
  return ( (P*) Env_cuda_shared_memory() )
    + NTHREAD_DEVICE_M * NU *
      Sweeper_thread_octant( sweeper )
  ;
#else
  return sweeper->vilocal_host__
    + NTHREAD_M * NU *
      ( Sweeper_thread_octant( sweeper ) + sweeper->nthread_octant *
        Sweeper_thread_e(      sweeper ) )
  ;
#endif
}

/*---------------------------------------------------------------------------*/

TARGET_HD static inline P* __restrict__ Sweeper_vslocal_this__(
                                                            Sweeper* sweeper )
{
#ifdef __CUDA_ARCH__
  return ( (P*) Env_cuda_shared_memory() )
    + ( NTHREAD_DEVICE_M *
        NU *
        sweeper->nthread_octant ) * 2
    + NTHREAD_A * NU *
      Sweeper_thread_octant( sweeper )
  ;
#else
  return sweeper->vslocal_host__
    + NTHREAD_A * NU *
      ( Sweeper_thread_octant( sweeper ) + sweeper->nthread_octant *
        Sweeper_thread_e(      sweeper ) )
  ;
#endif
}

/*---------------------------------------------------------------------------*/

TARGET_HD static inline P* __restrict__ Sweeper_volocal_this__(
                                                            Sweeper* sweeper )
{
#ifdef __CUDA_ARCH__
  return ( (P*) Env_cuda_shared_memory() )
    + ( NTHREAD_DEVICE_M *
        NU *
        sweeper->nthread_octant ) * 1
    + NTHREAD_M * NU *
      Sweeper_thread_octant( sweeper )
  ;
#else
  return sweeper->volocal_host__
    + NTHREAD_M * NU *
      ( Sweeper_thread_octant( sweeper ) + sweeper->nthread_octant *
        Sweeper_thread_e(      sweeper ) )
  ;
#endif
}

/*===========================================================================*/
/*---For kernel launch: CUDA thread/block counts---*/

static int Sweeper_nthreadblock( const Sweeper* sweeper, int axis )
{
  Assert( axis >= 0 && axis < 2 );

  return axis==0 ? sweeper->nthread_e :
         axis==1 ? 1 :
                   1;
}

/*---------------------------------------------------------------------------*/

static int Sweeper_nthread_in_threadblock( const Sweeper* sweeper, int axis )
{
  Assert( axis >= 0 && axis < 2 );

  return axis==0 ? NTHREAD_DEVICE_A :
         axis==1 ? sweeper->nthread_octant :
                   1;
}

/*===========================================================================*/
/*---Fof kernel launch: full size of CUDA shared memory---*/

static int Sweeper_shared_size__( Sweeper* sweeper,
                                  Env*     env )
{
  return ( Sweeper_nvilocal__( sweeper, env ) +
           Sweeper_nvslocal__( sweeper, env ) +
           Sweeper_nvolocal__( sweeper, env ) ) * sizeof( P );
}

/*===========================================================================*/
/*---Perform a sweep for a cell---*/

TARGET_HD void Sweeper_sweep_cell(
  Sweeper*               sweeper,
  P* __restrict__        vo_this,
  const P* __restrict__  vi_this,
  P* __restrict__        vilocal,
  P* __restrict__        vslocal,
  P* __restrict__        volocal,
  P* __restrict__        facexy,
  P* __restrict__        facexz,
  P* __restrict__        faceyz,
  const P* __restrict__  a_from_m,
  const P* __restrict__  m_from_a,
  const Quantities*      quan,
  const int              octant,
  const int              iz_base,
  const int              octant_in_block,
  const int              ie,
  const int              ix,
  const int              iy,
  const int              iz,
  const Bool_t           do_block_init_this,
  const Bool_t           is_cell_active );

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
