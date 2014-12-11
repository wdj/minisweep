/*---------------------------------------------------------------------------*/
/*!
 * \file   sweeper_kba_kernels.h
 * \author Wayne Joubert
 * \date   Tue Jan 28 16:37:41 EST 2014
 * \brief  sweeper_kba, code for device kernels.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _sweeper_kba_kernels_h_
#define _sweeper_kba_kernels_h_

#include "function_attributes.h"
#include "env_kernels.h"
#include "definitions_kernels.h"
#include "dimensions_kernels.h"
#include "pointer_kernels.h"
#include "quantities_kernels.h"

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

enum{ NTHREAD_DEVICE_U = VEC_LEN <= NU                  ? VEC_LEN :
                         ( NU * NM * 1 <= VEC_LEN * 1 ) ? NU :
                           NM - 0 == 16 && NU - 0 == 4  ?  2 :
                                                          NU };
enum{ NTHREAD_DEVICE_M = VEC_LEN / NTHREAD_DEVICE_U };
enum{ NTHREAD_DEVICE_A = NTHREAD_DEVICE_U * NTHREAD_DEVICE_M };

#ifdef __CUDA_ARCH__
  enum{ NTHREAD_A = NTHREAD_DEVICE_A };
  enum{ NTHREAD_M = NTHREAD_DEVICE_M };
  enum{ NTHREAD_U = NTHREAD_DEVICE_U };
#else
#ifdef __MIC__
  enum{ NTHREAD_A = VEC_LEN * 4 }; /*---tuning parameter---*/
  enum{ NTHREAD_M = NM };
  enum{ NTHREAD_U = NU };
#else
  enum{ NTHREAD_A = NTHREAD_DEVICE_A };
  enum{ NTHREAD_M = NTHREAD_DEVICE_M };
  enum{ NTHREAD_U = NTHREAD_DEVICE_U };
#endif
#endif

/*===========================================================================*/
/*---Lightweight version of Sweeper class for sending to device---*/

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
  int              nthread_y;
  int              nthread_z;

  int              nblock_z;
  int              nblock_octant;
  int              noctant_per_block;
  int              nsemiblock;
  int              nsubblock_x;
  int              nsubblock_y;
  int              nsubblock_z;
} Sweeper_Lite;

/*===========================================================================*/
/*---Thread indexers---*/

TARGET_HD static inline int Sweeper_thread_e( const Sweeper_Lite* sweeper )
{
#ifdef __CUDA_ARCH__
  return Env_cuda_threadblock( 0 );
#else
  Assert( sweeper->nthread_e *
          sweeper->nthread_octant *
          sweeper->nthread_y *
          sweeper->nthread_z == 1 || Env_omp_in_parallel() );
  return Env_omp_thread() % sweeper->nthread_e;
#endif
}

/*---------------------------------------------------------------------------*/

TARGET_HD static inline int Sweeper_thread_octant( const Sweeper_Lite* sweeper )
{
#ifdef __CUDA_ARCH__
  return Env_cuda_thread_in_threadblock( 1 );
#else
  Assert( sweeper->nthread_e *
          sweeper->nthread_octant *
          sweeper->nthread_y *
          sweeper->nthread_z == 1 || Env_omp_in_parallel() );
  return ( Env_omp_thread() / sweeper->nthread_e )
                            % sweeper->nthread_octant;
#endif
}

/*---------------------------------------------------------------------------*/

TARGET_HD static inline int Sweeper_thread_y( const Sweeper_Lite* sweeper )
{
#ifdef __CUDA_ARCH__
  return Env_cuda_threadblock( 1 );
#else
  Assert( sweeper->nthread_e *
          sweeper->nthread_octant *
          sweeper->nthread_y *
          sweeper->nthread_z == 1 || Env_omp_in_parallel() );
  return ( Env_omp_thread() / ( sweeper->nthread_e *
                                sweeper->nthread_octant )
                            %   sweeper->nthread_y );
#endif
}

/*---------------------------------------------------------------------------*/

TARGET_HD static inline int Sweeper_thread_z( const Sweeper_Lite* sweeper )
{
#ifdef __CUDA_ARCH__
  return Env_cuda_threadblock( 2 );
#else
  Assert( sweeper->nthread_e *
          sweeper->nthread_octant *
          sweeper->nthread_y *
          sweeper->nthread_z == 1 || Env_omp_in_parallel() );
  return ( Env_omp_thread() / ( sweeper->nthread_e *
                                sweeper->nthread_octant *
                                sweeper->nthread_y )
                            %   sweeper->nthread_z );
#endif
}

/*---------------------------------------------------------------------------*/

TARGET_HD static inline int Sweeper_thread_a( const Sweeper_Lite* sweeper )
{
#ifdef __CUDA_ARCH__
  return Env_cuda_thread_in_threadblock( 0 );
#else
  return 0;
#endif
}

/*---------------------------------------------------------------------------*/

TARGET_HD static inline int Sweeper_thread_m( const Sweeper_Lite* sweeper )
{
#ifdef __CUDA_ARCH__
  return Env_cuda_thread_in_threadblock( 0 ) / NTHREAD_U;
#else
  return 0;
#endif
}

/*---------------------------------------------------------------------------*/

TARGET_HD static inline int Sweeper_thread_u( const Sweeper_Lite* sweeper )
{
#ifdef __CUDA_ARCH__
  return Env_cuda_thread_in_threadblock( 0 ) % NTHREAD_U;
#else
  return 0;
#endif
}

/*===========================================================================*/
/*---Thread synchronization---*/

TARGET_HD static inline void Sweeper_sync_octant_threads( Sweeper_Lite* sweeper )
{
#ifdef __CUDA_ARCH__
  /*---NOTE: this may not be needed if these threads are mapped in-warp---*/
  Env_cuda_sync_threadblock();
#else
#ifdef USE_OPENMP_THREADS
#pragma omp barrier
#endif
#endif
}

/*---------------------------------------------------------------------------*/

TARGET_HD static inline void Sweeper_sync_yz_threads( Sweeper_Lite* sweeper )
{
#ifdef __CUDA_ARCH__
  /*---NOTE: this may not be needed if these threads are mapped in-warp---*/
  Env_cuda_sync_threadblock();
#else
#ifdef USE_OPENMP_THREADS
#pragma omp barrier
#endif
#endif
}

/*---------------------------------------------------------------------------*/

TARGET_HD static inline void Sweeper_sync_amu_threads( Sweeper_Lite* sweeper )
{
#ifdef __CUDA_ARCH__
  Env_cuda_sync_threadblock();
#endif
}

/*===========================================================================*/
/*---Select which part of v*local to use for current thread/block---*/

TARGET_HD static inline P* __restrict__ Sweeper_vilocal_this__(
                                                            Sweeper_Lite* sweeper )
{
#ifdef __CUDA_ARCH__
  return ( (P*) Env_cuda_shared_memory() )
    + NTHREAD_M * NU *
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
                                                            Sweeper_Lite* sweeper )
{
#ifdef __CUDA_ARCH__
  return ( (P*) Env_cuda_shared_memory() )
    + ( NTHREAD_M *
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
                                                            Sweeper_Lite* sweeper )
{
#ifdef __CUDA_ARCH__
  return ( (P*) Env_cuda_shared_memory() )
    + ( NTHREAD_M *
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
/*---Perform a sweep for a cell---*/

TARGET_HD inline void Sweeper_sweep_cell(
  Sweeper_Lite* __restrict__     sweeper,
  P* __restrict__                vo_this,
  const P* const __restrict__    vi_this,
  P* __restrict__                vilocal,
  P* __restrict__                vslocal,
  P* __restrict__                volocal,
  P* __restrict__                facexy,
  P* __restrict__                facexz,
  P* __restrict__                faceyz,
  const P* const __restrict__    a_from_m,
  const P* const __restrict__    m_from_a,
  const Quantities* __restrict__ quan,
  const int                      octant,
  const int                      iz_base,
  const int                      octant_in_block,
  const int                      ie,
  const int                      ix,
  const int                      iy,
  const int                      iz,
  const Bool_t                   do_block_init_this,
  const Bool_t                   is_cell_active );

/*===========================================================================*/
/*---Perform a sweep for a subblock---*/

TARGET_HD inline void Sweeper_sweep_subblock(
  Sweeper_Lite* __restrict__     sweeper,
  P* const __restrict__          vo_this,
  const P* const __restrict__    vi_this,
  P* const __restrict__          vilocal,
  P* const __restrict__          vslocal,
  P* const __restrict__          volocal,
  P* const __restrict__          facexy,
  P* const __restrict__          facexz,
  P* const __restrict__          faceyz,
  const P* const __restrict__    a_from_m,
  const P* const __restrict__    m_from_a,
  const Quantities* __restrict__ quan,
  const int                      octant,
  const int                      iz_base,
  const int                      octant_in_block,
  const int                      ie,
  const int                      ixmin,
  const int                      ixmax,
  const int                      iymin,
  const int                      iymax,
  const int                      izmin,
  const int                      izmax,
  const Bool_t                   do_block_init_this );

/*===========================================================================*/
/*---Perform a sweep for a semiblock---*/

TARGET_HD inline void Sweeper_sweep_semiblock(
  Sweeper_Lite*          sweeper,
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
  Sweeper_Lite*          sweeper,
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

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

#endif /*---_sweeper_kba_kernels_h_---*/

/*---------------------------------------------------------------------------*/
