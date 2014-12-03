/*---------------------------------------------------------------------------*/
/*!
 * \file   sweeper_kba_c_kernels.h
 * \author Wayne Joubert
 * \date   Tue Jan 28 16:37:41 EST 2014
 * \brief  sweeper_kba_c, code for device kernels.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _serial_c__sweeper_kba_c_kernels_h_
#define _serial_c__sweeper_kba_c_kernels_h_

#include "function_attributes.h"
#include "env_kernels.h"
#include "definitions_kernels.h"
#include "quantities_kernels.h"
#include "array_accessors_kernels.h"
#include "step_scheduler_kba_kernels.h"
#include "sweeper_kba_kernels.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*===========================================================================*/
/*---Apply boundary condition: xy face---*/

TARGET_HD static void Sweeper_set_boundary_xy(
  const Sweeper_Lite*   sweeper,
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
  const int iz_g = dir_z == DIR_UP ? -1 : sweeper->dims_g.nz;

  const int ie_min = ( sweeper->dims.ne *
                       ( Sweeper_thread_e( sweeper )     ) )
                   / sweeper->nthread_e;
  const int ie_max = ( sweeper->dims.ne *
                       ( Sweeper_thread_e( sweeper ) + 1 ) )
                   / sweeper->nthread_e;

  int ie = 0;

  for( ie=ie_min; ie<ie_max; ++ie )
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

TARGET_HD static void Sweeper_set_boundary_xz(
  const Sweeper_Lite*   sweeper,
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
  const int iy_g = dir_y == DIR_UP ? -1 : sweeper->dims_g.ny;

  const int ie_min = ( sweeper->dims.ne *
                       ( Sweeper_thread_e( sweeper )     ) )
                   / sweeper->nthread_e;
  const int ie_max = ( sweeper->dims.ne *
                       ( Sweeper_thread_e( sweeper ) + 1 ) )
                   / sweeper->nthread_e;

  int ie = 0;

  for( ie=ie_min; ie<ie_max; ++ie )
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

TARGET_HD static void Sweeper_set_boundary_yz(
  const Sweeper_Lite*   sweeper,
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
  const int ix_g = dir_x == DIR_UP ? -1 : sweeper->dims_g.nx;

  const int ie_min = ( sweeper->dims.ne *
                       ( Sweeper_thread_e( sweeper )     ) )
                   / sweeper->nthread_e;
  const int ie_max = ( sweeper->dims.ne *
                       ( Sweeper_thread_e( sweeper ) + 1 ) )
                   / sweeper->nthread_e;

  int ie = 0;

  for( ie=ie_min; ie<ie_max; ++ie )
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
/*---Perform a sweep for a cell---*/

TARGET_HD inline void Sweeper_sweep_cell(
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
  const int                      ix,
  const int                      iy,
  const int                      iz,
  const Bool_t                   do_block_init_this,
  const Bool_t                   is_cell_active )
{
  enum{ NU_PER_THREAD = NU / NTHREAD_U };

  int ia_base = 0;

  const int sweeper_thread_a = Sweeper_thread_a( sweeper );
  const int sweeper_thread_m = Sweeper_thread_m( sweeper );
  const int sweeper_thread_u = Sweeper_thread_u( sweeper );

  /*====================*/
  /*---Master loop over angle blocks---*/
  /*====================*/

  for( ia_base=0; ia_base<sweeper->dims_b.na; ia_base += NTHREAD_A )
  {
    int im_base = 0;

    for( im_base=0; im_base<NM; im_base += NTHREAD_M )
    {

      /*====================*/
      /*---If needed: reassign threads from angles to moments---*/
      /*====================*/

      if( im_base != 0 )
      {
        Sweeper_sync_amu_threads( sweeper );
      }

      /*====================*/
      /*---Load portion of vi---*/
      /*====================*/

      {
        __assume_aligned( vi_this, ( VEC_LEN < NTHREAD_M*NTHREAD_U ?
                                     VEC_LEN : NTHREAD_M*NTHREAD_U )
                                                                 * sizeof(P) );
        __assume_aligned( vilocal, ( VEC_LEN < NTHREAD_M*NTHREAD_U ?
                                     VEC_LEN : NTHREAD_M*NTHREAD_U )
                                                                 * sizeof(P) );
#ifndef __CUDA_ARCH__
        int sweeper_thread_m = 0;
        int sweeper_thread_u = 0;
        for( sweeper_thread_u=0; sweeper_thread_u<NTHREAD_U;
                                                          ++sweeper_thread_u )
#pragma ivdep
#pragma simd assert, vectorlengthfor( P )
        for( sweeper_thread_m=0; sweeper_thread_m<NTHREAD_M;
                                                          ++sweeper_thread_m )
#endif
        {
          const int im = im_base + sweeper_thread_m;
          if( ( NM % NTHREAD_M == 0 || im < NM ) &&
              is_cell_active )
          {
            if( ia_base == 0 ||
                NM*1 > NTHREAD_M*1 )
            {
              int iu_base = 0;
#pragma unroll
              for( iu_base=0; iu_base<NU; iu_base += NTHREAD_U )
              {
                const int iu = iu_base + sweeper_thread_u;

                if( NU % NTHREAD_U == 0 || iu < NU )
                {
                  *ref_vilocal( vilocal, sweeper->dims_b, NU, NTHREAD_M,
                                            sweeper_thread_m, iu ) =
                  *const_ref_state_flat( vi_this,
                                         sweeper->dims_b.nx,
                                         sweeper->dims_b.ny,
                                         sweeper->dims_b.nz,
                                         sweeper->dims_b.ne,
                                         NM,
                                         NU,
                                         ix, iy, iz, ie, im, iu );
                  /*---Can use this for non-MIC case:
                  /*--- *const_ref_state( vi_this, sweeper->dims_b, NU,
                  /*---                   ix, iy, iz, ie, im, iu );
                  ---*/
                }
              } /*---for iu---*/
            }
          }
        }
      }

      /*====================*/
      /*---Reassign threads from moments to angles---*/
      /*====================*/

      Sweeper_sync_amu_threads( sweeper );

      /*====================*/
      /*---Transform moments to angles---*/
      /*====================*/

      {
        /*
        WARNING!!!
        __assume( sweeper->dims_b.na % NTHREAD_A == 0 );
        */
        __assume_aligned( vslocal,  VEC_LEN * sizeof(P) );
        __assume_aligned( a_from_m, VEC_LEN * sizeof(P) );
        __assume_aligned( vilocal,  ( VEC_LEN < NTHREAD_M*NTHREAD_U ?
                                      VEC_LEN : NTHREAD_M*NTHREAD_U )
                                                                 * sizeof(P) );
#ifndef __CUDA_ARCH__
        int sweeper_thread_a = 0;
#pragma ivdep
#pragma simd assert, vectorlengthfor( P )
        for( sweeper_thread_a=0; sweeper_thread_a<NTHREAD_A;
                                                           ++sweeper_thread_a )
#endif
        {
          const int ia = ia_base + sweeper_thread_a;
          if( ia < sweeper->dims_b.na && is_cell_active )
          {
            int im_in_block = 0;
            int iu = 0;

            P v[NU];

#pragma unroll
            for( iu=0; iu<NU; ++iu )
            {
              v[iu] = ((P)0);
            }

            /*--------------------*/
            /*---Compute matvec in registers---*/
            /*--------------------*/

            for( im_in_block=0; im_in_block<NTHREAD_M; ++im_in_block )
            {
              const int im = im_base + im_in_block;

              if( NM % NTHREAD_M == 0 || im < NM )
              {
                const P a_from_m_this = *const_ref_a_from_m_flat(
                                             a_from_m,
                                             NM,
                                             sweeper->dims_b.na,
                                             im, ia, octant );
#pragma unroll
                for( iu=0; iu<NU; ++iu )
                {
                  v[iu] += a_from_m_this
                    * *const_ref_vilocal( vilocal, sweeper->dims_b,
                                        NU, NTHREAD_M, im_in_block, iu );

                }
              }
            } /*---for im_in_block---*/

            /*--------------------*/
            /*---Store/update to shared memory---*/
            /*--------------------*/

            if( im_base == 0 )
            {
#pragma unroll
              for( iu=0; iu<NU; ++iu )
              {
                vslocal[ ind_vslocal( sweeper->dims_b, NU, NTHREAD_A,
                                        sweeper_thread_a, iu ) ] = v[iu];
                /*---Can use this for non-MIC case
                *ref_vslocal( vslocal, sweeper->dims_b, NU, NTHREAD_A,
                                        sweeper_thread_a, iu )  = v[iu];
                */
              }
            }
            else
            {
#pragma unroll
              for( iu=0; iu<NU; ++iu )
              {
                vslocal[ ind_vslocal( sweeper->dims_b, NU, NTHREAD_A,
                                        sweeper_thread_a, iu ) ] += v[iu];
              }
            }
          } /*---if ia---*/
        }
      }
    } /*---for im_base---*/

    /*====================*/
    /*---Perform solve---*/
    /*====================*/

    /*
    WARNING!!!
    __assume( sweeper->dims_b.na % NTHREAD_A == 0 );
    */
    __assume_aligned( vslocal, VEC_LEN * sizeof(P) );
    __assume_aligned( facexy,  VEC_LEN * sizeof(P) );
    __assume_aligned( facexz,  VEC_LEN * sizeof(P) );
    __assume_aligned( faceyz,  VEC_LEN * sizeof(P) );
#ifndef __CUDA_ARCH__
    int sweeper_thread_a = 0;
#pragma ivdep
#pragma simd assert, vectorlengthfor( P )
    for( sweeper_thread_a=0; sweeper_thread_a<NTHREAD_A; ++sweeper_thread_a )
#endif
    {
      const int ia = ia_base + sweeper_thread_a;
      Quantities_solve( quan, vslocal,
                        ia, sweeper_thread_a, NTHREAD_A,
                        facexy, facexz, faceyz,
                        ix, iy, iz, ie,
                        ix+quan->ix_base, iy+quan->iy_base, iz+iz_base,
                        octant, octant_in_block,
                        sweeper->noctant_per_block,
                        sweeper->dims_b, sweeper->dims_g,
                        is_cell_active );
    }

    /*====================*/
    /*---Reassign threads from angles to moments---*/
    /*====================*/

    Sweeper_sync_amu_threads( sweeper );

    /*====================*/
    for( im_base=0; im_base<NM; im_base += NTHREAD_M )
    {
      {
        /*
        WARNING!!!
        __assume( sweeper->dims_b.na % NTHREAD_A == 0 );
        */
        __assume( sweeper->dims_b.nm == NM );
        __assume_aligned( vslocal,  VEC_LEN * sizeof(P) );
        __assume_aligned( m_from_a, VEC_LEN * sizeof(P) );
        __assume_aligned( vi_this,  ( VEC_LEN < NTHREAD_M*NTHREAD_U ?
                                      VEC_LEN : NTHREAD_M*NTHREAD_U )
                                                                 * sizeof(P) );
        __assume_aligned( vilocal,  ( VEC_LEN < NTHREAD_M*NTHREAD_U ?
                                      VEC_LEN : NTHREAD_M*NTHREAD_U )
                                                                 * sizeof(P) );
#ifndef __CUDA_ARCH__
        int sweeper_thread_u = 0;
        int sweeper_thread_m = 0;
        for( sweeper_thread_u=0; sweeper_thread_u<NTHREAD_U;
                                                           ++sweeper_thread_u )
#pragma ivdep
#pragma simd assert, vectorlengthfor( P )
        for( sweeper_thread_m=0; sweeper_thread_m<NTHREAD_M;
                                                           ++sweeper_thread_m )
#endif
        {
          const int im = im_base + sweeper_thread_m;

          P w[NU_PER_THREAD];

          int iu_per_thread = 0;
#pragma unroll
          for( iu_per_thread=0; iu_per_thread<NU_PER_THREAD; ++iu_per_thread )
          {
            w[iu_per_thread] = ((P)0);
          }

          /*====================*/
          /*---Transform angles to moments---*/
          /*====================*/

          if( ( NM % NTHREAD_M == 0 || im < NM ) &&
              is_cell_active )
          {
            int ia_in_block = 0;

            /*--------------------*/
            /*---Compute matvec in registers---*/
            /*--------------------*/

            /*---TODO: set up logic here to run fast for all cases---*/

#ifdef __MIC__
            if( ia_base + NTHREAD_A == sweeper->dims_b.na )
#else
            if( Bool_false )
#endif
            {
#ifdef __MIC__
/* "If applied to outer loop nests, the current implementation supports complete outer loop unrolling." */
#pragma unroll
#else
#pragma unroll 4
#endif
              for( ia_in_block=0; ia_in_block<NTHREAD_A; ++ia_in_block )
              {
                const int ia = ia_base + ia_in_block;

                const P m_from_a_this = m_from_a[
                                    ind_m_from_a_flat( sweeper->dims_b.nm,
                                                       sweeper->dims_b.na,
                                                       im, ia, octant ) ];
                {
#pragma unroll
                  for( iu_per_thread=0; iu_per_thread<NU_PER_THREAD;
                                                              ++iu_per_thread )
                  {
                    const int iu =  sweeper_thread_u + NTHREAD_U *
                                    iu_per_thread;

                    if( NU % NTHREAD_U == 0 || iu < NU )
                    {
                      w[ iu_per_thread ] +=
                          m_from_a_this
                          /*---Can use this for non-MIC case:
                          /*--- *const_ref_m_from_a( m_from_a, sweeper->dims_b,
                          /*---                     im, ia, octant )
                          ---*/
                        * *const_ref_vslocal( vslocal, sweeper->dims_b, NU,
                                              NTHREAD_A, ia_in_block, iu );
                    }
                  } /*---for iu_per_thread---*/
                }
              } /*---for ia_in_block---*/
            }
            else /*---ia_base---*/
            {
#ifdef __MIC__
/* "If applied to outer loop nests, the current implementation supports complete outer loop unrolling." */
#pragma unroll
#else
#pragma unroll 4
#endif
              for( ia_in_block=0; ia_in_block<NTHREAD_A; ++ia_in_block )
              {
                const int ia = ia_base + ia_in_block;
                const Bool_t mask = ia < sweeper->dims_b.na;

                const P m_from_a_this = mask ? m_from_a[
                                    ind_m_from_a_flat( sweeper->dims_b.nm,
                                                       sweeper->dims_b.na,
                                                       im, ia, octant ) ]
                                    : ((P)0);
                {
#pragma unroll
                  for( iu_per_thread=0; iu_per_thread<NU_PER_THREAD;
                                                              ++iu_per_thread )
                  {
                    const int iu =  sweeper_thread_u + NTHREAD_U *
                                    iu_per_thread;

                    if( NU % NTHREAD_U == 0 || iu < NU )
                    {
                      w[ iu_per_thread ] += mask ?
                          m_from_a_this
                          /*---Can use this for non-MIC case:
                          /*--- *const_ref_m_from_a( m_from_a, sweeper->dims_b,
                          /*---                     im, ia, octant )
                          ---*/
                        * *const_ref_vslocal( vslocal, sweeper->dims_b, NU,
                                              NTHREAD_A, ia_in_block, iu )
                        : ((P)0);
                    }
                  } /*---for iu_per_thread---*/
                }
              } /*---for ia_in_block---*/
            } /*---if ia_base---*/

            /*--------------------*/
            /*---Store/update to shared memory---*/
            /*--------------------*/

            if( ia_base == 0 ||
                NM*1 > NTHREAD_M*1 )
            {
#pragma unroll
              for( iu_per_thread=0; iu_per_thread<NU_PER_THREAD;
                                                              ++iu_per_thread )
              {
                const int iu =  sweeper_thread_u + NTHREAD_U *
                                iu_per_thread;

                if( NU % NTHREAD_U == 0 || iu < NU )
                {
                  *ref_volocal( volocal, sweeper->dims_b, NU, NTHREAD_M,
                            sweeper_thread_m, iu )  = w[ iu_per_thread ];
                }
              } /*---for iu_per_thread---*/
            }
            else
            {
#pragma unroll
              for( iu_per_thread=0; iu_per_thread<NU_PER_THREAD;
                                                              ++iu_per_thread )
              {
                const int iu =  sweeper_thread_u + NTHREAD_U *
                                iu_per_thread;

                if( (NU*1) % (NTHREAD_U*1) == 0 || iu < NU*1 )
                {
                  *ref_volocal( volocal, sweeper->dims_b, NU, NTHREAD_M,
                            sweeper_thread_m, iu ) += w[ iu_per_thread ];
                }
              } /*---for iu_per_thread---*/
            }
          } /*---if im---*/

          /*====================*/
          /*---Store/update portion of vo---*/
          /*====================*/

          if( ( (NM*1) % (NTHREAD_M*1) == 0 || im < NM*1 ) &&
              is_cell_active )
          {
            if( ia_base+NTHREAD_A >= sweeper->dims_b.na ||
                NM*1 > NTHREAD_M*1 )
            {
              int iu_base = 0;
#ifdef USE_OPENMP_VO_ATOMIC
#pragma unroll
              for( iu_base=0; iu_base<NU; iu += NTHREAD_U )
              {
                const int iu = iu_base + sweeper_thread_u;

                if( (NU*1) % (NTHREAD_U*1) == 0 || iu < (NU*1) )
                {
#pragma omp atomic update
                  *ref_state( vo_this, sweeper->dims_b, NU,
                              ix, iy, iz, ie, im, iu ) +=
                    *ref_volocal( volocal, sweeper->dims_b, NU, NTHREAD_M,
                                  sweeper_thread_m, iu );
                }
              }
#else /*---USE_OPENMP_VO_ATOMIC---*/
              if( ( ! do_block_init_this ) ||
                  ( NM*1 > NTHREAD_M*1 && ! ia_base==0 ) )
              {
#pragma unroll
                for( iu_base=0; iu_base<NU; iu_base += NTHREAD_U )
                {
                  const int iu = iu_base + sweeper_thread_u;

                  if( (NU*1) % (NTHREAD_U*1) == 0 || iu < NU*1 )
                  {
                    /*---Can use this for non-MIC case:
                    /*--- *ref_state( vo_this, sweeper->dims_b, NU,
                    /*---             ix, iy, iz, ie, im, iu ) +=
                    */
                    *ref_state_flat( vo_this,
                                     sweeper->dims_b.nx,
                                     sweeper->dims_b.ny,
                                     sweeper->dims_b.nz,
                                     sweeper->dims_b.ne,
                                     NM,
                                     NU,
                                     ix, iy, iz, ie, im, iu ) +=
                    *ref_volocal( volocal, sweeper->dims_b, NU, NTHREAD_M,
                                    sweeper_thread_m, iu );
                  }
                }
              }
              else
              {
#pragma unroll
                for( iu_base=0; iu_base<NU; iu_base += NTHREAD_U )
                {
                  const int iu = iu_base + sweeper_thread_u;

                  if( (NU*1) % (NTHREAD_U*1) == 0 || iu < NU*1 )
                  {
                    /*---Can use this for non-MIC case:
                    /*--- *ref_state( vo_this, sweeper->dims_b, NU,
                    /*---             ix, iy, iz, ie, im, iu ) =
                    */
                    *ref_state_flat( vo_this,
                                     sweeper->dims_b.nx,
                                     sweeper->dims_b.ny,
                                     sweeper->dims_b.nz,
                                     sweeper->dims_b.ne,
                                     NM,
                                     NU,
                                     ix, iy, iz, ie, im, iu )  =
                      *ref_volocal( volocal, sweeper->dims_b, NU, NTHREAD_M,
                                sweeper_thread_m, iu );
                  }
                }
              }
#endif /*---USE_OPENMP_VO_ATOMIC---*/
            }
          } /*---if im---*/
        }
      }

    } /*---for im_base---*/

  } /*---for ia_base---*/

}

/*===========================================================================*/
/*---Perform a sweep for a semiblock---*/

TARGET_HD inline void Sweeper_sweep_semiblock(
  Sweeper_Lite*          sweeper,
  P* __restrict__        vo_this,
  const P* __restrict__  vi_this,
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
  const Bool_t           do_block_init_this )
{
  /*---Declarations---*/

  const int octant  = step_info.octant;
  const int iz_base = step_info.block_z * sweeper->dims_b.nz;

  const int dir_x = Dir_x( octant );
  const int dir_y = Dir_y( octant );
  const int dir_z = Dir_z( octant );

  const int dir_inc_x = Dir_inc(dir_x);
  const int dir_inc_y = Dir_inc(dir_y);
  const int dir_inc_z = Dir_inc(dir_z);

  /*---Calculate v*local part to use---*/

  P* __restrict__ vilocal = Sweeper_vilocal_this__( sweeper );
  P* __restrict__ vslocal = Sweeper_vslocal_this__( sweeper );
  P* __restrict__ volocal = Sweeper_volocal_this__( sweeper );

  /*---Calculate loop extents---*/

  const int ixbeg = dir_x==DIR_UP ? ixmin : ixmax;
  const int iybeg = dir_y==DIR_UP ? iymin : iymax;
  const int izbeg = dir_z==DIR_UP ? izmin : izmax;

  const int ixend = dir_x==DIR_DN ? ixmin : ixmax;
  const int iyend = dir_y==DIR_DN ? iymin : iymax;
  const int izend = dir_z==DIR_DN ? izmin : izmax;

  const int ie_min = ( sweeper->dims.ne *
                       ( Sweeper_thread_e( sweeper )     ) )
                   / sweeper->nthread_e;
  const int ie_max = ( sweeper->dims.ne *
                       ( Sweeper_thread_e( sweeper ) + 1 ) )
                   / sweeper->nthread_e;

  int ie = 0;

  /*--------------------*/
  /*---Loop over energy groups---*/
  /*--------------------*/

  for( ie=ie_min; ie<ie_max; ++ie )
  {
    int ix = 0;
    int iy = 0;
    int iz = 0;

    /*--------------------*/
    /*---Loop over gridcells, in proper direction---*/
    /*--------------------*/

    for( iz=izbeg; iz!=izend+dir_inc_z; iz+=dir_inc_z )
    {
    for( iy=iybeg; iy!=iyend+dir_inc_y; iy+=dir_inc_y )
    {
    for( ix=ixbeg; ix!=ixend+dir_inc_x; ix+=dir_inc_x )
    {
      const Bool_t is_cell_active = ix < sweeper->dims_b.nx &&
                                    iy < sweeper->dims_b.ny &&
                                    iz < sweeper->dims_b.nz;
      /*--------------------*/
      /*---Sweep cell---*/
      /*--------------------*/
      Sweeper_sweep_cell( sweeper, vo_this, vi_this, vilocal, vslocal, volocal,
                          facexy, facexz, faceyz, a_from_m, m_from_a, quan,
                          octant, iz_base, octant_in_block, ie, ix, iy, iz,
                          do_block_init_this, is_cell_active );
    }
    }
    } /*---ix/iy/iz---*/

  } /*---ie---*/
}

/*===========================================================================*/
/*---Perform a sweep for a block, implementation---*/

TARGET_HD void Sweeper_sweep_block_impl(
  Sweeper_Lite*          sweeper,
        P* __restrict__  vo,
  const P* __restrict__  vi,
        P* __restrict__  facexy,
        P* __restrict__  facexz,
        P* __restrict__  faceyz,
  const P* __restrict__  a_from_m,
  const P* __restrict__  m_from_a,
  int                    step,
  const Quantities*      quan,
  Bool_t                 proc_x_min,
  Bool_t                 proc_x_max,
  Bool_t                 proc_y_min,
  Bool_t                 proc_y_max,
  Step_Info_Values       step_info_values,
  unsigned long int      do_block_init )
{
  /*---Declarations---*/

    const int noctant_per_block = sweeper->noctant_per_block;

    int semiblock = 0;

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
    =      making the update of vo at the end of Sweeper_sweep_semiblock atomic.
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
/*
    for( semiblock=0; semiblock<1; ++semiblock )
*/
    {
      /*--------------------*/
      /*---Loop over octants in octant block---*/
      /*--------------------*/

        const int octant_in_block_min =
                             ( sweeper->noctant_per_block *
                               ( Sweeper_thread_octant( sweeper )     ) )
                           / sweeper->nthread_octant;
        const int octant_in_block_max =
                             ( sweeper->noctant_per_block *
                               ( Sweeper_thread_octant( sweeper ) + 1 ) )
                           / sweeper->nthread_octant;

        int octant_in_block = 0;

      for( octant_in_block=octant_in_block_min;
           octant_in_block<octant_in_block_max; ++octant_in_block )
      {
        /*---Get step info---*/

        const Step_Info step_info = step_info_values.step_info[octant_in_block];

        /*--------------------*/
        /*---Begin compute section---*/
        /*--------------------*/

        if( step_info.is_active )
        {
          const int iz_base = step_info.block_z * sweeper->dims_b.nz;

          const P* vi_this
            = const_ref_state( vi, sweeper->dims, NU, 0, 0, iz_base, 0, 0, 0 );
          P* vo_this
                  = ref_state( vo, sweeper->dims, NU, 0, 0, iz_base, 0, 0, 0 );

          const int dir_x = Dir_x( step_info.octant );
          const int dir_y = Dir_y( step_info.octant );
          const int dir_z = Dir_z( step_info.octant );

          const int do_block_init_this = !! ( do_block_init &
                           ( ((unsigned long int)1) <<
                             ( octant_in_block + noctant_per_block *
                               semiblock ) ) );

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
                                           ( dir_x == DIR_UP );

          const Bool_t has_x_lo =     is_semiblock_x_lo   || ! is_x_semiblocked;
          const Bool_t has_x_hi = ( ! is_semiblock_x_lo ) || ! is_x_semiblocked;

          const int ixmin_b =     has_x_lo ? 0
                                             : (sweeper->dims_b.nx+1) / 2;
          const int ixmax_b = ( ! has_x_hi ) ? (sweeper->dims_b.nx+1) / 2 - 1
                                             :  sweeper->dims_b.nx        - 1;

          const int ixmax_b_up2 = ( sweeper->dims_b.nx % 2 && ixmax_b==sweeper->dims_b.nx-1 ) ? ( ixmax_b + 1 ) : ixmax_b;

          /*--------------------*/

          const Bool_t is_y_semiblocked = sweeper->nsemiblock > (1<<1);
          const Bool_t is_semiblock_y_lo = ( ( semiblock & (1<<1) ) == 0 ) ==
                                           ( dir_y == DIR_UP );

          const Bool_t has_y_lo =     is_semiblock_y_lo   || ! is_y_semiblocked;
          const Bool_t has_y_hi = ( ! is_semiblock_y_lo ) || ! is_y_semiblocked;

          const int iymin_b =     has_y_lo ? 0
                                             : (sweeper->dims_b.ny+1) / 2;
          const int iymax_b = ( ! has_y_hi ) ? (sweeper->dims_b.ny+1) / 2 - 1
                                             :  sweeper->dims_b.ny        - 1;

          const int iymax_b_up2 = ( sweeper->dims_b.ny % 2 && iymax_b==sweeper->dims_b.ny-1 ) ? ( iymax_b + 1 ) : iymax_b;

          /*--------------------*/

          const Bool_t is_z_semiblocked = sweeper->nsemiblock > (1<<2);
          const Bool_t is_semiblock_z_lo = ( ( semiblock & (1<<2) ) == 0 ) ==
                                           ( dir_z == DIR_UP );

          const Bool_t has_z_lo =     is_semiblock_z_lo   || ! is_z_semiblocked;
          const Bool_t has_z_hi = ( ! is_semiblock_z_lo ) || ! is_z_semiblocked;

          const int izmin_b =     has_z_lo ? 0
                                             : (sweeper->dims_b.nz+1) / 2;
          const int izmax_b = ( ! has_z_hi ) ? (sweeper->dims_b.nz+1) / 2 - 1
                                             :  sweeper->dims_b.nz        - 1;

          const int izmax_b_up2 = ( sweeper->dims_b.nz % 2 && izmax_b==sweeper->dims_b.nz-1 ) ? ( izmax_b + 1 ) : izmax_b;

          /*--------------------*/
          /*---Set physical boundary conditions if part of semiblock---*/
          /*--------------------*/

          if( ( dir_z == DIR_UP && step_info.block_z == 0 && has_z_lo ) ||
              ( dir_z == DIR_DN && step_info.block_z ==
                                      sweeper->nblock_z - 1 && has_z_hi ) )
          {
            Sweeper_set_boundary_xy( sweeper, facexy, quan,
                                     step_info.octant, octant_in_block,
                                     ixmin_b, ixmax_b, iymin_b, iymax_b );
          }

          /*--------------------*/

          if( ( dir_y == DIR_UP && proc_y_min && has_y_lo ) ||
              ( dir_y == DIR_DN && proc_y_max && has_y_hi ) )
          {
            Sweeper_set_boundary_xz( sweeper, facexz, quan, step_info.block_z,
                                     step_info.octant, octant_in_block,
                                     ixmin_b, ixmax_b, izmin_b, izmax_b );
          }

          /*--------------------*/

          if( ( dir_x == DIR_UP && proc_x_min && has_x_lo ) ||
              ( dir_x == DIR_DN && proc_x_max && has_x_hi ) )
          {
            Sweeper_set_boundary_yz( sweeper, faceyz, quan, step_info.block_z,
                                     step_info.octant, octant_in_block,
                                     iymin_b, iymax_b, izmin_b, izmax_b );
          }

          /*--------------------*/
          /*---Perform sweep on relevant block---*/
          /*--------------------*/

          Sweeper_sweep_semiblock( sweeper, vo_this, vi_this,
                                   facexy, facexz, faceyz,
                                   a_from_m, m_from_a,
                                   quan, step_info, octant_in_block,
                                   ixmin_b, ixmax_b_up2,
                                   iymin_b, iymax_b_up2,
                                   izmin_b, izmax_b_up2,
                                   do_block_init_this );

        }  /*---is_active---*/

      } /*---octant_in_block---*/

      /*---Sync between semiblock steps---*/

      Sweeper_sync_octant_threads( sweeper );

    } /*---semiblock---*/
}

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

#endif /*---_serial_c__sweeper_kba_c_kernels_h_---*/

/*---------------------------------------------------------------------------*/
