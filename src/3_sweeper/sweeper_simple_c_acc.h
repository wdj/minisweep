/*---------------------------------------------------------------------------*/
/*!
 * \file   sweeper_simple_c.h
 * \author Wayne Joubert
 * \date   Wed Jan 15 16:06:28 EST 2014
 * \brief  Definitions for performing a sweep, simple version.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _sweeper_simple_c_h_
#define _sweeper_simple_c_h_

#include "env.h"
#include "definitions.h"
#include "quantities.h"
#include "array_accessors.h"
#include "array_operations.h"
#include "sweeper_simple.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*===========================================================================*/
/*---Null object---*/

Sweeper Sweeper_null()
{
  Sweeper result;
  memset( (void*)&result, 0, sizeof(Sweeper) );
  return result;
}

/*===========================================================================*/
/*---Pseudo-constructor for Sweeper struct---*/

void Sweeper_create( Sweeper*          sweeper,
                     Dimensions        dims,
                     const Quantities* quan,
                     Env*              env,
                     Arguments*        args )
{
  Insist( Env_nproc( env ) == 1 && 
                             "This sweeper version runs only with one proc." );

  /*---Allocate arrays---*/

  sweeper->vslocal = malloc_host_P( dims.na * NU );
  sweeper->facexy  = malloc_host_P( dims.ncell_x * dims.ncell_y * dims.ne *
                         dims.na * NU * Sweeper_noctant_per_block( sweeper ) );
  sweeper->facexz  = malloc_host_P( dims.ncell_x * dims.ncell_z * dims.ne *
                         dims.na * NU * Sweeper_noctant_per_block( sweeper ) );
  sweeper->faceyz  = malloc_host_P( dims.ncell_y * dims.ncell_z * dims.ne *
                         dims.na * NU * Sweeper_noctant_per_block( sweeper ) );

  sweeper->dims = dims;
}

/*===========================================================================*/
/*---Pseudo-destructor for Sweeper struct---*/

void Sweeper_destroy( Sweeper* sweeper,
                      Env*     env )
{
  /*---Deallocate arrays---*/

  free_host_P( sweeper->vslocal );
  free_host_P( sweeper->facexy );
  free_host_P( sweeper->facexz );
  free_host_P( sweeper->faceyz );

  sweeper->vslocal = NULL;
  sweeper->facexy  = NULL;
  sweeper->facexz  = NULL;
  sweeper->faceyz  = NULL;
}

/*===========================================================================*/
/*---Quantities_init_face inline routine---*/

#pragma acc routine seq
P Quantities_init_face(int ia, int ie, int iu, int scalefactor_space, int octant)
{
  /*--- Quantities_init_facexy inline ---*/

  /*--- Quantities_affinefunction_ inline ---*/
  return ( (P) (1 + ia) ) 

    /*--- Quantities_scalefactor_angle_ inline ---*/
    * ( (P) (1 << (ia & ( (1<<3) - 1))) ) 

    /*--- Quantities_scalefactor_space_ inline ---*/
    * ( (P) scalefactor_space)

    /*--- Quantities_scalefactor_energy_ inline ---*/
    * ( (P) (1 << ((( (ie) * 1366 + 150889) % 714025) & ( (1<<2) - 1))) )

    /*--- Quantities_scalefactor_unknown_ inline ---*/
    * ( (P) (1 << ((( (iu) * 741 + 60037) % 312500) & ( (1<<2) - 1))) )

    /*--- Quantities_scalefactor_octant_ ---*/
    * ( (P) 1 + octant);

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
  int ix = 0;
  int iy = 0;
  int iz = 0;
  int ie = 0;
  int im = 0;
  int ia = 0;
  int iu = 0;
  int octant = 0;

  /*--- Dimensions ---*/
  int dims_ncell_x = sweeper->dims.ncell_x;
  int dims_ncell_y = sweeper->dims.ncell_y;
  int dims_ncell_z = sweeper->dims.ncell_z;
  int dims_ne = sweeper->dims.ne;
  int dims_nm = sweeper->dims.nm;
  int dims_na = sweeper->dims.na;
  Dimensions dims = sweeper->dims;

  /*---Initialize result array to zero---*/

  initialize_state_zero( Pointer_h( vo ), sweeper->dims, NU );

  /*---Loop over octants---*/

  for( octant=0; octant<NOCTANT; ++octant )
  {
    const int octant_in_block = 0;
    Assert( octant_in_block >= 0 &&
            octant_in_block < Sweeper_noctant_per_block( sweeper ) );

    /*---Decode octant directions from octant number---*/

    const int dir_x = Dir_x( octant );
    const int dir_y = Dir_y( octant );
    const int dir_z = Dir_z( octant );

    /*---Initialize faces---*/

    /*---The semantics of the face arrays are as follows.
         On entering a cell for a solve at the gridcell level,
         the face array is assumed to have a value corresponding to
         "one cell lower" in the relevant direction.
         On leaving the gridcell solve, the face has been updated
         to have the flux at that gridcell.
         Thus, the face is initialized at first to have a value
         "one cell" outside of the domain, e.g., for the XY face,
         either -1 or dims.ncell_x.
         Note also that the face initializer functions now take
         coordinates for all three spatial dimensions --
         the third dimension is used to denote whether it is the
         "lower" or "upper" face and also its exact location
         in that dimension.
    ---*/

    /*--- Array Pointers ---*/
    P* __restrict__ facexy = sweeper->facexy;
    P* __restrict__ facexz = sweeper->facexz;
    P* __restrict__ faceyz = sweeper->faceyz;
    P* v_a_from_m = (P*) Pointer_const_h( & quan->a_from_m);
    P* v_m_from_a = (P*) Pointer_const_h( & quan->m_from_a);
    P* vi_h = Pointer_h( vi );
    P* vo_h = Pointer_h( vo );
    P* vs_local = sweeper->vslocal;

    /*--- Array Sizes ---*/
    int facexy_size = dims_ncell_x * dims_ncell_y * 
      dims_ne * dims_na * NU * Sweeper_noctant_per_block( sweeper );
    int facexz_size = dims_ncell_x * dims_ncell_z * 
      dims_ne * dims_na * NU * Sweeper_noctant_per_block( sweeper );
    int faceyz_size = dims_ncell_y * dims_ncell_z * 
      dims_ne * dims_na * NU * Sweeper_noctant_per_block( sweeper );
      int v_size = dims.nm * dims.na * NOCTANT;
      int vi_h_size = dims.ncell_x * dims.ncell_y * dims.ncell_z * 
                 	dims.ne * dims.nm * NU;
      int vo_h_size = dims.ncell_x * dims.ncell_y * dims.ncell_z * 
                 	dims.ne * dims.nm * NU;
      int vs_local_size = dims.na * NU;

#pragma acc parallel copy(facexy[:facexy_size], facexz[:facexz_size], faceyz[:faceyz_size])
{
    {
      iz = dir_z == DIR_UP ? -1 : dims_ncell_z;

#pragma acc loop seq
      for( iu=0; iu<NU; ++iu )
#pragma acc loop seq
      for( iy=0; iy<dims_ncell_y; ++iy )
#pragma acc loop seq
      for( ix=0; ix<dims_ncell_x; ++ix )
#pragma acc loop independent gang
      for( ie=0; ie<dims_ne; ++ie )
#pragma acc loop independent vector
      for( ia=0; ia<dims_na; ++ia )
      {

	/*--- Quantities_scalefactor_space_ inline ---*/
	int scalefactor_space = 0;
	scalefactor_space = ( (scalefactor_space+(ix+2))*8121 + 28411 ) % 134456;
	scalefactor_space = ( (scalefactor_space+(iy+2))*8121 + 28411 ) % 134456;
	scalefactor_space = ( (scalefactor_space+(iz+2))*8121 + 28411 ) % 134456;
	scalefactor_space = ( (scalefactor_space+(ix+3*iy+7*iz+2))*8121 + 28411 ) % 134456;
	scalefactor_space = ix+3*iy+7*iz+2;
	scalefactor_space = scalefactor_space & ( (1<<2) - 1 );
	scalefactor_space = 1 << scalefactor_space;

	/*--- ref_facexy inline ---*/
	facexy[ia + dims_na      * (
			iu + NU           * (
                        ix + dims_ncell_x * (
                        iy + dims_ncell_y * (
                        ie + dims_ne      * (
                        octant_in_block ))))) ]

	  /*--- Quantities_init_face routine ---*/
	  = Quantities_init_face(ia, ie, iu, scalefactor_space, octant);

        /* *ref_facexy( sweeper->facexy, sweeper->dims, NU, */ 
	/*              Sweeper_noctant_per_block( sweeper ), */
        /*              ix, iy, ie, ia, iu, octant_in_block ) */
	/* = Quantities_init_facexy( */
	/*             quan, ix, iy, iz, ie, ia, iu, octant, sweeper->dims ); */
      }
      }

/*--- #pragma acc parallel ---*/
 }

    {
      iy = dir_y == DIR_UP ? -1 : sweeper->dims.ncell_y;

      for( iu=0; iu<NU; ++iu )
      for( iz=0; iz<sweeper->dims.ncell_z; ++iz )
      for( ix=0; ix<sweeper->dims.ncell_x; ++ix )
      for( ie=0; ie<sweeper->dims.ne; ++ie )
      for( ia=0; ia<sweeper->dims.na; ++ia )
      {
	/*--- Quantities_scalefactor_space_ inline ---*/
	int scalefactor_space = 0;
	scalefactor_space = ( (scalefactor_space+(ix+2))*8121 + 28411 ) % 134456;
	scalefactor_space = ( (scalefactor_space+(iy+2))*8121 + 28411 ) % 134456;
	scalefactor_space = ( (scalefactor_space+(iz+2))*8121 + 28411 ) % 134456;
	scalefactor_space = ( (scalefactor_space+(ix+3*iy+7*iz+2))*8121 + 28411 ) % 134456;
	scalefactor_space = ix+3*iy+7*iz+2;
	scalefactor_space = scalefactor_space & ( (1<<2) - 1 );
	scalefactor_space = 1 << scalefactor_space;

	/*--- ref_facexz inline ---*/
	facexz[ia + dims_na      * (
			iu + NU           * (
                        ix + dims_ncell_x * (
                        iz + dims_ncell_z * (
                        ie + dims_ne      * (
                        octant_in_block ))))) ]

	  /*--- Quantities_init_face routine ---*/
	  = Quantities_init_face(ia, ie, iu, scalefactor_space, octant);

        /* *ref_facexz( sweeper->facexz, sweeper->dims, NU, */
        /*              Sweeper_noctant_per_block( sweeper ), */
        /*              ix, iz, ie, ia, iu, octant_in_block ) */
        /*      = Quantities_init_facexz( */
        /*                  quan, ix, iy, iz, ie, ia, iu, octant, sweeper->dims ); */
      }
    }

    {
      ix = dir_x == DIR_UP ? -1 : sweeper->dims.ncell_x;
      for( iu=0; iu<NU; ++iu )
      for( iz=0; iz<sweeper->dims.ncell_z; ++iz )
      for( iy=0; iy<sweeper->dims.ncell_y; ++iy )
      for( ie=0; ie<sweeper->dims.ne; ++ie )
      for( ia=0; ia<sweeper->dims.na; ++ia )
      {
	/*--- Quantities_scalefactor_space_ inline ---*/
	int scalefactor_space = 0;
	scalefactor_space = ( (scalefactor_space+(ix+2))*8121 + 28411 ) % 134456;
	scalefactor_space = ( (scalefactor_space+(iy+2))*8121 + 28411 ) % 134456;
	scalefactor_space = ( (scalefactor_space+(iz+2))*8121 + 28411 ) % 134456;
	scalefactor_space = ( (scalefactor_space+(ix+3*iy+7*iz+2))*8121 + 28411 ) % 134456;
	scalefactor_space = ix+3*iy+7*iz+2;
	scalefactor_space = scalefactor_space & ( (1<<2) - 1 );
	scalefactor_space = 1 << scalefactor_space;

	/*--- ref_faceyz inline ---*/
	faceyz[ia + dims_na      * (
			iu + NU           * (
                        iy + dims_ncell_y * (
                        iz + dims_ncell_z * (
                        ie + dims_ne      * (
                        octant_in_block ))))) ]

	  /*--- Quantities_init_face routine ---*/
	  = Quantities_init_face(ia, ie, iu, scalefactor_space, octant);

        /* *ref_faceyz( sweeper->faceyz, sweeper->dims, NU, */
        /*              Sweeper_noctant_per_block( sweeper ), */
        /*              iy, iz, ie, ia, iu, octant_in_block ) */
        /*      = Quantities_init_faceyz( */
        /*                  quan, ix, iy, iz, ie, ia, iu, octant, sweeper->dims ); */
      }
    }

    /*---Loop over energy groups---*/
    for( ie=0; ie<sweeper->dims.ne; ++ie )
    {
      /*---Calculate spatial loop extents---*/

      int ixbeg = dir_x==DIR_UP ? 0 : sweeper->dims.ncell_x-1;
      int iybeg = dir_y==DIR_UP ? 0 : sweeper->dims.ncell_y-1;
      int izbeg = dir_z==DIR_UP ? 0 : sweeper->dims.ncell_z-1;

      int ixend = dir_x==DIR_DN ? 0 : sweeper->dims.ncell_x-1;
      int iyend = dir_y==DIR_DN ? 0 : sweeper->dims.ncell_y-1;
      int izend = dir_z==DIR_DN ? 0 : sweeper->dims.ncell_z-1;

      /*---Loop over cells, in proper direction---*/

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
      for( ia=0; ia<sweeper->dims.na; ++ia )
      { 
	// reset reduction
	P result = (P)0;
        for( im=0; im < dims.nm; ++im )
        {
	  /*--- const_ref_a_from_m inline ---*/
	  result += (v_a_from_m[ ia     + dims.na * (
	  	       im     +      NM * (
	  	       octant + NOCTANT * (
					   0 ))) ]) * 

	    /*--- const_ref_state inline ---*/
	    (vi_h[im + dims.nm      * (
	    		  iu + NU           * (
	    		  ix + dims.ncell_x * (
                          iy + dims.ncell_y * (
                          ie + dims.ne      * (
                          iz + dims.ncell_z * ( /*---NOTE: This axis MUST be slowest-varying---*/
	    				       0 ))))))]);

          /* result += *const_ref_a_from_m( Pointer_const_h( & quan->a_from_m ), */
          /*                                sweeper->dims, im, ia, octant ) * */
          /*           *const_ref_state(    Pointer_h( vi ), sweeper->dims, NU, */
          /*                                ix, iy, iz, ie, im, iu ); */
        }

	/*--- ref_vslocal inline ---*/
	vs_local[ ia + dims.na * (iu + NU  * (0) ) ] = result;

        /* *ref_vslocal( sweeper->vslocal, sweeper->dims, NU, */
        /*               sweeper->dims.na, ia, iu ) = result; */
      }

      /*--------------------*/
      /*---Perform solve---*/
      /*--------------------*/

      for( ia=0; ia<dims.na; ++ia )
      {
        Quantities_solve( quan, sweeper->vslocal, ia, ia, dims.na,
                          sweeper->facexy, sweeper->facexz, sweeper->faceyz,
                          ix, iy, iz, ie, ix, iy, iz,
                          octant, octant_in_block,
                          Sweeper_noctant_per_block( sweeper ),
                          dims, dims, Bool_true );
      }

      /*--------------------*/
      /*---Transform state vector from angles to moments---*/
      /*--------------------*/

      /*---Perform small dense matrix-vector products and store
           the result in the output state vector.
      ---*/

      for( iu=0; iu<NU; ++iu )
      for( im=0; im<dims.nm; ++im )
      {
        P result = (P)0;
        for( ia=0; ia<dims.na; ++ia )
        {
	  /*--- const_ref_m_from_a ---*/
	  result += v_m_from_a[ im     +      NM * (
         	  	       ia     + dims.na * (
                               octant + NOCTANT * (
                               0 ))) ] *

	    /*--- const_ref_vslocal ---*/
	    vs_local[ ia + dims.na * (
                     iu + NU    * (
                     0 )) ];

          /* result += *const_ref_m_from_a( Pointer_const_h( & quan->m_from_a ), */
          /*                                sweeper->dims, im, ia, octant )* */
	  /* *const_ref_vslocal(  sweeper->vslocal, sweeper->dims, NU, */
	  /*                      sweeper->dims.na, ia, iu ); */
        }

	/*--- ref_state inline ---*/
	vo_h[im + dims.nm      * (
             iu + NU           * (
             ix + dims.ncell_x * (
             iy + dims.ncell_y * (
             ie + dims.ne      * (
             iz + dims.ncell_z * ( /*---NOTE: This axis MUST be slowest-varying---*/
             0 ))))))] += result;

        /* *ref_state( Pointer_h( vo ), sweeper->dims, NU, */
        /*             ix, iy, iz, ie, im, iu ) += result; */
      }

    } /*---ix/iy/iz---*/

    } /*---ie---*/

  } /*---octant---*/

} /*---sweep---*/

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

#endif /*---_sweeper_simple_c_h_---*/

/*---------------------------------------------------------------------------*/
