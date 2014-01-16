/*---------------------------------------------------------------------------*/
/*!
 * \file   quantities_testing_c.h
 * \author Wayne Joubert
 * \date   Wed Jan 15 16:06:28 EST 2014
 * \brief  Definitions for physical quantities, testing case.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _serial_c__quantities_testing_c_h_
#define _serial_c__quantities_testing_c_h_

#include "dimensions.h"
#include "memory.h"
#include "array_accessors.h"
#include "quantities_testing.h"

/*---------------------------------------------------------------------------*/
/*---Pseudo-constructor for Quantities struct---*/

void Quantities_ctor( Quantities* quan, Dimensions dims )
{
  /*---Declarations---*/

  int im = 0;
  int ia = 0;
  int i  = 0;

  /*---Allocate arrays---*/

  quan->a_from_m = pmalloc( dims.nm * dims.na );
  quan->m_from_a = pmalloc( dims.nm * dims.na );

  /*-----------------------------*/
  /*---Set entries of a_from_m---*/
  /*-----------------------------*/

  /*---These two matrices are set in a special way so as to map a vector
       whose values satisfy a linear relationship v_i = a * i + b
       are mapped to another vector with same property.  This is to facilitate
       being able to have an analytical solution for the sweep.
  ---*/

  /*---First set to zero---*/

  for( im=0; im<dims.nm; ++im )
  for( ia=0; ia<dims.na; ++ia )
  {
    *ref_a_from_m( quan->a_from_m, dims, im, ia ) = P_zero();
  }

  /*---Map a linear vector of values to a similar linear vector of values---*/

  /*---The input state vector in the moment dimension is contrived to
       satisfy vi[im] = 1 + im, an affine function, possibly times a constant.
       The following matrix is artifically contrived to send this to
       a result satisfying, in angle, vl[ia] = 1 + ia, possibly times a
       constant.
  ---*/

  for( i=0; i<dims.na; ++i )
  {
    const int quot = ( i + 1 ) / dims.nm;
    const int rem  = ( i + 1 ) % dims.nm;

    *ref_a_from_m( quan->a_from_m, dims, dims.nm-1, i ) += quot;
    if( rem != 0 )
    {
      *ref_a_from_m( quan->a_from_m, dims, 0,   i ) += -P_one();
      *ref_a_from_m( quan->a_from_m, dims, rem, i ) +=  P_one();
    }
  }

  /*---Fill matrix with entries that leave linears unaffected---*/

  /*---This is to create a more dense, nontrivial matrix, with additions
       to the rows that are guaranteed to send affine functions to zero.
  ---*/

  for( im=0; im<dims.nm-2; ++im )
  for( ia=0; ia<dims.na;   ++ia )
  {
    const int randvalue = 21 + ( im + dims.nm * ia ) % 17;

    *ref_a_from_m( quan->a_from_m, dims, im+0, ia ) +=  -P_one() * randvalue;
    *ref_a_from_m( quan->a_from_m, dims, im+1, ia ) += 2*P_one() * randvalue;
    *ref_a_from_m( quan->a_from_m, dims, im+2, ia ) +=  -P_one() * randvalue;
  }

  /*-----------------------------*/
  /*---Set entries of m_from_a---*/
  /*-----------------------------*/

  /*---First set to zero---*/

  for( im=0; im<dims.nm; ++im )
  for( ia=0; ia<dims.na; ++ia )
  {
    *ref_m_from_a( quan->m_from_a, dims, im, ia ) = P_zero();
  }

  /*---Map a linear vector of values to a similar linear vector of values---*/

  /*---As previously, send functions vl[ia] = 1 + ia to functions
       vo[im] = 1 + im, subject to possible constant scalings
       and also to a power-of-two angle scalefactor adjustment
       designed to make the test more rigorous.
  ---*/

  for( i=0; i<dims.nm; ++i )
  {
    const int quot = ( i + 1 ) / dims.na;
    const int rem  = ( i + 1 ) % dims.na;

    *ref_m_from_a( quan->m_from_a, dims, i, dims.na-1 ) += quot;
    if( rem != 0 )
    {
      *ref_m_from_a( quan->m_from_a, dims, i, 0   ) += -P_one();
      *ref_m_from_a( quan->m_from_a, dims, i, rem ) +=  P_one();
    }
  }

  /*---Fill matrix with entries that leave linears unaffected---*/

  /*---As before, create more complex matrix by adding to rows
       entries that do not affect the scaled-affine input values expected.
  ---*/

  for( im=0; im<dims.nm;   ++im )
  for( ia=0; ia<dims.na-2; ++ia )
  {
    const int randvalue = 37 + ( im + dims.nm * ia ) % 19;

    *ref_m_from_a( quan->m_from_a, dims, im, ia+0 ) +=  -P_one() * randvalue;
    *ref_m_from_a( quan->m_from_a, dims, im, ia+1 ) += 2*P_one() * randvalue;
    *ref_m_from_a( quan->m_from_a, dims, im, ia+2 ) +=  -P_one() * randvalue;
  }

  /*---Scale matrix to compensate for 8 octants and also angle scale factor---*/

  for( im=0; im<dims.nm; ++im )
  for( ia=0; ia<dims.na; ++ia )
  {
    *ref_m_from_a( quan->m_from_a, dims, im, ia ) /= NOCTANT;
    *ref_m_from_a( quan->m_from_a, dims, im, ia ) /=
                                         Quantities_scalefactor_angle__( ia );
  }

} /*---Quantities_ctor---*/

/*---------------------------------------------------------------------------*/
/*---Pseudo-destructor for Quantities struct---*/

void Quantities_dtor( Quantities* quan )
{
  /*---Deallocate arrays---*/

  pfree( quan->a_from_m );
  pfree( quan->m_from_a );

} /*---Quantities_dtor---*/

/*---------------------------------------------------------------------------*/

#endif /*---_serial_c__quantities_testing_c_h_---*/

/*---------------------------------------------------------------------------*/
