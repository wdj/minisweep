/*---------------------------------------------------------------------------*/
/*!
 * \file   tester.c
 * \author Wayne Joubert
 * \date   Wed May 22 11:22:14 EDT 2013
 * \brief  Tester for sweep miniapp.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include <stdio.h>

#include "arguments.h"
#include "env.h"
#include "definitions.h"
#include "dimensions.h"
#include "memory.h"
#include "pointer.h"
#include "quantities.h"
#include "array_operations.h"
#include "sweeper.h"

#include "run_tools.h"

#define MAX_LINE_LEN 1024

/*===========================================================================*/
/*---Input a line from standard input---*/

static Bool_t get_line( char* line )
{
  int nchar = 0;
  int c = 0;
    
  while( (c = getchar()) != EOF )
  {
    if( c == '\n' )
    {
      break;
    }

    Assert( nchar + 2 <= MAX_LINE_LEN ? "Input line too long" : 0 );
 
    line[nchar] = c; 
    ++nchar;
  }

  if( c == EOF && nchar == 0 )
  {
    return Bool_false;
  }

  line[nchar] = '\0';
  return Bool_true;
}

/*===========================================================================*/
/*---Tester---*/

static void tester( Env* env )
{
  int ntest = 0;
  int ntest_passed = 0;

  Bool_t result1 = Bool_true;
  Bool_t result2 = Bool_true;

  /*---Loop over pairs of input arg strings to do runs---*/

  while( Bool_true )
  {
    char argstring1[MAX_LINE_LEN];
    char argstring2[MAX_LINE_LEN];

    if( Env_proc_this( env ) == 0 )
    {
      result1 = get_line( argstring1 );
      result2 = get_line( argstring2 );
    }
    Env_bcast_int( env, &result1, Env_proc_this( env ) );
    Env_bcast_int( env, &result2, Env_proc_this( env ) );
    Env_bcast_string( env, argstring1, MAX_LINE_LEN, Env_proc_this( env ) );
    Env_bcast_string( env, argstring2, MAX_LINE_LEN, Env_proc_this( env ) );

    if( ! ( result1 && result2 ) )
    {
      break;
    }

    Bool_t pass = compare_runs( env, argstring1, argstring2 );

    ++ntest;
    if( pass )
    {
      ++ntest_passed;
    }

  }

  if( Env_is_proc_master( env ) )
  {
    printf( "TESTS %i    PASSED %i    FAILED %i\n",
            ntest, ntest_passed, ntest-ntest_passed );
  }
}

/*===========================================================================*/
/*---Main---*/

int main( int argc, char** argv )
{
  /*---Declarations---*/
  Env env;

  /*---Initialize for execution---*/

  Env_initialize( &env, argc, argv );

  /*---Do testing---*/

  tester( &env );

  /*---Finalize execution---*/

  Env_finalize( &env );

} /*---main---*/

/*---------------------------------------------------------------------------*/
