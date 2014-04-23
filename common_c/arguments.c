/*---------------------------------------------------------------------------*/
/*!
 * \file   arguments.h
 * \author Wayne Joubert
 * \date   Wed Mar 12 14:19:26 EDT 2014
 * \brief  Definitions for Arguments struct.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>

#include "arguments.h"
#include "types.h"
#include "env_assert.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*===========================================================================*/
/* Pseudo-constructor for Arguments struct---*/

void Arguments_ctor( Arguments* args,
                     int        argc,
                     char**     argv )
{
  assert( argc > 0 );
  assert( argv != NULL );
  int i = 0;

  args->argc = argc;
  args->argv_unconsumed = (char**) malloc( argc * sizeof( char*) );

  for( i=0; i<argc; ++i )
  {
    assert( argv[i] != NULL && "Null command line argument encountered" );
    args->argv_unconsumed[i] = argv[i];
  }
} /*---Arguments_ctor---*/

/*===========================================================================*/
/* Pseudo-destructor for Arguments struct---*/

void Arguments_dtor( Arguments* args )
{
  free( (void*) args->argv_unconsumed );

} /*---Arguments_dtor---*/

/*===========================================================================*/
/* Determine whether an argument with a given name exists---*/

Bool_t Arguments_exists( Arguments*  args,
                         const char* arg_name )
{
  Bool_t result = Bool_false;
  int    i      = 0;

  for( i=0; i<args->argc; ++i )
  {
    if( args->argv_unconsumed[i] == NULL )
    {
      continue;
    }
    result = result || strcmp( args->argv_unconsumed[i], arg_name ) == 0;
  }

  return result;
}

/*===========================================================================*/
/* Process an argument of type int, remove from list---*/

int Arguments_consume_int( Arguments*  args,
                           const char* arg_name )
{
  int    result = 0;
  Bool_t found  = Bool_false;
  int    i      = 0;

  for( i=0; i<args->argc; ++i )
  {
    if( args->argv_unconsumed[i] == NULL )
    {
      continue;
    }
    if( strcmp( args->argv_unconsumed[i], arg_name ) == 0 )
    {
      found = Bool_true;
      args->argv_unconsumed[i] = NULL;
      ++i;
      Insist( i<args->argc );
      result = atoi( args->argv_unconsumed[i] );
      args->argv_unconsumed[i] = NULL;
    }
  }

  Insist( found && "Invalid use of argument." );
  return result;
}

/*===========================================================================*/
/* Consume an argument of type int, if not present then set to a default---*/

int Arguments_consume_int_or_default( Arguments*  args,
                                      const char* arg_name,
                                      int         default_value )
{
  return Arguments_exists( args, arg_name ) ?
                      Arguments_consume_int( args, arg_name ) : default_value;
}

/*===========================================================================*/
/* Determine whether all arguments have been consumed---*/

Bool_t Arguments_are_all_consumed( Arguments* args )
{
  Bool_t result = Bool_true;
  int    i      = 0;

  for( i=1; i<args->argc; ++i ) /*---Note: skip the zeroth element---*/
  {
    result = result && args->argv_unconsumed[i] == NULL;
  }

  return result;
}

/*===========================================================================*/
  
#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*===========================================================================*/
