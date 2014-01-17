/*---------------------------------------------------------------------------*/
/*!
 * \file   environment.h
 * \author Wayne Joubert
 * \date   Wed Jan 15 16:06:28 EST 2014
 * \brief  Environment settings specific to this programming API.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _serial_c__environment_h_
#define _serial_c__environment_h_

#include <stdio.h>
#include <stdlib.h>
#include "assert.h"

#include <sys/time.h>

#ifdef USE_MPI
#include "mpi.h"
#endif

#ifndef Insist
#define Insist( condition, message ) \
  (void)((condition) || (insist_ (#condition, __FILE__, __LINE__, message),0))
#endif

/*===========================================================================*/

#define PROGRAMMING_API serial

/*===========================================================================*/
/*---Insist: an always-on assert---*/

static void insist_( const char *condition_string, const char *file, int line,
                     const char *message )
{
  fprintf( stderr, "Insist error: \"%s\". At file %s, line %i. %s\n",
                   condition_string, file, line, message );
  exit( EXIT_FAILURE );
}

/*===========================================================================*/
/*---Initialize for execution---*/
  
static void initialize( int argc, char** argv )
{ 
#ifdef USE_MPI
  MPI_Init( &argc, &argv );
#endif
}
  
/*===========================================================================*/
/*---Finalize execution---*/

static void finalize() 
{ 
#ifdef USE_MPI
  MPI_Finalize(); 
#endif
} 
/*===========================================================================*/
/*---Indicate whether to do output---*/

static int do_output()
{
#ifdef USE_MPI
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  return ( rank == 0 );
#else
  return 1;
#endif
}

/*===========================================================================*/
/*---Timer utility---*/
  
static double get_time() 
{
    struct timeval tp;
    struct timezone tzp;
    int i = gettimeofday( &tp, &tzp );
    double result = ( (double) tp.tv_sec +
                      (double) tp.tv_usec * 1.e-6 );
    return result;
} 

/*===========================================================================*/

#endif /*---_serial_c__environment_h_---*/

/*---------------------------------------------------------------------------*/
