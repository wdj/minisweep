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

/*===========================================================================*/
/*---Assertions---*/

#ifndef Insist
#define Insist( condition, message ) \
  (void)((condition) || (insist_ (#condition, __FILE__, __LINE__, message),0))
#endif

static void insist_( const char *condition_string, const char *file, int line,
                     const char *message )
{
  fprintf( stderr, "Insist error: \"%s\". At file %s, line %i. %s\n",
                   condition_string, file, line, message );
  exit( EXIT_FAILURE );
}

/*===========================================================================*/
/*---Boolean type---*/

typedef int Bool_t;

/*===========================================================================*/
/*---Timer utility---*/
  
static double Env_get_time() 
{
    struct timeval tv;
    int i = gettimeofday( &tv, NULL );
    double result = ( (double) tv.tv_sec +
                      (double) tv.tv_usec * 1.e-6 );
    return result;
} 

/*===========================================================================*/
/*---Initialize for execution---*/
  
static void Env_initialize( int argc, char** argv )
{ 
#ifdef USE_MPI
  MPI_Init( &argc, &argv );
#endif
}
  
/*===========================================================================*/
/*---Finalize execution---*/

static void Env_finalize() 
{ 
#ifdef USE_MPI
  MPI_Finalize(); 
#endif
}

/*===========================================================================*/
/*---Default communicator---*/

#ifdef USE_MPI
typedef MPI_Comm Comm_t;
#else
typedef int Comm_t;
#endif

static Comm_t Env_default_comm()
{
#ifdef USE_MPI
  return MPI_COMM_WORLD;
#else
  return 0;
#endif
}

/*===========================================================================*/
/*---Number of procs---*/

static int Env_nproc()
{
  int result = 0;
#ifdef USE_MPI
  MPI_Comm_size( Env_default_comm(), &result );
#endif
  return result;
}

/*===========================================================================*/
/*---Proc number---*/

static int Env_proc()
{
  int result = 0;
#ifdef USE_MPI
  MPI_Comm_rank( Env_default_comm(), &result );
#endif
  return result;
}
 
/*===========================================================================*/
/*---Indicate whether to do output---*/

static Bool_t Env_do_output()
{
  return ( Env_proc() == 0 );
}

/*===========================================================================*/

#endif /*---_serial_c__environment_h_---*/

/*---------------------------------------------------------------------------*/
