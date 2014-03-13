/*---------------------------------------------------------------------------*/
/*!
 * \file   env.h
 * \author Wayne Joubert
 * \date   Wed Jan 15 16:06:28 EST 2014
 * \brief  Environment settings specific to this programming API.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _serial_c__env_h_
#define _serial_c__env_h_

#include <stdio.h>
#include <stdlib.h>
#include "assert.h"

#include <sys/time.h>

#ifdef USE_MPI
#include "mpi.h"
#endif

#ifdef USE_OPENMP
#include "omp.h"
#endif

/*===========================================================================*/
/*---Assertions---*/

#ifndef Insist
#define Insist( condition ) \
  (void)((condition) || (insist_ (#condition, __FILE__, __LINE__),0))
#endif

static void insist_( const char *condition_string, const char *file, int line )
{
  fprintf( stderr, "Insist error: \"%s\". At file %s, line %i.\n",
                   condition_string, file, line );
  exit( EXIT_FAILURE );
}

#ifndef NDEBUG
#define Static_Assert( condition ) { int a[ ( condition ) ? 1 : -1 ]; }
#else
#define Static_Assert( condition )
#endif

/*===========================================================================*/
/*---Basic types---*/

/*---Floating point type for sweep---*/

typedef double P;
enum{ P_IS_DOUBLE = 1 };

static inline P P_zero() { return (P)0.; }
static inline P P_one()  { return (P)1.; }

/*===========================================================================*/
/*---Boolean type---*/

typedef int Bool_t;

enum{ Bool_true = 1, Bool_false = 0 };

/*===========================================================================*/
/*---Struct with environment information---*/

typedef struct
{
  int nproc_x;    /*---Number of procs along x axis---*/
  int nproc_y;    /*---Number of procs along y axis---*/
  int tag;        /*---Next free message tag---*/
} Env;

/*===========================================================================*/
/*---Initialize for execution---*/

static void Env_initialize( Env *env, int argc, char** argv )
{
#ifdef USE_MPI
  MPI_Init( &argc, &argv );
#endif

  env->tag = 0;
}

/*===========================================================================*/
/*---Finalize execution---*/

static void Env_finalize( Env* env )
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

static int Env_nproc( const Env* env )
{
  int result = 1;
#ifdef USE_MPI
  MPI_Comm_size( Env_default_comm(), &result );
#endif
  return result;
}

/*---------------------------------------------------------------------------*/

static int Env_nproc_x( const Env* env )
{
  return env->nproc_x;
}

/*---------------------------------------------------------------------------*/

static int Env_nproc_y( const Env* env )
{
  return env->nproc_y;
}

/*===========================================================================*/
/*---Proc number---*/

static int Env_proc( const Env* env, int proc_x, int proc_y )
{
  assert( proc_x >= 0 && proc_x < Env_nproc_x( env ) );
  assert( proc_y >= 0 && proc_y < Env_nproc_y( env ) );
  int result = proc_x + Env_nproc_x( env ) * proc_y;
  return result;
}

/*---------------------------------------------------------------------------*/

static int Env_proc_x( const Env* env, int proc )
{
  assert( proc >= 0 && proc < Env_nproc( env ) );
  int result = proc % Env_nproc_x( env );
  return result;
}

/*---------------------------------------------------------------------------*/

static int Env_proc_y( const Env* env, int proc )
{
  assert( proc >= 0 && proc < Env_nproc( env ) );
  int result = proc / Env_nproc_x( env );
  return result;
}

/*===========================================================================*/
/*---Proc number this proc---*/

static int Env_proc_this( const Env* env )
{
  int result = 0;
#ifdef USE_MPI
  MPI_Comm_rank( Env_default_comm(), &result );
#endif
  return result;
}

/*---------------------------------------------------------------------------*/

static int Env_proc_x_this( const Env* env )
{
  return Env_proc_x( env, Env_proc_this( env ) );
}

/*---------------------------------------------------------------------------*/

static int Env_proc_y_this( const Env* env )
{
  return Env_proc_y( env, Env_proc_this( env ) );
}

/*===========================================================================*/
/*---MPI functions: global---*/

static void Env_barrier()
{
#ifdef USE_MPI
  MPI_Barrier( Env_default_comm() );
#endif
}

/*---------------------------------------------------------------------------*/

static double Env_sum_d( double value )
{
  double result = 0.;
#ifdef USE_MPI
  MPI_Allreduce( &value, &result, 1, MPI_DOUBLE, MPI_SUM, Env_default_comm() );
#else
  result = value;
#endif
  return result;
}

/*---------------------------------------------------------------------------*/

static P Env_sum_P( P value )
{
  Static_Assert( P_IS_DOUBLE );
  return Env_sum_d( value );
}

/*===========================================================================*/
/*---MPI functions: point to point---*/

#ifdef USE_MPI
typedef MPI_Request Request_t;
#else
typedef int Request_t;
#endif

static void Env_send_i( const int* data, size_t n, int proc, int tag )
{
  assert( data != NULL );

#ifdef USE_MPI
  MPI_Send( (void*)data, n, MPI_INT, proc, tag, Env_default_comm() );
#endif
}

/*---------------------------------------------------------------------------*/

static void Env_recv_i( int* data, size_t n, int proc, int tag )
{
  assert( data != NULL );

#ifdef USE_MPI
  MPI_Status status;
  MPI_Recv( (void*)data, n, MPI_INT, proc, tag, Env_default_comm(), &status );
#endif
}

/*---------------------------------------------------------------------------*/

static void Env_send_P( const P* data, size_t n, int proc, int tag )
{
  Static_Assert( P_IS_DOUBLE );
  assert( data != NULL );
  assert( n >= 0 );

#ifdef USE_MPI
  MPI_Send( (void*)data, n, MPI_DOUBLE, proc, tag, Env_default_comm() );
#endif
}

/*---------------------------------------------------------------------------*/

static void Env_recv_P( P* data, size_t n, int proc, int tag )
{
  Static_Assert( P_IS_DOUBLE );
  assert( data != NULL );
  assert( n >= 0 );

#ifdef USE_MPI
  MPI_Status status;
  MPI_Recv( (void*)data, n, MPI_DOUBLE, proc, tag,
                                                 Env_default_comm(), &status );
#endif
}

/*---------------------------------------------------------------------------*/

static void Env_asend_P( const P* data, size_t n, int proc, int tag,
                                                          Request_t* request )
{
  Static_Assert( P_IS_DOUBLE );
  assert( data != NULL );
  assert( n >= 0 );
  assert( request != NULL );

#ifdef USE_MPI
  MPI_Isend( (void*)data, n, MPI_DOUBLE, proc, tag, Env_default_comm(),
                                                                     request );
#endif
}

/*---------------------------------------------------------------------------*/

static void Env_arecv_P( const P* data, size_t n, int proc, int tag,
                                                          Request_t* request )
{
  Static_Assert( P_IS_DOUBLE );
  assert( data != NULL );
  assert( n >= 0 );
  assert( request != NULL );

#ifdef USE_MPI
  MPI_Irecv( (void*)data, n, MPI_DOUBLE, proc, tag, Env_default_comm(),
                                                                     request );
#endif
}

/*---------------------------------------------------------------------------*/

static void Env_wait( Request_t* request )
{
#ifdef USE_MPI
  MPI_Status status;
  MPI_Waitall( 1, request, &status );
#endif
}

/*===========================================================================*/
/*---Get openmp current thread number---*/

static inline int Env_thread_this( const Env* env )
{
  int result = 0;
#ifdef USE_OPENMP
  result = omp_get_thread_num();  
#endif
  return result;
}

/*===========================================================================*/
/*---Indicate whether to do output---*/

static Bool_t Env_do_output( const Env* env )
{
  return ( Env_proc_this( env ) == 0 );
}

/*===========================================================================*/
/*---Timer type---*/

typedef double Timer_t;

/*===========================================================================*/
/*---Timer utilities---*/

static Timer_t Env_get_time()
{
    struct timeval tv;
    int i = gettimeofday( &tv, NULL );
    Timer_t result = ( (Timer_t) tv.tv_sec +
                       (Timer_t) tv.tv_usec * 1.e-6 );
    return result;
}

/*---------------------------------------------------------------------------*/

static Timer_t Env_get_synced_time()
{
  Env_barrier();
  return Env_get_time();
}

/*===========================================================================*/

#endif /*---_serial_c__env_h_---*/

/*---------------------------------------------------------------------------*/
