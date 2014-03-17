/*---------------------------------------------------------------------------*/
/*!
 * \file   env_mpi.h
 * \author Wayne Joubert
 * \date   Wed Jan 15 16:06:28 EST 2014
 * \brief  Environment settings for MPI.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _mpi_c__env_mpi_h_
#define _mpi_c__env_mpi_h_

#include "env_assert.h"

#ifdef USE_MPI
#include "mpi.h"
#endif

/*===========================================================================*/
/*---Initialize mpi---*/

static void Env_initialize_mpi__( Env *env, int argc, char** argv )
{
#ifdef USE_MPI
  MPI_Init( &argc, &argv );
#endif
}

/*===========================================================================*/
/*---Finalize mpi---*/

static void Env_finalize_mpi__( Env* env )
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

#endif /*---_mpi_c__env_mpi_h_---*/

/*---------------------------------------------------------------------------*/
