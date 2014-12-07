/*---------------------------------------------------------------------------*/
/*!
 * \file   env_mpi.h
 * \author Wayne Joubert
 * \date   Wed Jan 15 16:06:28 EST 2014
 * \brief  Environment settings for MPI.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _env_mpi_h_
#define _env_mpi_h_

#include "types.h"
#include "arguments.h"
#include "env_assert.h"
#include "env_declarations.h"

#ifdef USE_MPI
#include "mpi.h"
#endif

#ifdef __cplusplus
extern "C"
{
#endif

/*===========================================================================*/
/*---Initialize mpi---*/

static void Env_mpi_initialize__( Env *env, int argc, char** argv )
{
#ifdef USE_MPI
  const int mpi_code = MPI_Init( &argc, &argv );
  Assert( mpi_code == MPI_SUCCESS );
  env->active_comm__ = MPI_COMM_NULL;
  env->tag__ = 0;
#endif
}

/*===========================================================================*/
/*---Finalize mpi---*/

static void Env_mpi_finalize__( Env* env )
{
#ifdef USE_MPI
  int mpi_code = 0;
  if( env->active_comm__ != MPI_COMM_NULL )
  {
    mpi_code = MPI_Comm_free( &env->active_comm__ );
    Assert( mpi_code == MPI_SUCCESS );
  }
  mpi_code = MPI_Finalize();
  Assert( mpi_code == MPI_SUCCESS );
#endif
}

/*===========================================================================*/
/*---Default communicator---*/

static Comm_t Env_mpi_active_comm__( const Env* env )
{
#ifdef USE_MPI
  return env->active_comm__;
#else
  return 0;
#endif
}

/*===========================================================================*/
/*---Number of procs---*/

static int Env_nproc_x( const Env* env )
{
  int result = 1;
#ifdef USE_MPI
  result = env->nproc_x__;
#endif
  return result;
}

/*---------------------------------------------------------------------------*/

static int Env_nproc_y( const Env* env )
{
  int result = 1;
#ifdef USE_MPI
  result = env->nproc_y__;
#endif
  return result;
}

/*---------------------------------------------------------------------------*/

static int Env_nproc( const Env* env )
{
  return Env_nproc_x( env ) * Env_nproc_y( env );
}

/*===========================================================================*/
/*---Is this proc within the communicator of procs in use---*/

static Bool_t Env_is_proc_active( const Env* env )
{
  Bool_t result = Bool_true;
#ifdef USE_MPI
  result = env->is_proc_active__;
#endif
  return result;
}

/*===========================================================================*/
/*---Set values from args---*/

static void Env_mpi_set_values__( Env *env, Arguments* args )
{
#ifdef USE_MPI
  int mpi_code = 0;

  env->nproc_x__ = Arguments_consume_int_or_default( args, "--nproc_x", 1 );
  env->nproc_y__ = Arguments_consume_int_or_default( args, "--nproc_y", 1);
  Insist( Env_nproc_x( env ) > 0 ? "Invalid nproc_x supplied." : 0 );
  Insist( Env_nproc_y( env ) > 0 ? "Invalid nproc_y supplied." : 0 );

  int nproc_requested = env->nproc_x__ * env->nproc_y__;
  int nproc_world = 0;
  mpi_code = MPI_Comm_size( MPI_COMM_WORLD, &nproc_world );
  Assert( mpi_code == MPI_SUCCESS );
  Insist( nproc_requested <= nproc_world ?
                                      "Not enough processors available." : 0 );

  int rank = 0;
  mpi_code = MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  Assert( mpi_code == MPI_SUCCESS );

  if( env->active_comm__ != MPI_COMM_NULL )
  {
    mpi_code = MPI_Comm_free( &env->active_comm__ );
    Assert( mpi_code == MPI_SUCCESS );
  }

  env->is_proc_active__ = rank < nproc_requested;
  mpi_code = MPI_Comm_split( MPI_COMM_WORLD, env->is_proc_active__,
                                                   rank, &env->active_comm__ );
  Assert( mpi_code == MPI_SUCCESS );

#endif
}

/*===========================================================================*/
/*---Tag manipulation---*/

static int Env_tag( const Env* env )
{
  int result = 0;
#ifdef USE_MPI
  result = env->tag__;
#endif
  return result;
}

/*---------------------------------------------------------------------------*/

static void Env_increment_tag( Env* env, int value )
{
#ifdef USE_MPI
  env->tag__ += value;
#endif
}

/*===========================================================================*/
/*---Proc number---*/

static int Env_proc( const Env* env, int proc_x, int proc_y )
{
  Assert( proc_x >= 0 && proc_x < Env_nproc_x( env ) );
  Assert( proc_y >= 0 && proc_y < Env_nproc_y( env ) );
  int result = proc_x + Env_nproc_x( env ) * proc_y;
  return result;
}

/*---------------------------------------------------------------------------*/

static int Env_proc_x( const Env* env, int proc )
{
  Assert( proc >= 0 && proc < Env_nproc( env ) );
  int result = proc % Env_nproc_x( env );
  return result;
}

/*---------------------------------------------------------------------------*/

static int Env_proc_y( const Env* env, int proc )
{
  Assert( proc >= 0 && proc < Env_nproc( env ) );
  int result = proc / Env_nproc_x( env );
  return result;
}

/*===========================================================================*/
/*---Proc number this proc---*/

static int Env_proc_this( const Env* env )
{
  int result = 0;
#ifdef USE_MPI
  const int mpi_code = MPI_Comm_rank( Env_mpi_active_comm__( env ), &result );
  Assert( mpi_code == MPI_SUCCESS );
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

static void Env_mpi_barrier( Env* env )
{
#ifdef USE_MPI
  const int mpi_code = MPI_Barrier( Env_mpi_active_comm__( env ) );
  Assert( mpi_code == MPI_SUCCESS );
#endif
}

/*---------------------------------------------------------------------------*/

static double Env_sum_d( double value, Env* env )
{
  double result = 0.;
#ifdef USE_MPI
  const int mpi_code = MPI_Allreduce( &value, &result, 1, MPI_DOUBLE, MPI_SUM,
                                                Env_mpi_active_comm__( env ) );
  Assert( mpi_code == MPI_SUCCESS );
#else
  result = value;
#endif
  return result;
}

/*---------------------------------------------------------------------------*/

static P Env_sum_P( P value, Env* env )
{
  Static_Assert( P_IS_DOUBLE );
  return Env_sum_d( value, env );
}

/*===========================================================================*/
/*---MPI functions: point to point---*/

static void Env_send_i( const int* data, size_t n, int proc, int tag, Env* env )
{
  Assert( data != NULL );

#ifdef USE_MPI
  const int mpi_code = MPI_Send( (void*)data, n, MPI_INT, proc, tag,
                                                Env_mpi_active_comm__( env ) );
  Assert( mpi_code == MPI_SUCCESS );
#endif
}

/*---------------------------------------------------------------------------*/

static void Env_recv_i( int* data, size_t n, int proc, int tag, Env* env )
{
  Assert( data != NULL );

#ifdef USE_MPI
  MPI_Status status;
  const int mpi_code = MPI_Recv( (void*)data, n, MPI_INT, proc, tag,
                                       Env_mpi_active_comm__( env ), &status );
  Assert( mpi_code == MPI_SUCCESS );
#endif
}

/*---------------------------------------------------------------------------*/

static void Env_send_P( const P* data, size_t n, int proc, int tag, Env* env )
{
  Static_Assert( P_IS_DOUBLE );
  Assert( data != NULL );
  Assert( n+1 >= 1 );

#ifdef USE_MPI
  const int mpi_code = MPI_Send( (void*)data, n, MPI_DOUBLE, proc, tag,
                                                Env_mpi_active_comm__( env ) );
  Assert( mpi_code == MPI_SUCCESS );
#endif
}

/*---------------------------------------------------------------------------*/

static void Env_recv_P( P* data, size_t n, int proc, int tag, Env* env )
{
  Static_Assert( P_IS_DOUBLE );
  Assert( data != NULL );
  Assert( n+1 >= 1 );

#ifdef USE_MPI
  MPI_Status status;
  const int mpi_code = MPI_Recv( (void*)data, n, MPI_DOUBLE, proc, tag,
                                       Env_mpi_active_comm__( env ), &status );
  Assert( mpi_code == MPI_SUCCESS );
#endif
}

/*---------------------------------------------------------------------------*/

static void Env_asend_P( const P* data, size_t n, int proc, int tag,
                                                 Request_t* request, Env* env )
{
  Static_Assert( P_IS_DOUBLE );
  Assert( data != NULL );
  Assert( n+1 >= 1 );
  Assert( request != NULL );

#ifdef USE_MPI
  const int mpi_code = MPI_Isend( (void*)data, n, MPI_DOUBLE, proc, tag,
                                       Env_mpi_active_comm__( env ), request );
  Assert( mpi_code == MPI_SUCCESS );
#endif
}

/*---------------------------------------------------------------------------*/

static void Env_arecv_P( const P* data, size_t n, int proc, int tag,
                                                 Request_t* request, Env* env )
{
  Static_Assert( P_IS_DOUBLE );
  Assert( data != NULL );
  Assert( n+1 >= 1 );
  Assert( request != NULL );

#ifdef USE_MPI
  const int mpi_code = MPI_Irecv( (void*)data, n, MPI_DOUBLE, proc, tag,
                                       Env_mpi_active_comm__( env ), request );
  Assert( mpi_code == MPI_SUCCESS );
#endif
}

/*---------------------------------------------------------------------------*/

static void Env_wait( Request_t* request )
{
#ifdef USE_MPI
  MPI_Status status;
  const int mpi_code = MPI_Waitall( 1, request, &status );
  Assert( mpi_code == MPI_SUCCESS );
#endif
}

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

#endif /*---_env_mpi_h_---*/

/*---------------------------------------------------------------------------*/
