/*---------------------------------------------------------------------------*/
/*!
 * \file   run_tools.h
 * \author Wayne Joubert
 * \date   Wed Jan 28 10:11:10 EST 2015
 * \brief  Declarations for tools to perform runs of sweeper.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _driver__run_h_
#define _driver__run_h_

#include "arguments.h"
#include "env.h"
#include "definitions.h"
#include "memory.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*===========================================================================*/
/*---Struct to hold run result data---*/

typedef struct
{
  P      normsq;
  P      normsqdiff;
  double flops;
  double floprate;
  Timer  time;
} RunData;

/*===========================================================================*/
/*---Perform run---*/

void run_case( Env* env, Arguments* args, RunData* rundata );

/*===========================================================================*/
/*---Perform two runs, compare results---*/

Bool_t compare_runs( Env* env, char* argstring1, char* argstring2 );

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

#endif /*---_driver__run_h_---*/

/*---------------------------------------------------------------------------*/
