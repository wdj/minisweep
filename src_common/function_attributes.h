/*---------------------------------------------------------------------------*/
/*!
 * \file   function_attributes.h
 * \author Wayne Joubert
 * \date   Fri Apr 25 10:41:43 EDT 2014
 * \brief  Function attribute definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _function_attributes_h_
#define _function_attributes_h_

/*===========================================================================*/

#ifdef __CUDACC__

#define TARGET_G  __global__
#define TARGET_HD __host__ __device__

#else

#define TARGET_G
#define TARGET_HD

#endif

/*===========================================================================*/

#endif /*---_function_attributes_h_---*/

/*---------------------------------------------------------------------------*/
