/*---------------------------------------------------------------------------*/
/*! 
 * \file   gauss_seidel.c
 * \author Wayne Joubert
 * \date   Mon Jul 27 11:38:40 EDT 2015
 * \brief  Code for Gauss-Seidel solve with block wavefront method.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

/*=============================================================================

This code performs Gauss-Seidel iteration to solve a specific linear system
of equations.  It is known that the Gauss-Seidel calculation has recursive
dependencies which are amenable to a wavefront method for parallel solves.

The problem solved here is a 5-point finite difference stencil problem on a
regular 2D grid.  Standard lexicographical ordering is applied to the
grid unknowns.

=============================================================================*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <sys/time.h>

typedef double Float_t;

/*===========================================================================*/
/*---Timer function---*/

double get_time()
{
    struct timeval tv;
    int i = gettimeofday (&tv, NULL);
    return ( (double) tv.tv_sec +
             (double) tv.tv_usec * 1.e-6 );
}

/*===========================================================================*/
/*---Utility functions---*/

int min(int i, int j)
{
    return i < j ? i : j;
}

/*---------------------------------------------------------------------------*/

int max(int i, int j)
{
    return i > j ? i : j;
}

/*---------------------------------------------------------------------------*/

int ceil_(int i, int j)
{
    return (i + j - 1 ) / j;
}

/*===========================================================================*/
/*---Accessor functions to access elements of array with multidim index---*/

Float_t* ref(Float_t* __restrict__ v, int ix, int iy, int ncellx, int ncelly)
{
    assert(ix >= 0 && ix < ncellx);
    assert(iy >= 0 && iy < ncelly);

    return &( v[ix+ncellx*iy] );
}

/*---------------------------------------------------------------------------*/

const Float_t* const_ref(const Float_t* __restrict__ v, int ix, int iy,
                         int ncellx, int ncelly)
{
    assert(ix >= 0 && ix < ncellx);
    assert(iy >= 0 && iy < ncelly);

    return &( v[ix+ncellx*iy] );
}

/*===========================================================================*/
/*---5-point stencil operation applied at a gridcell---*/

void process_cell(int ix, int iy, int ncellx, int ncelly,
                        Float_t* __restrict__ vo,
                  const Float_t* __restrict__ vi,
                  const Float_t* __restrict__ an,
                  const Float_t* __restrict__ as,
                  const Float_t* __restrict__ ae,
                  const Float_t* __restrict__ aw)
{
    /*---South, east connections use "old" values, north, west use "new"---*/
    *ref(      vo, ix,   iy,   ncellx, ncelly) =

    *const_ref(vo, ix-1, iy,   ncellx, ncelly) *
    *const_ref(aw, ix,   iy,   ncellx, ncelly) +

    *const_ref(vo, ix,   iy-1, ncellx, ncelly) *
    *const_ref(an, ix,   iy,   ncellx, ncelly) +

    *const_ref(vi, ix+1, iy,   ncellx, ncelly) *
    *const_ref(ae, ix,   iy,   ncellx, ncelly) +

    *const_ref(vi, ix,   iy+1, ncellx, ncelly) *
    *const_ref(as, ix,   iy,   ncellx, ncelly);
}

/*===========================================================================*/
/*---Perform Gauss-Seidel sweep with standard order of operations---*/

void solve_standard(int ncellx, int ncelly,
                          Float_t* __restrict__ vo,
                    const Float_t* __restrict__ vi,
                    const Float_t* __restrict__ an,
                    const Float_t* __restrict__ as,
                    const Float_t* __restrict__ ae,
                    const Float_t* __restrict__ aw)
{
    /*---Loop over interior gridcells---*/
    int iy = 0;
    for (iy=1; iy<ncelly-1; ++iy)
    {
        int ix = 0;
        for (ix=1; ix<ncellx-1; ++ix)
        {
            process_cell(ix, iy, ncellx, ncelly, vo, vi, an, as, ae, aw);
        }
    }
}

/*===========================================================================*/
/*---Perform Gauss-Seidel sweep with wavefront order of operations---*/
/*---Note all required recursion dependencies are satisfied---*/

void solve_wavefronts(int ncellx, int ncelly,
                            Float_t* __restrict__ vo,
                      const Float_t* __restrict__ vi,
                      const Float_t* __restrict__ an,
                      const Float_t* __restrict__ as,
                      const Float_t* __restrict__ ae,
                      const Float_t* __restrict__ aw)
{
    const int nwavefronts = ncellx + ncelly - 2;
    int wavefront = 0;

    /*---Loop over wavefronts---*/
    for (wavefront=0; wavefront<nwavefronts; ++wavefront)
    {
        const int ixmin = max(1, wavefront-ncelly+2);
        const int ixmax = min(ncellx-1, wavefront);
        int ix = 0;
        /*---Loop over gridcells in wavefront---*/
        for (ix=ixmin; ix<ixmax; ++ix)
        {
            const int iy = wavefront - ix;
            process_cell(ix, iy, ncellx, ncelly, vo, vi, an, as, ae, aw);
        }
    }
}

/*===========================================================================*/
/*---Perform Gauss-Seidel sweep with block-wavefront order of operations---*/
/*---Note all required recursion dependencies are satisfied---*/

void solve_block_wavefronts(int ncellx, int ncelly,
                                  Float_t* __restrict__ vo,
                            const Float_t* __restrict__ vi,
                            const Float_t* __restrict__ an,
                            const Float_t* __restrict__ as,
                            const Float_t* __restrict__ ae,
                            const Float_t* __restrict__ aw)
{
    /*---Set block sizes---*/
    const int nbsizex = 20;
    const int nbsizey = 20;

    const int nbx = ceil_(ncellx, nbsizex);
    const int nby = ceil_(ncelly, nbsizey);

    const int nbwavefronts = nbx + nby - 1;
    int bwavefront = 0;

    /*---Loop over block wavefronts---*/
    for (bwavefront=0; bwavefront<nbwavefronts; ++bwavefront)
    {
        const int ibxmin = max(0, bwavefront-nby);
        const int ibxmax = min(nbx, bwavefront+1);
        int ibx = 0;
        /*---Loop over blocks in wavefront---*/
        for (ibx=ibxmin; ibx<ibxmax; ++ibx)
        {
            const int iby = bwavefront - ibx;

            const int iymin = max(1, iby*nbsizey);
            const int iymax = min(ncelly-1, (iby+1)*nbsizey);
            int iy = 0;
            /*---Loop over gridcells in block---*/
            for (iy=iymin; iy<iymax; ++iy)
            {
                const int ixmin = max(1, ibx*nbsizex);
                const int ixmax = min(ncellx-1, (ibx+1)*nbsizex);
                int ix = 0;
                for (ix=ixmin; ix<ixmax; ++ix)
                {
                    process_cell(ix, iy, ncellx, ncelly, vo, vi,
                                 an, as, ae, aw);
                }
            }
        }
    }
}

#if 0

/*===========================================================================*/
/*---Experimental code to use OpenMP tasks---*/

void solve_block_wavefronts_tasks(int ncellx, int ncelly,
                                  Float_t* __restrict__ vo,
                            const Float_t* __restrict__ vi,
                            const Float_t* __restrict__ an,
                            const Float_t* __restrict__ as,
                            const Float_t* __restrict__ ae,
                            const Float_t* __restrict__ aw)
{
    const int nbsizex = 20;
    const int nbsizey = 20;

    const int nbx = ceil_(ncellx, nbsizex);
    const int nby = ceil_(ncelly, nbsizey);

    const int nbwavefronts = nbx + nby - 1;
    int bwavefront = 0;

#pragma omp parallel
#pragma omp single
    {

    for (ibx=0; ibx<nbx; ++ibx)
    {
    for (iby=0; iby<nby; ++iby)
    {
#pragma omp task \
  depend(in: *ref(vo, nbsizex*(ibx-1), nbsizey* iby,    ncellx, ncelly)) \
  depend(in: *ref(vo, nbsizex* ibx,    nbsizey*(iby-1), ncellx, ncelly)) \
  depend(out:*ref(vo, nbsizex* ibx,    nbsizey* iby,    ncellx, ncelly))
    {
#pragma omp taskgroup if(MULTILEVEL_TASKING)
    {
            const int iymin = max(1, iby*nbsizey);
            const int iymax = min(ncelly-1, (iby+1)*nbsizey);
            int iy = 0;
            for (iy=iymin; iy<iymax; ++iy)
            {
                const int ixmin = max(1, ibx*nbsizex);
                const int ixmax = min(ncellx-1, (ibx+1)*nbsizex);
                int ix = 0;
                for (ix=ixmin; ix<ixmax; ++ix)
                {
/*
May want to experiment with different action based on device used at runtime:
  - tasks vs. parallel for - via "if" clause
  - block sizes
*/
#pragma omp task \
  depend(in: *ref(vi, ix-1, iy,   ncellx, ncelly)) \
  depend(in: *ref(vi, ix,   iy-1, ncellx, ncelly)) \
  depend(out:*ref(vi, ix,   iy,   ncellx, ncelly)) \
  if(MULTILEVEL_TASKING)
    {
                    process_cell(ix, iy, ncellx, ncelly, vo, vi,
                                 an, as, ae, aw);
    } /*---task---*/
                }
            }
    } /*---taskgroup---*/

    } /*---task---*/

    } /*---for---*/
    } /*---for---*/

    } /*---single---*/
}


#endif

enum{ solve_method_standard = 1,
      solve_method_wavefronts = 2,
      solve_method_block_wavefronts = 3 };

/*===========================================================================*/
/*---Driver---*/

int main( int argc, char** argv )
{
    /*---Settings---*/

    const int niterations = 10;
    const int ncellx = 2000;
    const int ncelly = 2000;

    const int solve_method = solve_method_block_wavefronts;

    /*
    const int solve_method = solve_method_standard;
    const int solve_method = solve_method_wavefronts;
    const int solve_method = solve_method_block_wavefronts;
    */

    /*---Allocations---*/

    Float_t* __restrict__ v1 = malloc(ncellx*ncelly*sizeof(Float_t));
    Float_t* __restrict__ v2 = malloc(ncellx*ncelly*sizeof(Float_t));
    Float_t* __restrict__ an = malloc(ncellx*ncelly*sizeof(Float_t));
    Float_t* __restrict__ as = malloc(ncellx*ncelly*sizeof(Float_t));
    Float_t* __restrict__ ae = malloc(ncellx*ncelly*sizeof(Float_t));
    Float_t* __restrict__ aw = malloc(ncellx*ncelly*sizeof(Float_t));

    /*---Initializations: interior---*/

    int iy = 0;
    int ix = 0;
    for (iy=1; iy<ncelly-1; ++iy)
    {
        for (ix=1; ix<ncellx-1; ++ix)
        {
            /*---Initial guess is 1 on interior, 0 on bdry---*/
            /*---Assume right hand side is 0---*/
            *ref(v1, ix, iy, ncellx, ncelly) = 1.;
            *ref(v2, ix, iy, ncellx, ncelly) = 1.;
            /*---Constant coefficient Laplacian operator---*/
            *ref(an, ix, iy, ncellx, ncelly) = .25;
            *ref(as, ix, iy, ncellx, ncelly) = .25;
            *ref(ae, ix, iy, ncellx, ncelly) = .25;
            *ref(aw, ix, iy, ncellx, ncelly) = .25;
        }
    }

    /*---Initializations: boundary---*/
    /*---Assume Dirichlet zero boundaries---*/

    for (iy=0; iy<ncelly; ++iy)
    {
        *ref(v1, 0,        iy, ncellx, ncelly) = 0.;
        *ref(v2, 0,        iy, ncellx, ncelly) = 0.;
        *ref(v1, ncellx-1, iy, ncellx, ncelly) = 0.;
        *ref(v2, ncellx-1, iy, ncellx, ncelly) = 0.;
    }

    for (ix=0; ix<ncellx; ++ix)
    {
        *ref(v1, ix, 0,        ncellx, ncelly) = 0.;
        *ref(v2, ix, 0,        ncellx, ncelly) = 0.;
        *ref(v1, ix, ncelly-1, ncellx, ncelly) = 0.;
        *ref(v2, ix, ncelly-1, ncellx, ncelly) = 0.;
    }

    /*---Iteration loop---*/

    double tbeg = get_time();

    int iteration = 0;
    for (iteration=0; iteration<niterations; ++iteration)
    {
        const Float_t* const __restrict__ vi = iteration%2 ? v1 : v2;
              Float_t* const __restrict__ vo = iteration%2 ? v2 : v1;

        if (solve_method==solve_method_standard)
            solve_standard(ncellx, ncelly, vo, vi, an, as, ae, aw);

        if (solve_method==solve_method_wavefronts)
            solve_wavefronts(ncellx, ncelly, vo, vi, an, as, ae, aw);

        if (solve_method==solve_method_block_wavefronts)
            solve_block_wavefronts(ncellx, ncelly, vo, vi, an, as, ae, aw);
    }

    /*---Finish---*/

    double time = get_time() - tbeg;

    const Float_t* __restrict__ vfinal = niterations%2 ? v1 : v2;

    Float_t sumsq = 0.;

    for (iy=0; iy<ncelly; ++iy)
    {
        for (ix=0; ix<ncellx; ++ix)
        {
            Float_t value = *const_ref(vfinal, ix, iy, ncellx, ncelly);
            sumsq += value * value;
        }
    }

    printf("Cells %e sumsq %.16e time %.6f\n", (double)ncellx*ncelly,
           (double)sumsq, time);

    /*---Deallocations---*/

    free(v1);
    free(v2);
    free(an);
    free(as);
    free(ae);
    free(aw);
}

/*===========================================================================*/
