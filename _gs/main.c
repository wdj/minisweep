/*---------------------------------------------------------------------------*/
/*! 
 * \file   gauss_seidel.c
 * \author Wayne Joubert
 * \date   Mon Jul 27 11:38:40 EDT 2015
 * \brief  Code for applying Gauss-Seidel preconditioning with block
 *         wavefront method.  2D, lexicographical ordering, 5-point stencil.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <sys/time.h>

typedef double Float_t;

/*---------------------------------------------------------------------------*/

double get_time()
{
    struct timeval tv;
    int i = gettimeofday (&tv, NULL);
    return ( (double) tv.tv_sec +
             (double) tv.tv_usec * 1.e-6 );
}

/*---------------------------------------------------------------------------*/

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

/*---------------------------------------------------------------------------*/

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

/*---------------------------------------------------------------------------*/

void process_cell(int ix, int iy, int ncellx, int ncelly,
                        Float_t* __restrict__ vo,
                  const Float_t* __restrict__ vi,
                  const Float_t* __restrict__ an,
                  const Float_t* __restrict__ as,
                  const Float_t* __restrict__ ae,
                  const Float_t* __restrict__ aw)
{
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

/*---------------------------------------------------------------------------*/

void solve(int ncellx, int ncelly,
                 Float_t* __restrict__ vo,
           const Float_t* __restrict__ vi,
           const Float_t* __restrict__ an,
           const Float_t* __restrict__ as,
           const Float_t* __restrict__ ae,
           const Float_t* __restrict__ aw)
{
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

/*---------------------------------------------------------------------------*/

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

    for (wavefront=0; wavefront<nwavefronts; ++wavefront)
    {
        const int ixmin = max(1, wavefront-ncelly+2);
        const int ixmax = min(ncellx-1, wavefront);
        int ix = 0;
        for (ix=ixmin; ix<ixmax; ++ix)
        {
            const int iy = wavefront - ix;
            process_cell(ix, iy, ncellx, ncelly, vo, vi, an, as, ae, aw);
        }
    }
}

/*---------------------------------------------------------------------------*/

void solve_block_wavefronts(int ncellx, int ncelly,
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

    for (bwavefront=0; bwavefront<nbwavefronts; ++bwavefront)
    {
        const int ibxmin = max(0, bwavefront-nby);
        const int ibxmax = min(nbx, bwavefront+1);
        int ibx = 0;
        for (ibx=ibxmin; ibx<ibxmax; ++ibx)
        {
            const int iby = bwavefront - ibx;

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
                    process_cell(ix, iy, ncellx, ncelly, vo, vi,
                                 an, as, ae, aw);
                }
            }
        }
    }
}

/*---------------------------------------------------------------------------*/

int main( int argc, char** argv )
{
    /*---Settings---*/

#if 0 
    const int niterations = 1;
    const int ncellx = 5;
    const int ncelly = 5;
#endif
    const int niterations = 10;
    const int ncellx = 2000;
    const int ncelly = 2000;

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
            /*---Constant coefficient Laplacian---*/
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

#if 0
        solve(ncellx, ncelly, vo, vi, an, as, ae, aw);
        solve_wavefronts(ncellx, ncelly, vo, vi, an, as, ae, aw);
        solve_block_wavefronts(ncellx, ncelly, vo, vi, an, as, ae, aw);
#endif
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

/*---------------------------------------------------------------------------*/
