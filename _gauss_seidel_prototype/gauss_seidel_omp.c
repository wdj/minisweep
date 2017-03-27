/*---------------------------------------------------------------------------*/
/*! 
 * \file   gauss_seidel.c
 * \author Wayne Joubert
 * \date   Mon Jul 27 11:38:40 EDT 2015
 * \brief  Code for block Gauss-Seidel solve with block wavefront method.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

/*=============================================================================

This code performs block Gauss-Seidel iteration to solve a specific linear
system of equations.  It is known that the Gauss-Seidel calculation has
recursive dependencies which are amenable to a wavefront method for parallel
solves.

The problem solved here is a multiple-dof 5-point finite difference stencil
problem on a regular 2D grid.  The Gauss-Seidel ordering used is induced by
aStandard lexicographical ordering applied to the grid unknowns.

=============================================================================*/

#include "sys/time.h"
#include "unistd.h"
#include "omp.h>"

#include <assert.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

int opcount[16 * 24]; //Padding to avoid false cache sharing

#define BILLION 1000000000L

typedef double Float_t;

//int omp_get_thread_num(){
//  return 0;
//}

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

/*---------------------------------------------------------------------------*/

void get_block_config_param(int ncellx, int ncelly, int *nbsizex, int *nbsizey, int *nbx, int *nby){

	/*---Set block sizes---*/
#ifdef SIZE
#ifdef BLOCK_GRID_W
    *nbsizex = SIZE/BLOCK_GRID_W;
#else
    *nbsizex = 20;
#endif
#ifdef BLOCK_GRID_H
    *nbsizey = SIZE/BLOCK_GRID_H;
#else
    *nbsizey = 20;
#endif
#else
    *nbsizex = 20;
    *nbsizey = 20;
#endif

    *nbx = ceil_(ncellx, *nbsizex);
    *nby = ceil_(ncelly, *nbsizey);

    return;
}

/*===========================================================================*/
/*---Timer function---*/

double get_time()
{
    struct timeval tv;
    int i = gettimeofday (&tv, NULL);
    return ( (double) tv.tv_sec +
             (double) tv.tv_usec * 1.e-6 );
}

/*============================================================================*/
/*---Accessor functions to access elements of array with multidim index---*/
/*---For matrix entries---*/

Float_t* refm(Float_t* __restrict__ a, int iu, int ju, int ix, int iy,
              int nu, int ncellx, int ncelly)
{
    assert(iu >= 0 && iu < nu);
    assert(ju >= 0 && ju < nu);
    assert(ix >= 0 && ix < ncellx);
    assert(iy >= 0 && iy < ncelly);

    return &( a[ju+nu*(iu+nu*(ix+ncellx*iy))] );
}

/*---------------------------------------------------------------------------*/

const Float_t* const_refm(const Float_t* __restrict__ a,
                          int iu, int ju, int ix, int iy,
                          int nu, int ncellx, int ncelly)
{
    assert(iu >= 0 && iu < nu);
    assert(ju >= 0 && ju < nu);
    assert(ix >= 0 && ix < ncellx);
    assert(iy >= 0 && iy < ncelly);

    return &( a[ju+nu*(iu+nu*(ix+ncellx*iy))] );
}

/*===========================================================================*/
/*---Accessor functions to access elements of array with multidim index---*/
/*---For vector entries---*/

Float_t* refv(Float_t* __restrict__ a, int iu, int ix, int iy,
              int nu, int ncellx, int ncelly)
{
	//if (iu <0 || iu >= nu || ix <0 || ix >= ncellx || iy < 0 || iy >= ncelly){
	//	printf("iu=%d, ix=%d, iy = %d\n", iu, ix, iy);
	//}
	assert(iu >= 0 && iu < nu);
    assert(ix >= 0 && ix < ncellx);
    assert(iy >= 0 && iy < ncelly);

    return &( a[iu+nu*(ix+ncellx*iy)] );
}

/*---------------------------------------------------------------------------*/

const Float_t* const_refv(const Float_t* __restrict__ a, int iu, int ix, int iy,
                          int nu, int ncellx, int ncelly)
{
	//if (iu <0 || iu >= nu || ix <0 || ix >= ncellx || iy < 0 || iy >= ncelly){
	//	printf("iu=%d, ix=%d, iy = %d\n", iu, ix, iy);
	//}
	assert(iu >= 0 && iu < nu);
    assert(ix >= 0 && ix < ncellx);
    assert(iy >= 0 && iy < ncelly);

    return &( a[iu+nu*(ix+ncellx*iy)] );
}

/*===========================================================================*/
/*---5-point stencil operation applied at a gridcell---*/

void process_cell(int ix, int iy, int nu, int ncellx, int ncelly,
                        Float_t* __restrict__ vo,
                  const Float_t* __restrict__ vi,
                  const Float_t* __restrict__ an,
                  const Float_t* __restrict__ as,
                  const Float_t* __restrict__ ae,
                  const Float_t* __restrict__ aw)
{
    int iu = 0;
    int ju = 0;

    /*---South, east connections use "old" values, north, west use "new"---*/
    /*---The connection operation is matrix-vector product---*/

    //printf("cell (%d,%d) on thread %d\n", ix, iy, omp_get_thread_num());
    for (iu=0; iu<nu; ++iu)
    {
        *refv( vo, iu, ix, iy, nu, ncellx, ncelly) = 0.;
        for (ju=0; ju<nu; ++ju)
        {
            *refv(      vo, iu,     ix,   iy,   nu, ncellx, ncelly) +=

            *const_refm(aw, iu, ju, ix,   iy,   nu, ncellx, ncelly) *
            *const_refv(vo,     ju, ix-1, iy,   nu, ncellx, ncelly) +

            *const_refm(an, iu, ju, ix,   iy,   nu, ncellx, ncelly) *
            *const_refv(vo,     ju, ix,   iy-1, nu, ncellx, ncelly) +

            *const_refm(ae, iu, ju, ix,   iy,   nu, ncellx, ncelly) *
            *const_refv(vi,     ju, ix+1, iy,   nu, ncellx, ncelly) +

            *const_refm(as, iu, ju, ix,   iy,   nu, ncellx, ncelly) *
            *const_refv(vi,     ju, ix,   iy+1, nu, ncellx, ncelly);
        }
    }
    opcount[16 * omp_get_thread_num()]++;
}

/*===========================================================================*/
/*---Perform Gauss-Seidel sweep with standard order of operations---*/

void solve_standard(int nu, int ncellx, int ncelly,
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
            process_cell(ix, iy, nu, ncellx, ncelly, vo, vi, an, as, ae, aw);
        }
    }
}

/*===========================================================================*/
/*---Perform Gauss-Seidel sweep with wavefront order of operations---*/
/*---Note all required recursion dependencies are satisfied---*/

void solve_wavefronts(int nu, int ncellx, int ncelly,
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
#pragma omp parallel for
        for (ix=ixmin; ix<ixmax; ++ix)
        {
            const int iy = wavefront - ix;
            process_cell(ix, iy, nu, ncellx, ncelly, vo, vi, an, as, ae, aw);
        }
    }
}

/*===========================================================================*/
/*---Perform Gauss-Seidel sweep with block-wavefront order of operations---*/
/*---Note all required recursion dependencies are satisfied---*/

void solve_block_wavefronts(int nu, int ncellx, int ncelly,
                                  Float_t* __restrict__ vo,
                            const Float_t* __restrict__ vi,
                            const Float_t* __restrict__ an,
                            const Float_t* __restrict__ as,
                            const Float_t* __restrict__ ae,
                            const Float_t* __restrict__ aw)
{
    /*---Set block sizes---*/
    int nbsizex = 20;
    int nbsizey = 20;

    int nbx = ceil_(ncellx, nbsizex);
    int nby = ceil_(ncelly, nbsizey);

    get_block_config_param(ncellx,ncelly, &nbsizex, &nbsizey, &nbx, &nby);

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
            const int ixmin = max(1, ibx*nbsizex);
            const int ixmax = min(ncellx-1, (ibx+1)*nbsizex);
            //printf("Running Block (%d, %d) - - (%d,%d) -> (%d, %d)  on thrid %d \n", ibx, iby,
            //		ixmin, iymin, ixmax, iymax, omp_get_thread_num());

            int iy = 0;

            /*---Loop over gridcells in block---*/
            for (iy=iymin; iy<iymax; ++iy)
            {
                //const int ixmin = max(1, ibx*nbsizex);
                //const int ixmax = min(ncellx-1, (ibx+1)*nbsizex);
                int ix = 0;
                for (ix=ixmin; ix<ixmax; ++ix)
                {
                    process_cell(ix, iy, nu, ncellx, ncelly, vo, vi,
                                 an, as, ae, aw);
                }
            }
        }
    }
}

/*---------------------------------------------------------------------------*/

void solve_dep(int nu, int ncellx, int ncelly,
                                  Float_t* __restrict__ vo,
                            const Float_t* __restrict__ vi,
                            const Float_t* __restrict__ an,
                            const Float_t* __restrict__ as,
                            const Float_t* __restrict__ ae,
                            const Float_t* __restrict__ aw)
{
	int ixmin = 1;
	int ixmax = ncellx - 1;
	int iymin = 1;
	int iymax = ncelly-1;

    float **dataptr = (float **) vi; //Just need the memory, nothing is written

//	    printf("Block (%d,%d) - (%d,%d) -> (%d, %d) starting on thread %d\n",
//		   ibx, iby, ixmin, iymin, ixmax, iymax,  omp_get_thread_num());
#pragma omp parallel
#pragma omp single
    {
        //printf("Running Block (%d, %d) - - (%d,%d) -> (%d, %d) \n", ibx, iby,
        //		ixmin, iymin, ixmax, iymax);
        int j=0, i=0;
    	for (j = iymin; j < iymax; j++){
	       for (i = ixmin; i < ixmax; i++){
//        	       printf("Launching cell (%d, %d) from thread %d\n", i, j, omp_get_thread_num());
#pragma omp task  depend(in:dataptr[i][j-1]) depend(in:dataptr[i-1][j]) \
			  depend(out:dataptr[i][j])
	    	  process_cell(i, j, nu, ncellx, ncelly, vo, vi,
                             an, as, ae, aw);
	       }
       }
    } // end taskgroup
    //printf("Finished Block (%d, %d)\n", ibx, iby);
}

/*---------------------------------------------------------------------------*/

void process_block_dep(int ibx, int iby,int nbsizex, int nbsizey,
							int nu, int ncellx, int ncelly,
                                  Float_t* __restrict__ vo,
                            const Float_t* __restrict__ vi,
                            const Float_t* __restrict__ an,
                            const Float_t* __restrict__ as,
                            const Float_t* __restrict__ ae,
                            const Float_t* __restrict__ aw)
{
	int ixmin = max(1, ibx*nbsizex);
	int ixmax = min(ncellx-1, (ibx+1)*nbsizex);
	int iymin = max(1, iby*nbsizey);
	int iymax = min(ncelly-1, (iby+1)*nbsizey);

    float **dataptr = (float **) vi; //Just need the memory, nothing is written

//	    printf("Block (%d,%d) - (%d,%d) -> (%d, %d) starting on thread %d\n",
//		   ibx, iby, ixmin, iymin, ixmax, iymax,  omp_get_thread_num());
#pragma omp taskgroup
    {
        //printf("Running Block (%d, %d) - - (%d,%d) -> (%d, %d) \n", ibx, iby,
        //		ixmin, iymin, ixmax, iymax);
        int i=0, j=0;
    	for (j = iymin; j < iymax; j++){
	       for (i = ixmin; i < ixmax; i++){
//        	       printf("Launching cell (%d, %d) from thread %d\n", i, j, omp_get_thread_num());
#pragma omp task  depend(in:dataptr[i][j-1]) depend(in:dataptr[i-1][j]) \
			  depend(out:dataptr[i][j])
	    	  process_cell(i, j, nu, ncellx, ncelly, vo, vi,
                             an, as, ae, aw);
	       }
       }
    } // end taskgroup
    //printf("Finished Block (%d, %d)\n", ibx, iby);
}

/*---------------------------------------------------------------------------*/

void process_block_parloop(int ibx, int iby,int nbsizex, int nbsizey,
							int nu, int ncellx, int ncelly,
                                  Float_t* __restrict__ vo,
                            const Float_t* __restrict__ vi,
                            const Float_t* __restrict__ an,
                            const Float_t* __restrict__ as,
                            const Float_t* __restrict__ ae,
                            const Float_t* __restrict__ aw)
{
	int ixmin = max(1, ibx*nbsizex);
	int ixmax = min(ncellx-1, (ibx+1)*nbsizex);
	int iymin = max(1, iby*nbsizey);
	int iymax = min(ncelly-1, (iby+1)*nbsizey);

    printf("Running Block (%d, %d) - - (%d,%d) -> (%d, %d)  on thrid %d \n", ibx, iby,
    		ixmin, iymin, ixmax, iymax, omp_get_thread_num());
	int wavefront = 0;
	for (wavefront = ixmin+iymin; wavefront <= ixmax+iymax-1; ++wavefront){

	    const int ix_low = max(ixmin, wavefront - iymax + 1);
	    const int ix_high = min(ixmax - 1, wavefront - iymin);

        /*---Loop over gridcells in wavefront---*/
	    int ix = 0;
#pragma omp parallel for
        for (ix=ix_low; ix <= ix_high; ++ix)
        {
            const int iy = wavefront - ix;
            printf("%d  %d  ", ix, iy);
            process_cell(ix, iy, nu, ncellx, ncelly, vo, vi, an, as, ae, aw);
        }
        printf("\n");
    }
    //printf("Finished Block (%d, %d)\n", ibx, iby);
}

/*---------------------------------------------------------------------------*/

void solve_block_dep_dep(int nu, int ncellx, int ncelly,
                                  Float_t* __restrict__ vo,
                            const Float_t* __restrict__ vi,
                            const Float_t* __restrict__ an,
                            const Float_t* __restrict__ as,
                            const Float_t* __restrict__ ae,
                            const Float_t* __restrict__ aw)
{

	int nbsizex = 0, nbsizey=0, nbx=0, nby=0;

    get_block_config_param(ncellx, ncelly, &nbsizex, &nbsizey, &nbx, &nby);

    float **gridptr = (float **) vo; //Just need the memory, nothing is written

    //omp_set_nested(1); //Needed for parallel loop in blocks

#pragma omp parallel
#pragma omp single
{
    int ibx = 0, iby=0;
    for (ibx = 0; ibx < nbx; ibx++){
        for (iby = 0; iby < nby; iby++){
        //printf("submitting  block (%d,%d) thrid=%d\n", ibx, iby, omp_get_thread_num()); 
#pragma omp task depend(in:gridptr[ibx][iby-1]) \
	depend(in:gridptr[ibx-1][iby]) depend(out:gridptr[ibx][iby])
	    {
	    	process_block_dep(ibx, iby, nbsizex, nbsizey,
	    			nu, ncellx, ncelly, vo, vi, an, as, ae, aw);
	    } // end task
        }
    }
} //end omp single
}

/*---------------------------------------------------------------------------*/

void solve_block_dep_parloop(int nu, int ncellx, int ncelly,
                                  Float_t* __restrict__ vo,
                            const Float_t* __restrict__ vi,
                            const Float_t* __restrict__ an,
                            const Float_t* __restrict__ as,
                            const Float_t* __restrict__ ae,
                            const Float_t* __restrict__ aw)
{

	int nbsizex = 0, nbsizey=0, nbx=0, nby=0;

    get_block_config_param(ncellx, ncelly, &nbsizex, &nbsizey, &nbx, &nby);

    float **gridptr = (float **) vo; //Just need the memory, nothing is written

    //omp_set_nested(1); //Needed for parallel loop in blocks

#pragma omp parallel
#pragma omp single
{
    int ibx = 0, iby = 0;
    for (ibx = 0; ibx < nbx; ibx++){
        for (iby = 0; iby < nby; iby++){
        //printf("submitting  block (%d,%d) thrid=%d\n", ibx, iby, omp_get_thread_num());
#pragma omp task depend(in:gridptr[ibx][iby-1]) \
	depend(in:gridptr[ibx-1][iby]) depend(out:gridptr[ibx][iby])
	    {
	    	process_block_parloop(ibx, iby, nbsizex, nbsizey,
	    			nu, ncellx, ncelly, vo, vi, an, as, ae, aw);
	    } // end task
        }
    }
} //end omp single
}

/*---------------------------------------------------------------------------*/

void solve_block_parloop_dep(int nu, int ncellx, int ncelly,
                                  Float_t* __restrict__ vo,
                            const Float_t* __restrict__ vi,
                            const Float_t* __restrict__ an,
                            const Float_t* __restrict__ as,
                            const Float_t* __restrict__ ae,
                            const Float_t* __restrict__ aw)
{
	int nbsizex = 0, nbsizey=0, nbx=0, nby=0;

    get_block_config_param(ncellx, ncelly, &nbsizex, &nbsizey, &nbx, &nby);

    const int nbwavefronts = nbx + nby - 1;
    int bwavefront = 0;

    /*---Loop over block wavefronts---*/
    //omp_set_nested(1); //Needed for parallel loop in blocks
    for (bwavefront=0; bwavefront<nbwavefronts; ++bwavefront){
        const int ibxmin = max(0, bwavefront-nby);
        const int ibxmax = min(nbx, bwavefront+1);
        /*---Loop over blocks in wavefront---*/
        int ibx = 0;
#pragma omp parallel for
        for (ibx=ibxmin; ibx<ibxmax; ++ibx){
            const int iby = bwavefront - ibx;
	    	process_block_dep(ibx, iby, nbsizex, nbsizey,
	    			nu, ncellx, ncelly, vo, vi, an, as, ae, aw);
        }
    }
}

/*---------------------------------------------------------------------------*/

void solve_block_parloop_parloop(int nu, int ncellx, int ncelly,
                                  Float_t* __restrict__ vo,
                            const Float_t* __restrict__ vi,
                            const Float_t* __restrict__ an,
                            const Float_t* __restrict__ as,
                            const Float_t* __restrict__ ae,
                            const Float_t* __restrict__ aw)
{
	int nbsizex = 0, nbsizey=0, nbx=0, nby=0;

    get_block_config_param(ncellx, ncelly, &nbsizex, &nbsizey, &nbx, &nby);

    const int nbwavefronts = nbx + nby - 1;
    int bwavefront = 0;

    /*---Loop over block wavefronts---*/
    //omp_set_nested(1); //Needed for parallel loop in blocks
    for (bwavefront=0; bwavefront<nbwavefronts; ++bwavefront){
        const int ibxmin = max(0, bwavefront-nby);
        const int ibxmax = min(nbx, bwavefront+1);
        /*---Loop over blocks in wavefront---*/
        int ibx = 0;
#pragma omp parallel for
        for (ibx=ibxmin; ibx<ibxmax; ++ibx){
            const int iby = bwavefront - ibx;
#pragma omp task untied
	    	process_block_parloop(ibx, iby, nbsizex, nbsizey,
	    			nu, ncellx, ncelly, vo, vi, an, as, ae, aw);
        }
    }
}

enum SOLVE_METHOD {
	  solve_method_standard = 1,
      solve_method_parloop = 2,
      solve_method_block_parloop = 3,
      solve_method_block_parloop_dep = 4,
	  solve_method_block_parloop_parloop = 5,
	  solve_method_block_dep_dep = 6,
	  solve_method_block_dep_parloop = 7,
	  solve_method_dep = 8,
	  solve_method_block_wavefronts,           //Sequential block wavefront
};

/*===========================================================================*/
/*---Driver---*/

int main( int argc, char** argv )
{
    /*---Settings---*/

    const int niterations = 1;
#ifdef SIZE
    const int ncellx = SIZE;
    const int ncelly = SIZE;
#else
    const int ncellx = 200;
    const int ncelly = 200;
#endif

#ifdef SIZE_U
    const int nu = SIZE_U;
#else
    const int nu = 10;
#endif

     enum SOLVE_METHOD solve_method = solve_method_standard;

     void ( *solve_fun)(int nu, int ncellx, int ncelly,
                               Float_t* __restrict__ vo,
                         const Float_t* __restrict__ vi,
                         const Float_t* __restrict__ an,
                         const Float_t* __restrict__ as,
                         const Float_t* __restrict__ ae,
                         const Float_t* __restrict__ aw);
    
     if (argc < 2){
     	printf ("Format: %s STANDARD | PARLOOP | DEP | BLOCK_PARLOOP_DEP | BLOCK_PARLOOP_PARLOOP | BLOCK_DEP_PARLOOP | BLOCK_DEP_DEP \n", argv[0]);
     	exit(1);
     }

     solve_fun = solve_wavefronts;
     if (strcmp(argv[1], "STANDARD") == 0){
     	solve_method = solve_method_standard;
     	solve_fun = solve_standard;
     }
     else if (strcmp(argv[1], "PARLOOP") == 0){
     	solve_method = solve_method_parloop;
     	solve_fun = solve_wavefronts;
     }
     else if (strcmp(argv[1], "DEP") == 0){
     	solve_method = solve_method_dep;
     	solve_fun = solve_dep;
     }
     else if (strcmp(argv[1], "BLOCK_PARLOOP_DEP") == 0){
     	solve_method = solve_method_block_parloop_dep;
     	solve_fun = solve_block_parloop_dep;
     }
     else if (strcmp(argv[1], "BLOCK_PARLOOP_PARLOOP") == 0){
     	solve_method = solve_method_block_parloop_parloop;
     	solve_fun = solve_block_parloop_parloop;
     }
     else if (strcmp(argv[1], "BLOCK_DEP_PARLOOP") == 0){
     	solve_method = solve_method_block_dep_parloop;
     	solve_fun = solve_block_dep_parloop;
     }
     else if (strcmp(argv[1], "BLOCK_DEP_DEP") == 0){
     	solve_method = solve_method_block_dep_dep;
     	solve_fun = solve_block_dep_dep;
     }
     else if (strcmp(argv[1], "BLOCK_WAVEFRONT") == 0){
     	solve_method = solve_method_block_wavefronts; //Sequential block wavefront
     	solve_fun = solve_block_wavefronts;
     }
     else {
     	printf("Invalid solution method %s\n", argv[1]);
     	exit(1);
     }
     int io = 0;
     for (io = 0; io < 24; io++){
    	 opcount[io] = 0;
     }


    /*
    const int solve_method = solve_method_standard;
    const int solve_method = solve_method_wavefronts;
    const int solve_method = solve_method_block_wavefronts;
    */

    /*---Allocations---*/

    Float_t* __restrict__ v1 = malloc(nu   *ncellx*ncelly*sizeof(Float_t));
    Float_t* __restrict__ v2 = malloc(nu   *ncellx*ncelly*sizeof(Float_t));
    Float_t* __restrict__ an = malloc(nu*nu*ncellx*ncelly*sizeof(Float_t));
    Float_t* __restrict__ as = malloc(nu*nu*ncellx*ncelly*sizeof(Float_t));
    Float_t* __restrict__ ae = malloc(nu*nu*ncellx*ncelly*sizeof(Float_t));
    Float_t* __restrict__ aw = malloc(nu*nu*ncellx*ncelly*sizeof(Float_t));

    /*---Initializations: interior---*/

    int iy = 0;
    int ix = 0;
    int iu = 0;
    int ju = 0;
    for (iy=1; iy<ncelly-1; ++iy)
    {
        for (ix=1; ix<ncellx-1; ++ix)
        {
            for (iu=0; iu<nu; ++iu)
            {
                for (ju=0; ju<nu; ++ju)
                {
                    *refm(an, iu, ju, ix, iy, nu, ncellx, ncelly) = 0.;
                    *refm(as, iu, ju, ix, iy, nu, ncellx, ncelly) = 0.;
                    *refm(ae, iu, ju, ix, iy, nu, ncellx, ncelly) = 0.;
                    *refm(aw, iu, ju, ix, iy, nu, ncellx, ncelly) = 0.;
                }
                /*---Initial guess is 1 on interior, 0 on bdry---*/
                /*---Assume right hand side is 0---*/
                *refv(v1, iu, ix, iy, nu, ncellx, ncelly) = 1.;
                *refv(v2, iu, ix, iy, nu, ncellx, ncelly) = 1.;
                /*---Constant coefficient Laplacian operator---*/
                *refm(an, iu, iu, ix, iy, nu, ncellx, ncelly) = .25;
                *refm(as, iu, iu, ix, iy, nu, ncellx, ncelly) = .25;
                *refm(ae, iu, iu, ix, iy, nu, ncellx, ncelly) = .25;
                *refm(aw, iu, iu, ix, iy, nu, ncellx, ncelly) = .25;
            }
        }
    }

    /*---Initializations: boundary---*/
    /*---Assume Dirichlet zero boundaries---*/

    for (iy=0; iy<ncelly; ++iy)
    {
        for (iu=0; iu<nu; ++iu)
        {
            *refv(v1, iu, 0,        iy, nu, ncellx, ncelly) = 0.;
            *refv(v2, iu, 0,        iy, nu, ncellx, ncelly) = 0.;
            *refv(v1, iu, ncellx-1, iy, nu, ncellx, ncelly) = 0.;
            *refv(v2, iu, ncellx-1, iy, nu, ncellx, ncelly) = 0.;
        }
    }

    for (ix=0; ix<ncellx; ++ix)
    {
        for (iu=0; iu<nu; ++iu)
        {
            *refv(v1, iu, ix, 0,        nu, ncellx, ncelly) = 0.;
            *refv(v2, iu, ix, 0,        nu, ncellx, ncelly) = 0.;
            *refv(v1, iu, ix, ncelly-1, nu, ncellx, ncelly) = 0.;
            *refv(v2, iu, ix, ncelly-1, nu, ncellx, ncelly) = 0.;
        }
    }

    /*---Iteration loop---*/


    struct timespec start, end;

    double tbeg = get_time();

    clock_gettime(CLOCK_MONOTONIC, &start);
    int iteration = 0;
    for (iteration=0; iteration<niterations; ++iteration)
    {
        const Float_t* const __restrict__ vi = iteration%2 ? v1 : v2;
              Float_t* const __restrict__ vo = iteration%2 ? v2 : v1;

        solve_fun(nu, ncellx, ncelly, vo, vi, an, as, ae, aw);
    }

    /*---Finish---*/

    clock_gettime(CLOCK_MONOTONIC, &end);

    double time = get_time() - tbeg;

    const Float_t* __restrict__ vfinal = niterations%2 ? v1 : v2;

    Float_t sumsq = 0.;

    for (iy=0; iy<ncelly; ++iy)
    {
        for (ix=0; ix<ncellx; ++ix)
        {
            for (iu=0; iu<nu; ++iu)
            {
                Float_t value = *const_refv(vfinal, iu, ix, iy,
                                            nu, ncellx, ncelly);
                sumsq += value * value;
            }
        }
    }

    double elapsed_time = (BILLION * (end.tv_sec - start.tv_sec) + end.tv_nsec - start.tv_nsec) * 1.e-9;

    printf("SIZE = %d  SIZE_U = %d  Unknowns %e sumsq %.12e time %.6e Sec\n",
    		ncellx, nu, (double)nu*(double)ncellx*(double)ncelly, (double)sumsq, elapsed_time);
    for (iy = 0; iy < 16 * 12; iy+=16){
    	printf ("%d  ", opcount[iy]);
    }
    printf("\n");

    /*---Deallocations---*/

    free(v1);
    free(v2);
    free(an);
    free(as);
    free(ae);
    free(aw);
}

/*===========================================================================*/
