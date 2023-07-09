/************************************************************************
 * By Amirreza Rastegari                                                *
 *   arstgri@gmail.com                                                  *
 *                                                                      *
 * Assumes a Multigrid data input                                       *
 * calculates rms, vorticity and mean velocities                        *
 * Calculates the energy budget2D                                       *
 * Skewness and Kurtosis                                                *
 * Assumes a 2D mean velocity i.e. U(y,z)                               *
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include "mpi.h"


#define zdv  197
#define ydv  256
#define xdv  512
#define CFL 0.1
#define Rebs 7200.
#define ub 0.0866034

/**************************************************************************
 * fine grid's parameters                                                 *
 * gratio : dx^c/dx^f ratio of the coarse to fine grid's resolution       *
 * ***********************************************************************/
#define gratio 4

// Number of z grid point on fine grid 
// a +1 is added inside the code to z for the wall points
#define zds 56


// This is the factor that defines the vorticity, either 1. or 0.5
#define fac 1.

#define TSTR 920
#define TEND 1220
#define DT 1

#define Fs(a,i,j,k,xd,yd,zd) (*(a+(i*yd*zd)+(j*zd)+k))
#define Fb(a,i,j,k,b,xd,yd,zd,d) (*(a+((i+2)*yd*zd*d)+(j*zd*d)+(k*d)+b))
#define Ts(i,j,k,xd,yd,zd) ((i*yd*zd)+(j*zd)+k)

// Definition of Derivatives
// F: Forward difference
// B: Backward difference

//1st ord
#define DDZ1F(a,k) ((*(a+k+1))-(*(a+k)))
#define DDZ1B(a,k) ((*(a+k))-(*(a+k-1)))
#define D2DZ1F(a,k) ((1.*(*(a+k))-2.*(*(a+k+1))+1.*(*(a+k+2)))/1.)
#define D2DZ1B(a,k) ((-1.*(*(a+k))+2.*(*(a+k-1))-1.*(*(a+k-2)))/1.)

//2nd ord
#define DDZ2F(a,k) ((-3.*(*(a+k))+4.*(*(a+k+1))-(*(a+k+2)))*0.5)
#define DDZ2B(a,k) ((3.*(*(a+k))-4.*(*(a+k-1))+(*(a+k-2)))*0.5)
#define DDZ2C(a,k) (((*(a+k+1))-(*(a+k-1)))*0.5)
#define D2DZ2F(a,k) ((2.*(*(a+k))-5.*(*(a+k+1))+4.*(*(a+k+2))-(*(a+k+3)))/1.)
#define D2DZ2B(a,k) ((-2.*(*(a+k))+5.*(*(a+k-1))-4.*(*(a+k-2))+(*(a+k-3)))/1.)

//3rd ord
#define DDZ3F(a,k) ((-11.*(*(a+k))+18.*(*(a+k+1))-9.*(*(a+k+2))+2.*(*(a+k+3)))/6.)
#define DDZ3B(a,k) ((11.*(*(a+k))-18.*(*(a+k-1))+9.*(*(a+k-2))-2.*(*(a+k-3)))/6.)
#define D2DZ3F(a,k) ((35.*(*(a+k))-108.*(*(a+k+1))+114.*(*(a+k+2))-56.*(*(a+k+3))+11.*(*(a+k+4)))/12.)
#define D2DZ3B(a,k) ((-35.*(*(a+k))+108.*(*(a+k-1))-114.*(*(a+k-2))+56.*(*(a+k-3))-11.*(*(a+k-4)))/12.)

//4th ord
#define DDZ4F(a,k) ((-25.*(*(a+k))+48.*(*(a+k+1))-36.*(*(a+k+2))+(48./3.)*(*(a+k+3))-3.*(*(a+k+4)))/12.)
#define DDZ4B(a,k) ((25.*(*(a+k))-48.*(*(a+k-1))+36.*(*(a+k-2))-(48./3.)*(*(a+k-3))+3.*(*(a+k-4)))/12.)
#define DDZ4C(a,k) ((2./3.)*((*(a+k+1))-(*(a+k-1)))-(1./12.)*((*(a+k+2))-(*(a+k-2))))
#define D2DZ4F(a,k) ((15.*(*(a+k))-(154./3.)*(*(a+k+1))+(214./3.)*(*(a+k+2))-52.*(*(a+k+3))+(61./3.)*(*(a+k+4))-(10./3.)*(*(a+k+5)))/4.)
#define D2DZ4B(a,k) ((-15.*(*(a+k))+(154./3.)*(*(a+k-1))-(214./3.)*(*(a+k-2))+52.*(*(a+k-3))-(61./3.)*(*(a+k-4))+(10./3.)*(*(a+k-5)))/4.)

//5th ord
#define DDZ5F(a,k) ((-137.*(*(a+k))+300.*(*(a+k+1))-300.*(*(a+k+2))+200.*(*(a+k+3))-75.*(*(a+k+4))+12.*(*(a+k+5)))/60.)
#define DDZ5B(a,k) ((137.*(*(a+k))-300.*(*(a+k-1))+300.*(*(a+k-2))-200.*(*(a+k-3))+75.*(*(a+k-4))-12.*(*(a+k-5)))/60.)

//non uniform l
#define DDZ2CNUl(a,k,r) ((1./(1.+r))*((1./r)*((*(a+k+1))-(*(a+k)))+r*((*(a+k))-(*(a+k-1)))))

//non-uniform u
#define DDZ2CNUu(a,k,r) ((1./(1.+r))*((1./r)*((*(a+k))-(*(a+k-1)))+r*((*(a+k+1))-(*(a+k)))))

#define cs 7
#define bs 12

#define MPDP MPI_DOUBLE

typedef double DP;
#define PI 3.1415

typedef struct pointer_to_arrays
{
	int * s;
	DP * velc;
	DP * denc;
	DP * velfu;
	DP * denfu;
	DP * velfl;
	DP * denfl;
	
	DP * vel;
	DP * den;
} POINTER;

typedef struct parallel_data
{
	int myrank;
	int numproc;
	int right;
	int left;
} PDATA;

extern DP w[19];
extern int e0[19];
extern int e1[19];
extern int e2[19];
extern DP E0[19];
extern DP E1[19];
extern DP E2[19];
extern int opposite[19];
