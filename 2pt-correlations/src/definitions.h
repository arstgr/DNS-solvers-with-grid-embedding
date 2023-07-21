/********************************************************************************
 * Amirreza Rastegari							                             	*
 * arstgri@gmail.com						                             		*
 * This code is used to calculate the 2 point velocity correlations in parallel *
 * Formulation:								                                 	*
 * Mathiu & Scott, An introduction to turbulent flows, Cambridge university	    *
 * press, 2000, pp. 62								                            *
 ********************************************************************************/

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

#define TSTR 916
#define TEND 1216
#define DT 1

#define Fs(a,i,j,k,xd,yd,zd) (*(a+(i*yd*zd)+(j*zd)+k))
#define Fb(a,i,j,k,b,xd,yd,zd,d) (*(a+(i*yd*zd*d)+(j*zd*d)+(k*d)+b))
#define Fy(a,j,k,b,yd,zd,d) (*(a+(k*yd*d)+(j*d)+b))
#define Fx(a,i,k,b,xd,zd,d) (*(a+(k*xd*d)+(i*d)+b))

#define MPDP MPI_DOUBLE

typedef double DP;
#define PI 4.*atan(1.0)

typedef struct pointer_to_arrays
{
	DP * velfl;
	DP * velfu;
	DP * velc;
	
	DP * uavgc;
	DP * vavgc;
	DP * wavgc;
	
	DP * uavgfl;
	DP * vavgfl;
	DP * wavgfl;
	
	DP * uavgfu;
	DP * vavgfu;
	DP * wavgfu;
	
	DP * velrecvc;
	DP * velrecvfl;
	DP * velrecvfu;
	
	DP * xcorlc;
	DP * xcorlfl;
	DP * xcorlfu;
	
	DP * ycorlc;
	DP * ycorlfl;
	DP * ycorlfu;
} POINTER;

typedef struct parallel_data
{
	int myrank;
	int numproc;
	int right;
	int left;
} PDATA;