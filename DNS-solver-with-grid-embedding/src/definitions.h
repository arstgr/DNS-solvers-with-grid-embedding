/***************************************************************************
 *   Amirreza Rastegari                                                    *
 *   arstgri@gmail.com                                                     *
 *                                                                         *
 *                                                                         *
 *   Parallel LBM-BGK code with constant flow rate Q for DNS of turbulent  *
 *   channel flows                                                         *
 *   This program uses a new formulation for the constant flow rate        *
 *   implementation                                                        *
 *   A multigrid implementations is used for increasing the accuracy       *
 *   in the near wall reagion                                              *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <string.h>
#include "mpi.h"

/**************************************************************************
 *   Input Values of the program are listed here                          *
 *   Reb : The bulk Reynolds number of the flow Re_b = ubxh/nu : h = Lz/2 *
 *   ub  : The bulk flow velocity ub/c : c =dx/dt                         *
 *   zdv : Number of points in Z direction (Wall normall direction)       *
 *         a +2 is added to this value to consider the overlap region     *
 *   ydv : Number of points in Y direction (spanwise direction)           *
 *   xdv : Number of points in X direction (streamwise direction)         *
 *   CFL : The value of the CFL number for the program                    *
 *   TSTR: Time at which the computations start                           *
 *   TEND: Time at which computations finish                              *
 **************************************************************************/

#define Rebs 7200.000
#define ub 0.0866034
#define zdv  220
#define ydv  256
#define xdv  512
#define TSTR 400
#define TEND 800
#define CFL  0.1

/**************************************************************************
 * fine grid's parameters                                                 *
 * gratio : dx^c/dx^f ratio of the coarse to fine grid's resolution       *
 * ***********************************************************************/
#define gratio 4

// Number of z grid point on fine grid 
// a +1 is added inside the code to z for the wall points
#define zds 10

/**************************************************************************
 * gplussize: width of the gas grooves                                    *
 * wplussize: width of the no slip ridge                                  *
 *************************************************************************/
#define gyplussize 56
#define wyplussize 8
//if posts or spanwise ridges
#define gxplussize 4
#define wxplussize 4

/**************************************************************************
 * Defining the slip pattern on the walls                                 *
 * PSLIP: slip with posts                                                 *
 * RSLIP: slip with streamwise ridges                                     *
 * SSLIP: slip with transverse ridges                                     *
 **************************************************************************/
//#define PSLIP
#define RSLIP
//#define SSLIP

/**************************************************************************
 * Uses a 2D cartesian topology of a x b processors (x,y)                 *
 * Note: a x b = Number of Cores                                          *
 **************************************************************************/
#define XDIM 32
#define YDIM 16

/**************************************************************************
 * Filtering the small scales when going to the coarse grid               *
 * A gaussian box filter is used based on:                                *
 * Touil, H., Ricot, D. and Emmanuel, L. "Direct and Large-Eddy           *
 * Simulation of Turbulent Flows on Composite Multi-Resolution Grids by   *
 * the Lattice Boltzmann Method", 2013, French ministry of industry       *
 **************************************************************************/
#define NOFLTR
//#define FLTR
//#define BFLTR
//#define GFLTR

/**************************************************************************
 * Defining the number of iterations                                      *
 * This value must be used whenever a predetermined number of iterations  *
 * is desired, otherwise it must be commented out                         *
 * These are used in debug mode to detect the errors                      *
 **************************************************************************/
//#define ITER 1

/**************************************************************************
 * Input parameters for initialization from the lb-dns code with uniform  *
 * grid                                                                   *
 * xdu: Number of grid points in streamwise direction                     *
 * ydu: Number of grid points in spanwise direction                       *
 * zdu + 2: Number of grid points in wall normal direction                *
 **************************************************************************/
#define xdu 512
#define ydu 256
#define zdu 221

/**************************************************************************
 * Cahnnel Sizes                                                          *
 * Lx/h = 2xpi/alpha is the channel size in the streamwise direction      *
 * Ly/h = 2xpi/beta is the size in the spanwise direction                 *
 * Lz/h=2 is the size in wall normal direction                            *
 **************************************************************************/

#define ALPHA 0.339
#define BETA 0.678

/**************************************************************************
 * Code starts from here                                                  *
 * DON'T TOUCH                                                            *
 **************************************************************************/

#define STRIPE_COUNT "64"       /* must be an ascii string */
#define STRIPE_SIZE "1048576"    /* must be an ascii string */
#define CB_NODES_FCT "64"       /* max number of cores used in collective write operations */

#define Fs(a,i,j,k,xd,yd,zd) (*(a+(i*yd*zd)+(j*zd)+k))
#define Fb(a,i,j,k,b,xd,yd,zd,d) (*(a+(i*yd*zd*d)+(j*zd*d)+(k*d)+b))
//#define Fbi(a,i,j,k,b,xd,yd,zd,d) (*(a+((i+2)*yd*zd*d)+(j*zd*d)+(k*d)+b))
#define Fbi(a,i,j,t,b,xd,yd,td,d) (*(a+(t*(xd+5)*(yd+5)*d)+((i+2)*(yd+5)*d)+((j+2)*d)+b))
//#define Fb(a,i,j,k,b,xd,yd,zd,d) (*(a+(b*xd*yd*zd)+(i*yd*zd)+(j*zd)+k))
#define Ts(i,j,k,xd,yd,zd) ((i*yd*zd)+(j*zd)+k)

#define Vsp(a,i,j,k,b) (*(a+j*spx*(spz+1)*3+k*spx*3+i*3+b))
#define Vlb(a,i,j,k,b) (*(a+j*xdt*zd*3+i*zd*3+k*3+b))
#define Vspp(a,i,j,k) (*(a+j*spx*(spz+1)+k*spx+i))
#define Vlbp(a,i,j,k) (*(a+j*xdt*zd+i*zd+k))
#define VLB(a,i,j,k,b) (*(a+j*xd*zd*3+i*zd*3+k*3+b))
#define PLB(a,i,j,k) (*(a+j*xd*zd+i*zd+k))

#define N1(x,y,z) ((1./8.)*(1-x)*(1-y)*(1-z))
#define N2(x,y,z) ((1./8.)*(1+x)*(1-y)*(1-z))
#define N3(x,y,z) ((1./8.)*(1+x)*(1+y)*(1-z))
#define N4(x,y,z) ((1./8.)*(1-x)*(1+y)*(1-z))
#define N5(x,y,z) ((1./8.)*(1-x)*(1-y)*(1+z))
#define N6(x,y,z) ((1./8.)*(1+x)*(1-y)*(1+z))
#define N7(x,y,z) ((1./8.)*(1+x)*(1+y)*(1+z))
#define N8(x,y,z) ((1./8.)*(1-x)*(1+y)*(1+z))

#define csc 7
#define csf 7

#define MPDP MPI_DOUBLE

typedef double DP;
DP PI;

typedef struct pointer_to_arrays
{
	DP rho[3];
	DP force[3];
	int * s;
	DP * f;
	DP * ftemp;
	DP * rightbufs;
	DP * rightbufr;
	DP * leftbufs;
	DP * leftbufr;
	DP * upbufs;
	DP * upbufr;
	DP * dwbufs;
	DP * dwbufr;
	
	DP * urbufs;
	DP * dlbufr;
	DP * dlbufs;
	DP * urbufr;
	MPI_Request req[12];
// Fine grid top: u
	int * sus;
	DP * suf;
	DP * suftemp;
	DP * surightbufs;
	DP * surightbufr;
	DP * suleftbufs;
	DP * suleftbufr;
	DP * suupbufs;
	DP * suupbufr;
	DP * sudwbufs;
	DP * sudwbufr;
	
	DP * suurbufs;
	DP * sudlbufr;
	DP * sudlbufs;
	DP * suurbufr;
	MPI_Request sureq[12];
// Fine grid buttom: l
	int * sls;
	DP * slf;
	DP * slftemp;
	DP * slrightbufs;
	DP * slrightbufr;
	DP * slleftbufs;
	DP * slleftbufr;
	DP * slupbufs;
	DP * slupbufr;
	DP * sldwbufs;
	DP * sldwbufr;
	
	DP * slurbufs;
	DP * sldlbufr;
	DP * sldlbufs;
	DP * slurbufr;
	MPI_Request slreq[12];
// Interpolation
	DP * uin;
	DP * uinp;
	DP * uinpp;
	DP * lin;
	DP * linp;
	DP * linpp;

// Time interpolation
	DP * uint;
	DP * lint;
} POINTER;

typedef struct parallel_data
{
	int myrank;
	int numproc;
	int right;
	int left;
	int up;
	int dw;
	int ur;
	int ul;
	int dr;
	int dl;
} PDATA;

extern const DP w[19];
extern const DP BI[19];

extern int e0[19];
extern int e1[19];
extern int e2[19];
extern DP E0[19];
extern DP E1[19];
extern DP E2[19];
extern int opposite[19];
extern int mirror[19];

extern int mpi_errors,mpi_errors_length,mpi_errors_class;
extern char mpi_errors_message[MPI_MAX_ERROR_STRING];
extern MPI_Comm cart_grid;

#ifndef isnan
__inline int isnan(double var)
{
	return var!=var;
}
#endif