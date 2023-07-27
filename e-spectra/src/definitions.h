/* ****************************************************************************
 * by Amirreza Rastegari                                                      *
 * arstgri@gmail.com                                                          *
 *                                                                            *
 * To be used with postproc-spectra.c                                         *
 *                                                                            *
 * last time tested: 11/20/2014                                               *
 *                                                                            *
 * Computes 1-dimensional energy spectra from DNS of turbulent channel flow   *
 *                                                                            *
 * Inputs:                                                                    *
 *    xd: dimension of velocity array in streamwise direction                 *
 *    yd: dimension of velocity array in spanwise direction                   *
 *    zd: dimension of velocity array in wall-normal direction                *
 *    TSTR: starting time for the calculations                                *
 *    TEND: final time for the calculations                                   *
 *                                                                            *
 * Input files:                                                               *
 *    vel.time: velocity file with array indices running from                 *
 *    (u_x,u_y,u_z), then z index, then y index, then x index                 *
 *    example: vel.0001                                                       *
 *                                                                            *
 * Output files:                                                              *
 *    Euux.time, Evvx.time, Ewwx.time                                         *
 *    Euuy.time, Evvy.time, Ewwy.time                                         *
 *                                                                            *
 * This program uses the fast fourier transform, you need to load the         *
 * necessary modules before the compilation.                                  *
 * To compile on stampede:                                                    *
 * -> module load fftw3                                                       *
 * -> mpicc 1d-spectra.c -I$TACC_FFTW3_INC -L$TACC_FFTW3_LIB -lfftw3_mpi      *
 *    -lfftw3                                                                 *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fftw3.h>
#include <fftw3-mpi.h>
#include "mpi.h"

#define xdv 512
#define ydv 256
#define zdv 197

/**************************************************************************
 * fine grid's parameters                                                 *
 * gratio : dx^c/dx^f ratio of the coarse to fine grid's resolution       *
 * ***********************************************************************/
#define gratio 4

// Number of z grid point on fine grid 
// a +1 is added inside the code to z for the wall points
#define zds 56


#define TSTR 916
#define TEND 1216
#define DT 1

#define Fb(a,i,j,k,b,xd,yd,zd,d) (*(a+(i*yd*zd*d)+(j*zd*d)+(k*d)+b))

typedef double DP;
#define PI 4.*atan(1.0)

typedef struct pointer_to_arrays
{
	DP * velc;
	DP * velfl;
	DP * velfu;
} POINTER;

typedef struct parallel_data
{
	int myrank;
	int numproc;
	int left;
	int right;
} PDATA;

