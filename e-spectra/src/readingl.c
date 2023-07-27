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

#include "definitions.h"

POINTER readingl(int xdf, int ydf, int zdf, POINTER V, PDATA pd, int ts)
{
	char fn[40];
	int cntr,b;
	MPI_File fh;
	MPI_Status status;
	MPI_Offset offset,off1,off2;

	cntr=ts;
	if (pd.myrank==0)
		fprintf(stderr,"TRYING TO READ DATA FROM FILES\n");
	sprintf(fn,"../vel-slf.%.4d",cntr);
	fprintf(stderr,"reading data from %s\n",fn);
	MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);

	b=pd.myrank;

	off1 = (MPI_Offset)(zdf);
	off2 = (MPI_Offset)(ydf)*off1;
	offset=((MPI_Offset)b)*((MPI_Offset)(xdf-1))*((MPI_Offset)3)*off2*sizeof(double);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.velfl),((xdf-1)*ydf*zdf*3),MPI_DOUBLE,&status);
	
	MPI_File_close(&fh);

  return V;
}