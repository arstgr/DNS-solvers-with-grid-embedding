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
#include "definitions.h"

/* fluctuation field in a constant z plane near the lower wall 
*  Note: 0 =< position < zdf                                   */
POINTER uv_fluc(int xdf, int ydf, int zdf, POINTER V, PDATA pd, long ts, int position)
{
	int i,j,k,a,q,c,upn;
	double ufl[xdf-1][ydf],vfl[xdf-1][ydf];
	char fl[30];
	MPI_File fh;
	MPI_Offset offset;
	MPI_Status status;

	offset=pd.myrank*(xdf-1)*ydf*sizeof(double);
	
	k=position;
	for (i=0;i<(xdf-1);i++)
	{
		for (j=0;j<ydf;j++)
		{
			ufl[i][j] = Fb(V.velfl,i,j,k,0,xdf,ydf,zdf,3);
			vfl[i][j] = Fb(V.velfl,i,j,k,1,xdf,ydf,zdf,3);
		}
	}

	sprintf(fl,"ufl4.%.4d",((int)ts));
	MPI_File_open(MPI_COMM_WORLD,fl,MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_write_all(fh,&ufl[0][0],(xdf-1)*ydf,MPI_DOUBLE,&status);
	MPI_File_close(&fh);

	sprintf(fl,"vfl4.%.4d",((int)ts));
	MPI_File_open(MPI_COMM_WORLD,fl,MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_write_all(fh,&vfl[0][0],(xdf-1)*ydf,MPI_DOUBLE,&status);
	MPI_File_close(&fh);

return V;
}
