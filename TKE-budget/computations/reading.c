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

/* Reads the data in a row; first section by rank 0 and so on */
POINTER reading(int xdc, int ydc, int zdc, int xdf, int ydf, int zdf, POINTER V, PDATA pd, long ts, DP ubulk)
{
	int i,j,k,nn,mm;
	int a,q,b,kk;
	FILE *sv;
	int exty,extx;
	char fn[40];
	long cntr;
	double ux,uy,uz,dens;
	DP rl,rr,dq,kb;
	const int L[5]={2,8,9,13,14};
	const int R[5]={1,7,10,11,12};
	double *af,*sf;
	MPI_File fh;
        MPI_Status status;
        MPI_Offset offset;

	q=512*223*3*sizeof(double);
	cntr=ts;
	if (pd.myrank==0)
		fprintf(stderr,"TRYING TO READ DATA FROM FILES\n");
	sprintf(fn,"../vel-c.%.4d",cntr);
        fprintf(stderr,"reading data from %s\n",fn);
	MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);

	b=pd.myrank;

	offset=b*(xdc-1)*3*ydc*zdc*sizeof(double);
        MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.velc+2*(ydc*zdc*3)),((xdc-1)*ydc*zdc*3),MPI_DOUBLE,&status);

	MPI_Barrier(MPI_COMM_WORLD);	
	offset=(pd.left*(xdc-1)*3*ydc*zdc+(xdc-3)*3*ydc*zdc)*sizeof(double);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.velc),(ydc*zdc*3),MPI_DOUBLE,&status);
	
	MPI_Barrier(MPI_COMM_WORLD);	
	offset=(pd.left*(xdc-1)*3*ydc*zdc+(xdc-2)*3*ydc*zdc)*sizeof(double);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.velc+ydc*zdc*3),(ydc*zdc*3),MPI_DOUBLE,&status);

	MPI_Barrier(MPI_COMM_WORLD);
	offset=(pd.right*(xdc-1)*3*ydc*zdc)*sizeof(double);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.velc+(xdc+1)*ydc*zdc*3),(ydc*zdc*3),MPI_DOUBLE,&status);
	
	MPI_Barrier(MPI_COMM_WORLD);
	offset=(pd.right*(xdc-1)*3*ydc*zdc+ydc*zdc*3)*sizeof(double);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.velc+(xdc+2)*ydc*zdc*3),(ydc*zdc*3),MPI_DOUBLE,&status);
	
	MPI_File_close(&fh);
	
	sprintf(fn,"../vel-slf.%.4d",cntr);
        fprintf(stderr,"reading data from %s\n",fn);
	MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);

	b=pd.myrank;

	offset=b*(xdf-1)*3*ydf*zdf*sizeof(double);
        MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.velfl+2*(ydf*zdf*3)),((xdf-1)*ydf*zdf*3),MPI_DOUBLE,&status);

	MPI_Barrier(MPI_COMM_WORLD);	
	offset=(pd.left*(xdf-1)*3*ydf*zdf+(xdf-3)*3*ydf*zdf)*sizeof(double);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.velfl),(ydf*zdf*3),MPI_DOUBLE,&status);
	
	MPI_Barrier(MPI_COMM_WORLD);	
	offset=(pd.left*(xdf-1)*3*ydf*zdf+(xdf-2)*3*ydf*zdf)*sizeof(double);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.velfl+ydf*zdf*3),(ydf*zdf*3),MPI_DOUBLE,&status);

	MPI_Barrier(MPI_COMM_WORLD);
	offset=(pd.right*(xdf-1)*3*ydf*zdf)*sizeof(double);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.velfl+(xdf+1)*ydf*zdf*3),(ydf*zdf*3),MPI_DOUBLE,&status);
	
	MPI_Barrier(MPI_COMM_WORLD);
	offset=(pd.right*(xdf-1)*3*ydf*zdf+ydf*zdf*3)*sizeof(double);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.velfl+(xdf+2)*ydf*zdf*3),(ydf*zdf*3),MPI_DOUBLE,&status);
	
	MPI_File_close(&fh);
	
	sprintf(fn,"../vel-suf.%.4d",cntr);
        fprintf(stderr,"reading data from %s\n",fn);
	MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);

	b=pd.myrank;

	offset=b*(xdf-1)*3*ydf*zdf*sizeof(double);
        MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.velfu+2*(ydf*zdf*3)),((xdf-1)*ydf*zdf*3),MPI_DOUBLE,&status);

	MPI_Barrier(MPI_COMM_WORLD);	
	offset=(pd.left*(xdf-1)*3*ydf*zdf+(xdf-3)*3*ydf*zdf)*sizeof(double);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.velfu),(ydf*zdf*3),MPI_DOUBLE,&status);
	
	MPI_Barrier(MPI_COMM_WORLD);	
	offset=(pd.left*(xdf-1)*3*ydf*zdf+(xdf-2)*3*ydf*zdf)*sizeof(double);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.velfu+ydf*zdf*3),(ydf*zdf*3),MPI_DOUBLE,&status);

	MPI_Barrier(MPI_COMM_WORLD);
	offset=(pd.right*(xdf-1)*3*ydf*zdf)*sizeof(double);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.velfu+(xdf+1)*ydf*zdf*3),(ydf*zdf*3),MPI_DOUBLE,&status);
	
	MPI_Barrier(MPI_COMM_WORLD);
	offset=(pd.right*(xdf-1)*3*ydf*zdf+ydf*zdf*3)*sizeof(double);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.velfu+(xdf+2)*ydf*zdf*3),(ydf*zdf*3),MPI_DOUBLE,&status);
	
	MPI_File_close(&fh);

	sprintf(fn,"../den-c.%.4d",cntr);
        fprintf(stderr,"reading data from %s\n",fn);
	MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);

	b=pd.myrank;

	offset=b*(xdc-1)*ydc*zdc*sizeof(double);
        MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.denc+2*(ydc*zdc)),((xdc-1)*ydc*zdc),MPI_DOUBLE,&status);

	MPI_Barrier(MPI_COMM_WORLD);	
	offset=(pd.left*(xdc-1)*ydc*zdc+(xdc-3)*ydc*zdc)*sizeof(double);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.denc),(ydc*zdc),MPI_DOUBLE,&status);
	
	MPI_Barrier(MPI_COMM_WORLD);	
	offset=(pd.left*(xdc-1)*ydc*zdc+(xdc-2)*ydc*zdc)*sizeof(double);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.denc+ydc*zdc),(ydc*zdc),MPI_DOUBLE,&status);

	MPI_Barrier(MPI_COMM_WORLD);
	offset=(pd.right*(xdc-1)*ydc*zdc)*sizeof(double);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.denc+(xdc+1)*ydc*zdc),(ydc*zdc),MPI_DOUBLE,&status);
	
	MPI_Barrier(MPI_COMM_WORLD);
	offset=(pd.right*(xdc-1)*ydc*zdc+ydc*zdc)*sizeof(double);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.denc+(xdc+2)*ydc*zdc),(ydc*zdc),MPI_DOUBLE,&status);
	
	MPI_File_close(&fh);
	
	sprintf(fn,"../den-slf.%.4d",cntr);
        fprintf(stderr,"reading data from %s\n",fn);
	MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);

	b=pd.myrank;

	offset=b*(xdf-1)*ydf*zdf*sizeof(double);
        MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.denfl+2*(ydf*zdf)),((xdf-1)*ydf*zdf),MPI_DOUBLE,&status);

	MPI_Barrier(MPI_COMM_WORLD);	
	offset=(pd.left*(xdf-1)*ydf*zdf+(xdf-3)*ydf*zdf)*sizeof(double);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.denfl),(ydf*zdf),MPI_DOUBLE,&status);
	
	MPI_Barrier(MPI_COMM_WORLD);	
	offset=(pd.left*(xdf-1)*ydf*zdf+(xdf-2)*ydf*zdf)*sizeof(double);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.denfl+ydf*zdf),(ydf*zdf),MPI_DOUBLE,&status);

	MPI_Barrier(MPI_COMM_WORLD);
	offset=(pd.right*(xdf-1)*ydf*zdf)*sizeof(double);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.denfl+(xdf+1)*ydf*zdf),(ydf*zdf),MPI_DOUBLE,&status);
	
	MPI_Barrier(MPI_COMM_WORLD);
	offset=(pd.right*(xdf-1)*ydf*zdf+ydf*zdf)*sizeof(double);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.denfl+(xdf+2)*ydf*zdf),(ydf*zdf),MPI_DOUBLE,&status);
	
	MPI_File_close(&fh);
	
	sprintf(fn,"../den-suf.%.4d",cntr);
        fprintf(stderr,"reading data from %s\n",fn);
	MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);

	b=pd.myrank;

	offset=b*(xdf-1)*ydf*zdf*sizeof(double);
        MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.denfu+2*(ydf*zdf)),((xdf-1)*ydf*zdf),MPI_DOUBLE,&status);

	MPI_Barrier(MPI_COMM_WORLD);	
	offset=(pd.left*(xdf-1)*ydf*zdf+(xdf-3)*ydf*zdf)*sizeof(double);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.denfu),(ydf*zdf),MPI_DOUBLE,&status);
	
	MPI_Barrier(MPI_COMM_WORLD);	
	offset=(pd.left*(xdf-1)*ydf*zdf+(xdf-2)*ydf*zdf)*sizeof(double);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.denfu+ydf*zdf),(ydf*zdf),MPI_DOUBLE,&status);

	MPI_Barrier(MPI_COMM_WORLD);
	offset=(pd.right*(xdf-1)*ydf*zdf)*sizeof(double);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.denfu+(xdf+1)*ydf*zdf),(ydf*zdf),MPI_DOUBLE,&status);
	
	MPI_Barrier(MPI_COMM_WORLD);
	offset=(pd.right*(xdf-1)*ydf*zdf+ydf*zdf)*sizeof(double);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.denfu+(xdf+2)*ydf*zdf),(ydf*zdf),MPI_DOUBLE,&status);
	
	MPI_File_close(&fh);

return V;
}
