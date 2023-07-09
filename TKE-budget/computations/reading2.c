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

POINTER reading2(int xdc, int ydc, int zdc, int xdf, int ydf, int zdf, POINTER V, PDATA pd, int ts)
{
	int i,j,k,nn,mm;
	int a,q,b,kk;
	FILE *sv;
	int exty,extx;
	char fn[80];
	int cntr;
	MPI_File fh;
        MPI_Status status;
        MPI_Offset offset;
	
	cntr=ts;
	
	sprintf(fn,"../vel-c.%.4d",cntr); 
	MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
	MPI_File_set_view(fh,0,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read_at_all(fh,(pd.myrank*(xdc-1)*ydc*zdc*3),(V.velc+2*ydc*zdc*3),((xdc-1)*ydc*zdc*3),MPI_DOUBLE,&status);
	MPI_File_close(&fh);
	
	MPI_Get_count(&status,MPI_DOUBLE,&i);
	printf("rank=%d V.velc %d\n",pd.myrank,i);
	
	MPI_Sendrecv((V.velc+2*ydc*zdc*3),(2*ydc*zdc*3),MPI_DOUBLE,pd.left,1129,(V.velc+(xdc-1+2)*3*ydc*zdc),(2*ydc*zdc*3),MPI_DOUBLE,pd.right,1129,MPI_COMM_WORLD,&status);
	MPI_Sendrecv((V.velc+(xdc-1-2+2)*ydc*zdc*3),(2*ydc*zdc*3),MPI_DOUBLE,pd.right,1128,(V.velc),(2*ydc*zdc*3),MPI_DOUBLE,pd.left,1128,MPI_COMM_WORLD,&status);
	
	sprintf(fn,"../den-c.%.4d",cntr); 
	MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
	MPI_File_set_view(fh,0,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read_at_all(fh,(pd.myrank*(xdc-1)*ydc*zdc),(V.denc+2*ydc*zdc),((xdc-1)*ydc*zdc),MPI_DOUBLE,&status);
	MPI_File_close(&fh);
	
	MPI_Sendrecv((V.denc+2*ydc*zdc),(2*ydc*zdc),MPI_DOUBLE,pd.left,1129,(V.denc+(xdc-1+2)*ydc*zdc),(2*ydc*zdc),MPI_DOUBLE,pd.right,1129,MPI_COMM_WORLD,&status);
	MPI_Sendrecv((V.denc+(xdc-1-2+2)*ydc*zdc),(2*ydc*zdc),MPI_DOUBLE,pd.right,1128,(V.denc),(2*ydc*zdc),MPI_DOUBLE,pd.left,1128,MPI_COMM_WORLD,&status);
	
	sprintf(fn,"../vel-slf.%.4d",cntr); 
	MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
	MPI_File_set_view(fh,0,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read_at_all(fh,(pd.myrank*(xdf-1)*ydf*zdf*3),(V.velfl+2*ydf*zdf*3),((xdf-1)*ydf*zdf*3),MPI_DOUBLE,&status);
	MPI_File_close(&fh);
	
	MPI_Get_count(&status,MPI_DOUBLE,&i);
	printf("rank=%d V.velfl %d\n",pd.myrank,i);
	
	MPI_Sendrecv((V.velfl+2*ydf*zdf*3),(2*ydf*zdf*3),MPI_DOUBLE,pd.left,1129,(V.velfl+(xdf-1+2)*3*ydf*zdf),(2*ydf*zdf*3),MPI_DOUBLE,pd.right,1129,MPI_COMM_WORLD,&status);
	MPI_Sendrecv((V.velfl+(xdf-1-2+2)*ydf*zdf*3),(2*ydf*zdf*3),MPI_DOUBLE,pd.right,1128,(V.velfl),(2*ydf*zdf*3),MPI_DOUBLE,pd.left,1128,MPI_COMM_WORLD,&status);
	
	sprintf(fn,"../den-slf.%.4d",cntr); 
	MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
	MPI_File_set_view(fh,0,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read_at_all(fh,(pd.myrank*(xdf-1)*ydf*zdf),(V.denfl+2*ydf+zdf),((xdf-1)*ydf*zdf),MPI_DOUBLE,&status);
	MPI_File_close(&fh);
	
	MPI_Sendrecv((V.denfl+2*ydf*zdf),(2*ydf*zdf),MPI_DOUBLE,pd.left,1129,(V.denfl+(xdf-1+2)*ydf*zdf),(2*ydf*zdf),MPI_DOUBLE,pd.right,1129,MPI_COMM_WORLD,&status);
	MPI_Sendrecv((V.denfl+(xdf-1-2+2)*ydf*zdf),(2*ydf*zdf),MPI_DOUBLE,pd.right,1128,(V.denfl),(2*ydf*zdf),MPI_DOUBLE,pd.left,1128,MPI_COMM_WORLD,&status);
	
	sprintf(fn,"../vel-suf.%.4d",cntr); 
	MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
	MPI_File_set_view(fh,0,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read_at_all(fh,(pd.myrank*(xdf-1)*ydf*zdf*3),(V.velfu+2*ydf*zdf*3),((xdf-1)*ydf*zdf*3),MPI_DOUBLE,&status);
	MPI_File_close(&fh);
	
	MPI_Get_count(&status,MPI_DOUBLE,&i);
	printf("rank=%d V.velfu %d\n",pd.myrank,i);
	
	MPI_Sendrecv((V.velfu+2*ydf*zdf*3),(2*ydf*zdf*3),MPI_DOUBLE,pd.left,1129,(V.velfu+(xdf-1+2)*3*ydf*zdf),(2*ydf*zdf*3),MPI_DOUBLE,pd.right,1129,MPI_COMM_WORLD,&status);
	MPI_Sendrecv((V.velfu+(xdf-1-2+2)*ydf*zdf*3),(2*ydf*zdf*3),MPI_DOUBLE,pd.right,1128,(V.velfu),(2*ydf*zdf*3),MPI_DOUBLE,pd.left,1128,MPI_COMM_WORLD,&status);
	
	sprintf(fn,"../den-suf.%.4d",cntr); 
	MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
	MPI_File_set_view(fh,0,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read_at_all(fh,(pd.myrank*(xdf-1)*ydf*zdf),(V.denfu+2*ydf*zdf),((xdf-1)*ydf*zdf),MPI_DOUBLE,&status);
	MPI_File_close(&fh);
	
	MPI_Sendrecv((V.denfu+2*ydf*zdf),(2*ydf*zdf),MPI_DOUBLE,pd.left,1129,(V.denfu+(xdf-1+2)*ydf*zdf),(2*ydf*zdf),MPI_DOUBLE,pd.right,1129,MPI_COMM_WORLD,&status);
	MPI_Sendrecv((V.denfu+(xdf-1-2+2)*ydf*zdf),(2*ydf*zdf),MPI_DOUBLE,pd.right,1128,(V.denfu),(2*ydf*zdf),MPI_DOUBLE,pd.left,1128,MPI_COMM_WORLD,&status);
  
  return V;
}