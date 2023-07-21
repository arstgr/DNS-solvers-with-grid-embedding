/********************************************************************************
 * Amirreza Rastegari							                             	*
 * arstgri@gmail.com						                             		*
 * This code is used to calculate the 2 point velocity correlations in parallel *
 * Formulation:								                                 	*
 * Mathiu & Scott, An introduction to turbulent flows, Cambridge university	    *
 * press, 2000, pp. 62								                            *
 ********************************************************************************/
#include "definitions.h"

POINTER reading(int xdc, int ydc, int zdc, int xdf, int ydf, int zdf, POINTER V, PDATA pd, long ts)
{
	int i,j,k,nn,mm;
	int a,q,b,kk;
	FILE *sv;
	int exty,extx;
	char fn[40];
	long cntr;
	MPI_File fh;
        MPI_Status status;
        MPI_Offset offset;

	cntr=ts;
	if (pd.myrank==0)
		fprintf(stderr,"TRYING TO READ DATA FROM FILES, I'm root\n");
	sprintf(fn,"../vel-c.%.4d",cntr);
//        fprintf(stderr,"reading data from %s myrank is%d\n",fn,pd.myrank);
	MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);

	b=pd.myrank;

	offset=b*xdc*3*ydc*zdc*sizeof(double);
        MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,V.velc,(xdc*ydc*zdc*3),MPI_DOUBLE,&status);
	
	MPI_File_close(&fh);
	
	sprintf(fn,"../vel-slf.%.4d",cntr);
//        fprintf(stderr,"reading data from %s myrank is%d\n",fn,pd.myrank);
	MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
	offset=b*xdf*3*ydf*zdf*sizeof(double);
        MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,V.velfl,(xdf*ydf*zdf*3),MPI_DOUBLE,&status);
	MPI_File_close(&fh);
	
	sprintf(fn,"../vel-suf.%.4d",cntr);
//        fprintf(stderr,"reading data from %s myrank is%d\n",fn,pd.myrank);
	MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
	offset=b*xdf*3*ydf*zdf*sizeof(double);
        MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,V.velfu,(xdf*ydf*zdf*3),MPI_DOUBLE,&status);
	MPI_File_close(&fh);

//	fprintf(stderr,"reading finished i'm %d\n",pd.myrank);
return V;
}