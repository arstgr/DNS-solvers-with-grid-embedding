/********************************************************************************
 * Amirreza Rastegari							                             	*
 * arstgri@gmail.com						                             		*
 * This code is used to calculate the 2 point velocity correlations in parallel *
 * Formulation:								                                 	*
 * Mathiu & Scott, An introduction to turbulent flows, Cambridge university	    *
 * press, 2000, pp. 62								                            *
 ********************************************************************************/
#include "definitions.h"

POINTER xcorrel(int xdc,int ydc,int zdc, int xdf,int ydf,int zdf, POINTER V,PDATA pd)
{
	int i,j,k,a,b,c;
	double *xcorc,*xcorfl,*xcorfu;
	int ito,iito;
	int num;
	MPI_Status status;
	int stag=100,rtag=100;

	memcpy(V.velrecvc,V.velc,xdc*ydc*zdc*3*sizeof(double));
	memcpy(V.velrecvfl,V.velfl,xdf*ydf*zdf*3*sizeof(double));
	memcpy(V.velrecvfu,V.velfu,xdf*ydf*zdf*3*sizeof(double));
	
	xcorc = (double *)calloc(3*xdv*zdc,sizeof(double));
	xcorfl = (double *)calloc(3*xdv*gratio*zdf,sizeof(double));
	xcorfu = (double *)calloc(3*xdv*gratio*zdf,sizeof(double));

//	MPI_Barrier(MPI_COMM_WORLD);
//	printf("before the loop 1\n");
//	MPI_Barrier(MPI_COMM_WORLD);
	
	num=pd.myrank;
	
	for (a=0;a<pd.numproc;a++)
	{
		for (k=0;k<zdc;k++)
		{
			for (j=0;j<ydc;j++)
			{
				for (i=0;i<xdc;i++)
				{
					for (ito=0;ito<xdc;ito++)
					{
						iito = (num*xdc + ito) - (pd.myrank*xdc+i);
						if (iito<0)
							iito += xdv;
						 
						Fx(xcorc,iito,k,0,xdv,zdc,3) += Fb(V.velc,i,j,k,0,xdc,ydc,zdc,3)*Fb(V.velrecvc,ito,j,k,0,xdc,ydc,zdc,3);
						Fx(xcorc,iito,k,1,xdv,zdc,3) += Fb(V.velc,i,j,k,1,xdc,ydc,zdc,3)*Fb(V.velrecvc,ito,j,k,1,xdc,ydc,zdc,3);
						Fx(xcorc,iito,k,2,xdv,zdc,3) += Fb(V.velc,i,j,k,2,xdc,ydc,zdc,3)*Fb(V.velrecvc,ito,j,k,2,xdc,ydc,zdc,3);
					}
				}
			}
		}
		MPI_Sendrecv_replace(V.velrecvc, (3*xdc*ydc*zdc), MPI_DOUBLE,pd.right, stag, pd.left, rtag, MPI_COMM_WORLD,&status);
		MPI_Sendrecv_replace(&num, 1, MPI_INT,pd.right, stag, pd.left, rtag, MPI_COMM_WORLD,&status);
	}
//	MPI_Barrier(MPI_COMM_WORLD);
//        printf("before the loop 2\n");
//        MPI_Barrier(MPI_COMM_WORLD);
	num=pd.myrank;
	for (a=0;a<pd.numproc;a++)
	{
		for (k=0;k<zdf;k++)
		{
			for (j=0;j<ydf;j++)
			{
				for (i=0;i<xdf;i++)
				{
					for (ito=0;ito<xdf;ito++)
					{
						iito = (num*xdf + ito) - (pd.myrank*xdf+i);
						if (iito<0)
							iito += gratio*xdv;
						 
						Fx(xcorfl,iito,k,0,xdv*gratio,zdf,3) += Fb(V.velfl,i,j,k,0,xdf,ydf,zdf,3)*Fb(V.velrecvfl,ito,j,k,0,xdf,ydf,zdf,3);
						Fx(xcorfl,iito,k,1,xdv*gratio,zdf,3) += Fb(V.velfl,i,j,k,1,xdf,ydf,zdf,3)*Fb(V.velrecvfl,ito,j,k,1,xdf,ydf,zdf,3);
						Fx(xcorfl,iito,k,2,xdv*gratio,zdf,3) += Fb(V.velfl,i,j,k,2,xdf,ydf,zdf,3)*Fb(V.velrecvfl,ito,j,k,2,xdf,ydf,zdf,3);
					}
				}
			}
		}
		MPI_Sendrecv_replace(V.velrecvfl, (3*xdf*ydf*zdf), MPI_DOUBLE,pd.right, stag, pd.left, rtag, MPI_COMM_WORLD,&status);
		MPI_Sendrecv_replace(&num, 1, MPI_INT,pd.right, stag, pd.left, rtag, MPI_COMM_WORLD,&status);
	}
//	MPI_Barrier(MPI_COMM_WORLD);
//        printf("before the loop 1\n");
//        MPI_Barrier(MPI_COMM_WORLD);
	num=pd.myrank;
	for (a=0;a<pd.numproc;a++)
	{
		for (k=0;k<zdf;k++)
		{
			for (j=0;j<ydf;j++)
			{
				for (i=0;i<xdf;i++)
				{
					for (ito=0;ito<xdf;ito++)
					{
						iito = (num*xdf + ito) - (pd.myrank*xdf+i);
						if (iito<0)
							iito += gratio*xdv;
						 
						Fx(xcorfu,iito,k,0,xdv*gratio,zdf,3) += Fb(V.velfu,i,j,k,0,xdf,ydf,zdf,3)*Fb(V.velrecvfu,ito,j,k,0,xdf,ydf,zdf,3);
						Fx(xcorfu,iito,k,1,xdv*gratio,zdf,3) += Fb(V.velfu,i,j,k,1,xdf,ydf,zdf,3)*Fb(V.velrecvfu,ito,j,k,1,xdf,ydf,zdf,3);
						Fx(xcorfu,iito,k,2,xdv*gratio,zdf,3) += Fb(V.velfu,i,j,k,2,xdf,ydf,zdf,3)*Fb(V.velrecvfu,ito,j,k,2,xdf,ydf,zdf,3);
					}
				}
			}
		}
		MPI_Sendrecv_replace(V.velrecvfu, (3*xdf*ydf*zdf), MPI_DOUBLE,pd.right, stag, pd.left, rtag, MPI_COMM_WORLD,&status);
		MPI_Sendrecv_replace(&num, 1, MPI_INT,pd.right, stag, pd.left, rtag, MPI_COMM_WORLD,&status);
	}
//	MPI_Barrier(MPI_COMM_WORLD);
//        printf("loops ended\n");
//        MPI_Barrier(MPI_COMM_WORLD);
/*
	for (k=0;k<zd;k++)
	{
		for (i=0;i<xdv;i++)
		{
			Fx(xcor,i,k,0,xdv,zd,3) /= ((double)(yd*xdv));
			Fx(xcor,i,k,1,xdv,zd,3) /= ((double)(yd*xdv));
			Fx(xcor,i,k,2,xdv,zd,3) /= ((double)(yd*xdv));
		}
	}
*/			
	for (i=0;i<(xdv*zdc*3);i++)
		*(V.xcorlc+i) = 0.;
	for (i=0;i<(xdv*gratio*zdf*3);i++){
		*(V.xcorlfl+i) = 0.;
		*(V.xcorlfu+i) = 0.;
	}

	MPI_Reduce(xcorc,V.xcorlc,(3*xdv*zdc),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(xcorfl,V.xcorlfl,(3*xdv*gratio*zdf),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(xcorfu,V.xcorlfu,(3*xdv*gratio*zdf),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	
	free(xcorc);
	free(xcorfl);
	free(xcorfu);
return V;
}
