/********************************************************************************
 * Amirreza Rastegari							                             	*
 * arstgri@gmail.com						                             		*
 * This code is used to calculate the 2 point velocity correlations in parallel *
 * Formulation:								                                 	*
 * Mathiu & Scott, An introduction to turbulent flows, Cambridge university	    *
 * press, 2000, pp. 62								                            *
 ********************************************************************************/
#include "definitions.h"
#include "fun_defs.h"

POINTER mean(int xdc,int ydc,int zdc,int xdf,int ydf,int zdf, POINTER V,PDATA pd)
{
	int i,j,k,a,q,b,c;
	int cntr;
	int num=0,temp=0,upn=0;
	DP uavgc[zdc],vavgc[zdc],wavgc[zdc];
	DP uavgfl[zdf],vavgfl[zdf],wavgfl[zdf];
	DP uavgfu[zdf],vavgfu[zdf],wavgfu[zdf];
	double temp1,temp2;
	DP *DUMMYYZc,*DUMMYYZf;
	MPI_Status status;
	
//	fprintf(stderr,"starting mean subtract i'm %d\n",pd.myrank);
	
	DUMMYYZc = (double *)calloc(zdc,sizeof(double));
	DUMMYYZf = (double *)calloc(zdf,sizeof(double));
	
	for (cntr=TSTR;cntr<=TEND;cntr+=DT)
	{
		V=reading(xdc, ydc, zdc, xdf,ydf,zdf, V, pd, cntr);
//		MPI_Barrier(MPI_COMM_WORLD);
//		printf("read %d\n",cntr);
//		MPI_Barrier(MPI_COMM_WORLD);
		
		for (i=0;i<xdc;i++)
		{
			for (j=0;j<ydc;j++)
			{
				for (k=0;k<zdc;k++)
				{
					V.uavgc[k] += Fb(V.velc,i,j,k,0,xdc,ydc,zdc,3);
					V.vavgc[k] += Fb(V.velc,i,j,k,1,xdc,ydc,zdc,3);
					V.wavgc[k] += Fb(V.velc,i,j,k,2,xdc,ydc,zdc,3);
				}
			}
		}
		for (i=0;i<xdf;i++)
		{
			for (j=0;j<ydf;j++)
			{
				for (k=0;k<zdf;k++)
				{
					V.uavgfl[k] += Fb(V.velfl,i,j,k,0,xdf,ydf,zdf,3);
					V.vavgfl[k] += Fb(V.velfl,i,j,k,1,xdf,ydf,zdf,3);
					V.wavgfl[k] += Fb(V.velfl,i,j,k,2,xdf,ydf,zdf,3);
					
					V.uavgfu[k] += Fb(V.velfu,i,j,k,0,xdf,ydf,zdf,3);
					V.vavgfu[k] += Fb(V.velfu,i,j,k,1,xdf,ydf,zdf,3);
					V.wavgfu[k] += Fb(V.velfu,i,j,k,2,xdf,ydf,zdf,3);
				}
			}
		}
	}
	num = (TEND-TSTR+1)/DT;
	for (k=0;k<zdc;k++)
	{
		V.uavgc[k] /= (xdc*ydc*pd.numproc*num);
		V.vavgc[k] /= (xdc*ydc*pd.numproc*num);
		V.wavgc[k] /= (xdc*ydc*pd.numproc*num);
	}
	for (k=0;k<zdf;k++)
	{
		V.uavgfl[k] /= (xdf*ydf*pd.numproc*num);
		V.vavgfl[k] /= (xdf*ydf*pd.numproc*num);
		V.wavgfl[k] /= (xdf*ydf*pd.numproc*num);
		
		V.uavgfu[k] /= (xdf*ydf*pd.numproc*num);
		V.vavgfu[k] /= (xdf*ydf*pd.numproc*num);
		V.wavgfu[k] /= (xdf*ydf*pd.numproc*num);
	}
//	MPI_Barrier(MPI_COMM_WORLD);
//	printf("start of communication\n");
//	MPI_Barrier(MPI_COMM_WORLD);	
	MPI_Allreduce(V.uavgc,DUMMYYZc,zdc,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	memcpy(V.uavgc,DUMMYYZc,zdc*sizeof(double));
	
	MPI_Allreduce(V.vavgc,DUMMYYZc,zdc,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	memcpy(V.vavgc,DUMMYYZc,zdc*sizeof(double));
	
	MPI_Allreduce(V.wavgc,DUMMYYZc,zdc,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	memcpy(V.wavgc,DUMMYYZc,zdc*sizeof(double));
	
	MPI_Allreduce(V.uavgfl,DUMMYYZf,zdf,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	memcpy(V.uavgfl,DUMMYYZf,zdf*sizeof(double));
	
	MPI_Allreduce(V.vavgfl,DUMMYYZf,zdf,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	memcpy(V.vavgfl,DUMMYYZf,zdf*sizeof(double));
	
	MPI_Allreduce(V.wavgfl,DUMMYYZf,zdf,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	memcpy(V.wavgfl,DUMMYYZf,zdf*sizeof(double));
	
	MPI_Allreduce(V.uavgfu,DUMMYYZf,zdf,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	memcpy(V.uavgfu,DUMMYYZf,zdf*sizeof(double));
	
	MPI_Allreduce(V.vavgfu,DUMMYYZf,zdf,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	memcpy(V.vavgfu,DUMMYYZf,zdf*sizeof(double));
	
	MPI_Allreduce(V.wavgfu,DUMMYYZf,zdf,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	memcpy(V.wavgfu,DUMMYYZf,zdf*sizeof(double));
	
//	MPI_Barrier(MPI_COMM_WORLD);
//        printf("end of communication\n");
//        MPI_Barrier(MPI_COMM_WORLD);

	for (k=0;k<(zdc+1)/2;k++)
	{
		temp1 = V.uavgc[k];
		temp2 = V.uavgc[zdc-k-1];
		
		V.uavgc[k] = 0.5*(temp1+temp2);
		V.uavgc[zdc-k-1] = 0.5*(temp1+temp2);
		
		temp1 = V.vavgc[k];
		temp2 = V.vavgc[zdc-k-1];
		
		V.vavgc[k] = 0.5*(temp1+temp2);
		V.vavgc[zdc-k-1] = 0.5*(temp1+temp2);
		
		temp1 = V.wavgc[k];
		temp2 = V.wavgc[zdc-k-1];
		
		V.wavgc[k] = 0.5*(temp1-temp2);
		V.wavgc[zdc-k-1] = -0.5*(temp1-temp2);
	}
//	MPI_Barrier(MPI_COMM_WORLD);
//        printf("end of first loop\n");
//        MPI_Barrier(MPI_COMM_WORLD);
	for (k=0;k<zdf;k++)
	{
		temp1 = V.uavgfl[k];
		temp2 = V.uavgfu[zdf-k-1];
		
		V.uavgfl[k] = 0.5*(temp1+temp2);
		V.uavgfu[zdf-k-1] = 0.5*(temp1+temp2);
		
		
		temp1 = V.vavgfl[k];
		temp2 = V.vavgfu[zdf-k-1];
		
		V.vavgfl[k] = 0.5*(temp1+temp2);
		V.vavgfu[zdf-k-1] = 0.5*(temp1+temp2);
		
		
		temp1 = V.wavgfl[k];
		temp2 = V.wavgfu[zdf-k-1];
		
		V.wavgfl[k] = 0.5*(temp1-temp2);
		V.wavgfu[zdf-k-1] = -0.5*(temp1-temp2);
	}
	
	free(DUMMYYZc);
	free(DUMMYYZf);
//	MPI_Barrier(MPI_COMM_WORLD);
//	fprintf(stderr,"end of mean subtract i'm %d\n",pd.myrank);
//  	MPI_Barrier(MPI_COMM_WORLD);
  return V;
}