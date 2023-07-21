/********************************************************************************
 * Amirreza Rastegari							                             	*
 * arstgri@gmail.com						                             		*
 * This code is used to calculate the 2 point velocity correlations in parallel *
 * Formulation:								                                 	*
 * Mathiu & Scott, An introduction to turbulent flows, Cambridge university	    *
 * press, 2000, pp. 62								                            *
 ********************************************************************************/
#include "definitions.h"

POINTER ycorrel(int xdc,int ydc,int zdc, int xdf,int ydf,int zdf, POINTER V,PDATA pd)
{
	int i,j,k,a,b,c;
	double *ycorc,*ycorfl,*ycorfu;
	int jto,jjto;

	ycorc = (double *)calloc(3*ydc*zdc,sizeof(double));
	ycorfl = (double *)calloc(3*ydf*zdf,sizeof(double));
	ycorfu = (double *)calloc(3*ydf*zdf,sizeof(double));
	
	for (k=0;k<zdc;k++)
	{
		for (j=0;j<ydc;j++)
		{
			for (jto=0;jto<ydc;jto++)
			{
				jjto = jto -j;
				if (jjto < 0)
					jjto += ydc;

				for (i=0;i<xdc;i++)
				{
					Fy(ycorc,jjto,k,0,ydc,zdc,3) += Fb(V.velc,i,j,k,0,xdc,ydc,zdc,3)*Fb(V.velc,i,jto,k,0,xdc,ydc,zdc,3);
					Fy(ycorc,jjto,k,1,ydc,zdc,3) += Fb(V.velc,i,j,k,1,xdc,ydc,zdc,3)*Fb(V.velc,i,jto,k,1,xdc,ydc,zdc,3);
					Fy(ycorc,jjto,k,2,ydc,zdc,3) += Fb(V.velc,i,j,k,2,xdc,ydc,zdc,3)*Fb(V.velc,i,jto,k,2,xdc,ydc,zdc,3);
				}
			}
		}
	}
	for (k=0;k<zdf;k++)
	{
		for (j=0;j<ydf;j++)
		{
			for (jto=0;jto<ydf;jto++)
			{
				jjto = jto -j;
				if (jjto < 0)
					jjto += ydf;

				for (i=0;i<xdf;i++)
				{
					Fy(ycorfl,jjto,k,0,ydf,zdf,3) += Fb(V.velfl,i,j,k,0,xdf,ydf,zdf,3)*Fb(V.velfl,i,jto,k,0,xdf,ydf,zdf,3);
					Fy(ycorfl,jjto,k,1,ydf,zdf,3) += Fb(V.velfl,i,j,k,1,xdf,ydf,zdf,3)*Fb(V.velfl,i,jto,k,1,xdf,ydf,zdf,3);
					Fy(ycorfl,jjto,k,2,ydf,zdf,3) += Fb(V.velfl,i,j,k,2,xdf,ydf,zdf,3)*Fb(V.velfl,i,jto,k,2,xdf,ydf,zdf,3);
					
					Fy(ycorfu,jjto,k,0,ydf,zdf,3) += Fb(V.velfu,i,j,k,0,xdf,ydf,zdf,3)*Fb(V.velfu,i,jto,k,0,xdf,ydf,zdf,3);
					Fy(ycorfu,jjto,k,1,ydf,zdf,3) += Fb(V.velfu,i,j,k,1,xdf,ydf,zdf,3)*Fb(V.velfu,i,jto,k,1,xdf,ydf,zdf,3);
					Fy(ycorfu,jjto,k,2,ydf,zdf,3) += Fb(V.velfu,i,j,k,2,xdf,ydf,zdf,3)*Fb(V.velfu,i,jto,k,2,xdf,ydf,zdf,3);
				}
			}
		}
	}
/*
	for (k=0;k<zd;k++)
	{
		for (j=0;j<yd;j++)
		{
			Fy(ycor,j,k,0,yd,zd,3) /= ((double)(yd*xdv));
			Fy(ycor,j,k,1,yd,zd,3) /= ((double)(yd*xdv));
			Fy(ycor,j,k,2,yd,zd,3) /= ((double)(yd*xdv));
		}
	}
*/	
	for (j=0;j<(ydc*zdc*3);j++)
		V.ycorlc[j] = 0.;
	
	for (j=0;j<(ydf*zdf*3);j++){
		V.ycorlfl[j] = 0.;
		V.ycorlfu[j] = 0.;
	}
	
			
	MPI_Reduce(ycorc,V.ycorlc,(3*ydc*zdc),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(ycorfl,V.ycorlfl,(3*ydf*zdf),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(ycorfu,V.ycorlfu,(3*ydf*zdf),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
			
	free(ycorc);
	free(ycorfl);
	free(ycorfu);
return V;
}