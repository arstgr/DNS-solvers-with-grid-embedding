/********************************************************************************
 * Amirreza Rastegari							                             	*
 * arstgri@gmail.com						                             		*
 * This code is used to calculate the 2 point velocity correlations in parallel *
 * Formulation:								                                 	*
 * Mathiu & Scott, An introduction to turbulent flows, Cambridge university	    *
 * press, 2000, pp. 62								                            *
 ********************************************************************************/
#include "definitions.h"

POINTER mean_subtract(int xdc,int ydc,int zdc,int xdf,int ydf,int zdf, POINTER V,PDATA pd)
{
	int i,j,k;
	for (i=0;i<xdc;i++)
	{
		for (j=0;j<ydc;j++)
		{
			for (k=0;k<zdc;k++)
			{
				Fb(V.velc,i,j,k,0,xdc,ydc,zdc,3) -= V.uavgc[k];
				Fb(V.velc,i,j,k,1,xdc,ydc,zdc,3) -= V.vavgc[k];
				Fb(V.velc,i,j,k,2,xdc,ydc,zdc,3) -= V.wavgc[k];
			}
		}
	}
	for (i=0;i<xdf;i++)
	{
		for (j=0;j<ydf;j++)
		{
			for (k=0;k<zdf;k++)
			{
				Fb(V.velfl,i,j,k,0,xdf,ydf,zdf,3) -= V.uavgfl[k];
				Fb(V.velfl,i,j,k,1,xdf,ydf,zdf,3) -= V.vavgfl[k];
				Fb(V.velfl,i,j,k,2,xdf,ydf,zdf,3) -= V.wavgfl[k];
				
				Fb(V.velfu,i,j,k,0,xdf,ydf,zdf,3) -= V.uavgfu[k];
				Fb(V.velfu,i,j,k,1,xdf,ydf,zdf,3) -= V.vavgfu[k];
				Fb(V.velfu,i,j,k,2,xdf,ydf,zdf,3) -= V.wavgfu[k];
			}
		}
	}
}