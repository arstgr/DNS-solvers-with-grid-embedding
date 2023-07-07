/***************************************************************************
 *   Amirreza Rastegari                                                    *
 *   arstgri@gmail.com                                                     *
 *                                                                         *
 *   Parallel LBM-BGK code with constant flow rate Q for DNS of turbulent  *
 *   channel flows                                                         *
 *   This program uses a new formulation for the constant flow rate        *
 *   implementation                                                        *
 *   A multigrid implementations is used for increasing the accuracy       *
 *   in the near wall reagion                                              *
 ***************************************************************************/

#include "definitions.h"

POINTER init(int xdc, int ydc, int zdc, int xdf, int ydf, int zdf, POINTER V)
{ 
	int i,j,k,mm;
	for (i=0; i<xdc; i++)
		for (j=0; j<ydc; j++)
			for (k=0; k<zdc; k++)
				if (Fs(V.s,i,j,k,xdc,ydc,zdc)>-1)
				{
					for (mm=0;mm<19;mm++)
					{
						Fb(V.f,i,j,k,mm,xdc,ydc,zdc,19)=w[mm];
						Fb(V.ftemp,i,j,k,mm,xdc,ydc,zdc,19)=w[mm];
					}
				}

	for (i=0; i<xdf; i++)
		for (j=0; j<ydf; j++)
			for (k=0; k<zdf; k++)
				if (Fs(V.sus,i,j,k,xdf,ydf,zdf)>-1)
				{
					for (mm=0;mm<19;mm++)
					{
						Fb(V.suf,i,j,k,mm,xdf,ydf,zdf,19)=w[mm];
						Fb(V.suftemp,i,j,k,mm,xdf,ydf,zdf,19)=w[mm];
					}
				}

	for (i=0; i<xdf; i++)
		for (j=0; j<ydf; j++)
			for (k=0; k<zdf; k++)
				if (Fs(V.sls,i,j,k,xdf,ydf,zdf)>-1)
				{
					for (mm=0;mm<19;mm++)
					{
						Fb(V.slf,i,j,k,mm,xdf,ydf,zdf,19)=w[mm];
						Fb(V.slftemp,i,j,k,mm,xdf,ydf,zdf,19)=w[mm];
					}
				}
				
//	fprintf(stderr,"in function init, c error handler is %s\n",strerror(errno));
return V;
}