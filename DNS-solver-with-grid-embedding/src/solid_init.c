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

/* Initializing the solid surfaces 
   solid no-slip     s = -1
   solid slip        s = -2
   fluid             s = 0
   periodic boundary s = 1
*/

POINTER solid_init(int xdc, int ydc, int zdc, int xdf, int ydf, int zdf, POINTER V)
{
	int i,j,k;
	int solid=-1,periodic=1,fluid=0;
	
	for (i=0;i<xdc;i++)
		for (j=0;j<ydc;j++)
			for (k=0;k<zdc;k++)
				Fs(V.s,i,j,k,xdc,ydc,zdc)=fluid;

	for (i=0;i<xdf;i++)
		for (j=0;j<ydf;j++)
			for (k=0;k<zdf;k++)
				Fs(V.sus,i,j,k,xdf,ydf,zdf)=fluid;

	for (i=0;i<xdf;i++)
		for (j=0;j<ydf;j++)
			for (k=0;k<zdf;k++)
				Fs(V.sls,i,j,k,xdf,ydf,zdf)=fluid;
	
	k=0;
	for (i=0;i<xdc;i++)
		for (j=0;j<ydc;j++)
			Fs(V.s,i,j,k,xdc,ydc,zdc)=periodic;

	k=0;
	for (i=0;i<xdf;i++)
		for (j=0;j<ydf;j++)
			Fs(V.sls,i,j,k,xdf,ydf,zdf)=solid;

	k=0;
	for (i=0;i<xdf;i++)
		for (j=0;j<ydf;j++)
			Fs(V.sus,i,j,k,xdf,ydf,zdf)=periodic;

	k=zdc-1;
	for (i=0;i<xdc;i++)
		for (j=0;j<ydc;j++)
			Fs(V.s,i,j,k,xdc,ydc,zdc)=periodic;

	k=zdf-1;
	for (i=0;i<xdf;i++)
		for (j=0;j<ydf;j++)
			Fs(V.sus,i,j,k,xdf,ydf,zdf)=solid;

	k=zdf-1;
	for (i=0;i<xdf;i++)
		for (j=0;j<ydf;j++)
			Fs(V.sls,i,j,k,xdf,ydf,zdf)=periodic;
	
	i=0;
	for (j=0;j<ydc;j++)
		for (k=0;k<zdc;k++)
			Fs(V.s,i,j,k,xdc,ydc,zdc)=periodic;

	i=0;
	for (j=0;j<ydf;j++)
		for (k=0;k<(zdf-1);k++)
			Fs(V.sus,i,j,k,xdf,ydf,zdf)=periodic;

	i=0;
	for (j=0;j<ydf;j++)
		for (k=1;k<zdf;k++)
			Fs(V.sls,i,j,k,xdf,ydf,zdf)=periodic;

	i=xdc-1;
	for (j=0;j<ydc;j++)
		for (k=0;k<zdc;k++)
			Fs(V.s,i,j,k,xdc,ydc,zdc)=periodic;

	i=xdf-1;
	for (j=0;j<ydf;j++)
		for (k=0;k<(zdf-1);k++)
			Fs(V.sus,i,j,k,xdf,ydf,zdf)=periodic;

	i=xdf-1;
	for (j=0;j<ydf;j++)
		for (k=1;k<zdf;k++)
			Fs(V.sls,i,j,k,xdf,ydf,zdf)=periodic;
	
	j=0;
	for (i=0;i<xdc;i++)
		for (k=0;k<zdc;k++)
			Fs(V.s,i,j,k,xdc,ydc,zdc)=periodic;

	j=0;
	for (i=0;i<xdf;i++)
		for (k=0;k<(zdf-1);k++)
			Fs(V.sus,i,j,k,xdf,ydf,zdf)=periodic;

	j=0;
	for (i=0;i<xdf;i++)
		for (k=1;k<zdf;k++)
			Fs(V.sls,i,j,k,xdf,ydf,zdf)=periodic;

	j=ydc-1;
	for (i=0;i<xdc;i++)
		for (k=0;k<zdc;k++)
			Fs(V.s,i,j,k,xdc,ydc,zdc)=periodic;

	j=ydf-1;
	for (i=0;i<xdf;i++)
		for (k=0;k<(zdf-1);k++)
			Fs(V.sus,i,j,k,xdf,ydf,zdf)=periodic;

	j=ydf-1;
	for (i=0;i<xdf;i++)
		for (k=1;k<zdf;k++)
			Fs(V.sls,i,j,k,xdf,ydf,zdf)=periodic;
	return V;
	
//	fprintf(stderr,"in function solid init, c error handler is %s\n",strerror(errno));
}