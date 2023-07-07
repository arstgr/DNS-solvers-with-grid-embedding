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

/* Data Exchange between the coarse and fine grids */
/* From coarse to fine grids: first layer */
POINTER ctfdtransfer(int xdc, int ydc, int zdc, int xdf, int ydf, int zdf, int itr, double tauc, double tauf, int GR, POINTER V, PDATA pd)
{
	int i,j,k,a,q;
	double ux,uy,uz,dens;
	double ueq[3],feq,fneq,v,u;
	int ic,jc,kc;
	double tmp;
	
// Coarse to fine, lower side to zdf-1
	for (i=0;i<xdf;i++)
	{
		for (j=0;j<ydf;j++)
		{
			ic=i/GR;
			jc=j/GR;
			kc=1;
			if (i==(xdf-1))
				ic=xdc-1;
			if (j==(ydf-1))
				jc=ydc-1;

			dens = bcuint( xdc, ydc, GR, i, j, V.lint, 0);
			ueq[0]= bcuint( xdc, ydc, GR, i, j, V.lint, 1);
			ueq[1]= bcuint( xdc, ydc, GR, i, j, V.lint, 2);
			ueq[2]= bcuint( xdc, ydc, GR, i, j, V.lint, 3);
			v = ueq[0]*ueq[0]+ueq[1]*ueq[1]+ueq[2]*ueq[2];
			
			for (a=0;a<19;a++)
			{
				u=E0[a]*ueq[0]+E1[a]*ueq[1]+E2[a]*ueq[2];
				feq=w[a]*dens*(1.+3.*u+4.5*u*u-1.5*v);
				fneq = bcuint(xdc,ydc,GR,i,j,V.lint,(a+4));
				Fb(V.slf,i,j,(zdf-1),a,xdf,ydf,zdf,19)=feq + (tauc/((double)GR*tauf))*fneq;
				Fb(V.slftemp,i,j,(zdf-1),a,xdf,ydf,zdf,19)=feq + (tauc/((double)GR*tauf))*fneq;
			}
		}
	}

// Coarse to fine, upper side to 0
	for (i=0;i<xdf;i++)
	{
		for (j=0;j<ydf;j++)
		{
			ic=i/GR;
			jc=j/GR;
			kc=zdc-2;
			if (i==(xdf-1))
			  ic=xdc-1;
			if (j==(ydf-1))
			  jc=ydc-1;
			
			dens = bcuint( xdc, ydc, GR, i, j, V.uint, 0);
			ueq[0]= bcuint( xdc, ydc, GR, i, j, V.uint, 1);
			ueq[1]= bcuint( xdc, ydc, GR, i, j, V.uint, 2);
			ueq[2]= bcuint( xdc, ydc, GR, i, j, V.uint, 3);
			v = ueq[0]*ueq[0]+ueq[1]*ueq[1]+ueq[2]*ueq[2];

			for (a=0;a<19;a++)
			{
				u=E0[a]*ueq[0]+E1[a]*ueq[1]+E2[a]*ueq[2];
				feq=w[a]*dens*(1.+3.*u+4.5*u*u-1.5*v);
				fneq = bcuint(xdc,ydc,GR,i,j,V.uint,(a+4));
				Fb(V.suf,i,j,0,a,xdf,ydf,zdf,19)= feq + (tauc/((double)GR*tauf))*fneq;
				Fb(V.suftemp,i,j,0,a,xdf,ydf,zdf,19)= feq + (tauc/((double)GR*tauf))*fneq;
			}
		}
	}
	
//	fprintf(stderr,"in function ctfdtransfer, c error handler is %s\n",strerror(errno));
	return V;
}