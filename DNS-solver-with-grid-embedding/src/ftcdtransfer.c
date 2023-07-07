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

/* Data transfer fom fine to coarse grid: first layer */
POINTER ftcdtransfer(int xdc, int ydc, int zdc, int xdf, int ydf, int zdf, int itr, double tauc, double tauf, int GR, POINTER V, PDATA pd)
{
	int i,j,k,a,q,d;
	double ux,uy,uz,dens;
	double cnt,edg,fce;
	double ueq[3],feq,fneq,v,u;
	int ic,jc,kc;
	double tmp;
	const int UP[5]={5,11,14,15,18};
	const int DW[5]={6,12,13,16,17};
	
// Fine to coarse, lower side
	for (ic=0;ic<xdc;ic++)
	{
		for (jc=0;jc<ydc;jc++)
		{
			i=ic*GR;
			j=jc*GR;
			k= zdf-1-GR;
			if (ic==(xdc-1))
			  i=xdf-1;
			if (jc==(ydc-1))
			  j=ydf-1;
			v = 0.;
			ux=0.;uy=0.;uz=0.;dens=0.;
			for (a=0;a<19;a++)
			{
				ux += E0[a]*Fb(V.slf,i,j,k,a,xdf,ydf,zdf,19);
				uy += E1[a]*Fb(V.slf,i,j,k,a,xdf,ydf,zdf,19);
				uz += E2[a]*Fb(V.slf,i,j,k,a,xdf,ydf,zdf,19);
				dens += Fb(V.slf,i,j,k,a,xdf,ydf,zdf,19);
			}
			ueq[0]=ux/dens;
			ueq[1]=uy/dens;
			ueq[2]=uz/dens;
			v = ueq[0]*ueq[0]+ueq[1]*ueq[1]+ueq[2]*ueq[2];

/*			for (a=0;a<5;a++)
			{
				u=E0[UP[a]]*ueq[0]+E1[UP[a]]*ueq[1]+E2[UP[a]]*ueq[2];
				feq=w[UP[a]]*dens*(1.+3.*u+4.5*u*u-1.5*v);
				fneq=Fb(V.slf,i,j,k,UP[a],xdf,ydf,zdf,19)-feq;
				Fb(V.f,ic,jc,0,UP[a],xdc,ydc,zdc,19)= feq + GR*(tauf/tauc)*fneq ;
				Fb(V.ftemp,ic,jc,0,UP[a],xdc,ydc,zdc,19)= feq + GR*(tauf/tauc)*fneq ;
			}*/					
			for (a=0;a<19;a++)
			{
				u=E0[a]*ueq[0]+E1[a]*ueq[1]+E2[a]*ueq[2];
				feq=w[a]*dens*(1.+3.*u+4.5*u*u-1.5*v);
				fneq=Fb(V.slf,i,j,k,a,xdf,ydf,zdf,19)-feq;
				Fb(V.f,ic,jc,0,a,xdc,ydc,zdc,19)= feq + GR*(tauf/tauc)*fneq ;
				Fb(V.ftemp,ic,jc,0,a,xdc,ydc,zdc,19)= feq + GR*(tauf/tauc)*fneq ;
			}
		}
	}

// Fine to coarse, upper side
	for (ic=0;ic<xdc;ic++)
	{
		for (jc=0;jc<ydc;jc++)
		{
			i=ic*GR;
			j=jc*GR;
			k= GR;
			if (ic==(xdc-1))
			  i=xdf-1;
			if (jc==(ydc-1))
			  j=ydf-1;
			v = 0.;
			ux=0.;uy=0.;uz=0.;dens=0.;
			for (a=0;a<19;a++)
			{
				ux += E0[a]*Fb(V.suf,i,j,k,a,xdf,ydf,zdf,19);
				uy += E1[a]*Fb(V.suf,i,j,k,a,xdf,ydf,zdf,19);
				uz += E2[a]*Fb(V.suf,i,j,k,a,xdf,ydf,zdf,19);
				dens += Fb(V.suf,i,j,k,a,xdf,ydf,zdf,19);
			}
			ueq[0]=ux/dens;
			ueq[1]=uy/dens;
			ueq[2]=uz/dens;
			v = ueq[0]*ueq[0]+ueq[1]*ueq[1]+ueq[2]*ueq[2];
			
// 			for (a=0;a<5;a++)
// 			{
// 				u=E0[DW[a]]*ueq[0]+E1[DW[a]]*ueq[1]+E2[DW[a]]*ueq[2];
// 				feq=w[DW[a]]*dens*(1.+3.*u+4.5*u*u-1.5*v);
// 				fneq=Fb(V.suf,i,j,k,DW[a],xdf,ydf,zdf,19)-feq;
// 				Fb(V.f,ic,jc,(zdc-1),DW[a],xdc,ydc,zdc,19) = feq + GR*(tauf/tauc)*fneq;
// 				Fb(V.ftemp,ic,jc,(zdc-1),DW[a],xdc,ydc,zdc,19) = feq + GR*(tauf/tauc)*fneq;
// 			}
			for (a=0;a<19;a++)
			{
				u=E0[a]*ueq[0]+E1[a]*ueq[1]+E2[a]*ueq[2];
				feq=w[a]*dens*(1.+3.*u+4.5*u*u-1.5*v);
				fneq=Fb(V.suf,i,j,k,a,xdf,ydf,zdf,19)-feq;
				Fb(V.f,ic,jc,(zdc-1),a,xdc,ydc,zdc,19) = feq + GR*(tauf/tauc)*fneq;
				Fb(V.ftemp,ic,jc,(zdc-1),a,xdc,ydc,zdc,19) = feq + GR*(tauf/tauc)*fneq;
			}
		}
	}
//	fprintf(stderr,"in function ftcdtransfer, c error handler is %s\n",strerror(errno));
	return V;
}