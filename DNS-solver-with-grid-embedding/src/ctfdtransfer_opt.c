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
POINTER ctfdtransfer_opt(int xdc, int ydc, int zdc, int xdf, int ydf, int zdf, int itr, double tauc, double tauf, int GR, POINTER V, PDATA pd)
{
	int i,j,k,a,q;
	int ii, jj,nn;
	double ux,uy,uz,dens;
	double ueq[3],feq,fneq,v,u;
	int ic,jc,kc;
	double tt,uu;
	double c[23][4][4]={0.};
	double ans=0.;
	int ict,jct;
	double tmp;
	
// Coarse to fine, lower side to zdf-1
	for (ic=0;ic<xdc;ic++)
	{
		for (jc=0;jc<ydc;jc++)
		{
			bcucof_opt(xdc, ydc, ic, jc, GR, V.lint, &c[0][0][0]);
			
			kc=1;
			if (ic<(xdc-1))
				ict = GR;
			else
				ict = 1;
			
			if (jc<(ydc-1))
				jct = GR;
			else
				jct = 1;
			
			for (ii=0;ii<ict;ii++)
			{
				for (jj=0;jj<jct;jj++)
				{
					i = ic * GR + ii;
					j = jc * GR + jj;
					
					tmp = i%GR; tt = ((double)tmp)/((double)GR);
					tmp = j%GR; uu = ((double)tmp)/((double)GR);
					
					dens = 0.;
					for (nn=3;nn>-1;nn--)
					{
					  dens=tt*dens+((c[0][nn][3]*uu+c[0][nn][2])*uu+c[0][nn][1])*uu+c[0][nn][0];
					}
					
					ueq[0] = 0.;
					for (nn=3;nn>-1;nn--)
					{
					  ueq[0]=tt*ueq[0]+((c[1][nn][3]*uu+c[1][nn][2])*uu+c[1][nn][1])*uu+c[1][nn][0];
					}
					
					ueq[1] = 0.;
					for (nn=3;nn>-1;nn--)
					{
					  ueq[1]=tt*ueq[1]+((c[2][nn][3]*uu+c[2][nn][2])*uu+c[2][nn][1])*uu+c[2][nn][0];
					}
					
					ueq[2] = 0.;
					for (nn=3;nn>-1;nn--)
					{
					  ueq[2]=tt*ueq[2]+((c[3][nn][3]*uu+c[3][nn][2])*uu+c[3][nn][1])*uu+c[3][nn][0];
					}
					v = ueq[0]*ueq[0]+ueq[1]*ueq[1]+ueq[2]*ueq[2];
					
					for (a=0;a<19;a++)
					{
						u=E0[a]*ueq[0]+E1[a]*ueq[1]+E2[a]*ueq[2];
						feq=w[a]*dens*(1.+3.*u+4.5*u*u-1.5*v);
						fneq=0.;
						for (nn=3;nn>-1;nn--)
						{
						  fneq=tt*fneq+((c[a+4][nn][3]*uu+c[a+4][nn][2])*uu+c[a+4][nn][1])*uu+c[a+4][nn][0];
						}
						
						Fb(V.slf,i,j,(zdf-1),a,xdf,ydf,zdf,19)=feq + (tauc/((double)GR*tauf))*fneq;
						Fb(V.slftemp,i,j,(zdf-1),a,xdf,ydf,zdf,19)=feq + (tauc/((double)GR*tauf))*fneq;
					}
				}
			}
		}
	}

// Coarse to fine, upper side to 0
	for (ic=0;ic<xdc;ic++)
	{
		for (jc=0;jc<ydc;jc++)
		{
			kc=zdc-2;
			
			bcucof_opt(xdc, ydc, ic, jc, GR, V.uint, &c[0][0][0]);
			
			if (ic<(xdc-1))
				ict = GR;
			else
				ict = 1;
			
			if (jc<(ydc-1))
				jct = GR;
			else
				jct = 1;
			
			for (ii=0;ii<ict;ii++)
			{
				for (jj=0;jj<jct;jj++)
				{
					i = ic * GR + ii;
					j = jc * GR + jj;
					
					tmp = i%GR; tt = ((double)tmp)/((double)GR);
					tmp = j%GR; uu = ((double)tmp)/((double)GR);
					
					dens = 0.;
					for (nn=3;nn>-1;nn--)
					{
					  dens=tt*dens+((c[0][nn][3]*uu+c[0][nn][2])*uu+c[0][nn][1])*uu+c[0][nn][0];
					}
					
					ueq[0] = 0.;
					for (nn=3;nn>-1;nn--)
					{
					  ueq[0]=tt*ueq[0]+((c[1][nn][3]*uu+c[1][nn][2])*uu+c[1][nn][1])*uu+c[1][nn][0];
					}
					
					ueq[1] = 0.;
					for (nn=3;nn>-1;nn--)
					{
					  ueq[1]=tt*ueq[1]+((c[2][nn][3]*uu+c[2][nn][2])*uu+c[2][nn][1])*uu+c[2][nn][0];
					}
					
					ueq[2] = 0.;
					for (nn=3;nn>-1;nn--)
					{
					  ueq[2]=tt*ueq[2]+((c[3][nn][3]*uu+c[3][nn][2])*uu+c[3][nn][1])*uu+c[3][nn][0];
					}
					v = ueq[0]*ueq[0]+ueq[1]*ueq[1]+ueq[2]*ueq[2];
					
					for (a=0;a<19;a++)
					{
						u=E0[a]*ueq[0]+E1[a]*ueq[1]+E2[a]*ueq[2];
						feq=w[a]*dens*(1.+3.*u+4.5*u*u-1.5*v);
						fneq = 0.;
						for (nn=3;nn>-1;nn--)
						{
						  fneq=tt*fneq+((c[a+4][nn][3]*uu+c[a+4][nn][2])*uu+c[a+4][nn][1])*uu+c[a+4][nn][0];
						}
						
						Fb(V.suf,i,j,0,a,xdf,ydf,zdf,19)= feq + (tauc/((double)GR*tauf))*fneq;
						Fb(V.suftemp,i,j,0,a,xdf,ydf,zdf,19)= feq + (tauc/((double)GR*tauf))*fneq;
					}
				}
			}
		}
	}
	
//	fprintf(stderr,"in function ctfdtransfer, c error handler is %s\n",strerror(errno));
	return V;
}