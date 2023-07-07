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

POINTER fine_force(int xdf, int ydf , int zdf, int zdc, int GR, POINTER V, PDATA pd, DP ubulk)
{
	int i,j,k,a;
	double coef;
	int NZs;
	double momentum=0.,tmomentum=0.,tdens=0.,denst=0.;
	double ux,uy,uz,dens;
	NZs=2*zdf+GR*(zdc-4)+GR-1;
	
	for (i=0;i<(xdf-1);i++)
	{
		for (j=0;j<(ydf-1);j++)
		{
			for (k=1;k<zdf;k++)
			{
				ux=0.;uy=0.;uz=0.;dens=0.;
				for (a=0;a<19;a++)
				{
					dens+=Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
					ux+= E0[a]*Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
					uy+= E1[a]*Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
					uz+= E2[a]*Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
				}
				
				if (k==(zdf-1))
				  coef=0.5;
				else
				  coef=1.;
				momentum += coef*ux;
				denst += coef*dens;
			}
		}
	}
	
	for (i=0;i<(xdf-1);i++)
	{
		for (j=0;j<(ydf-1);j++)
		{
			for (k=0;k<(zdf-1);k++)
			{
				ux=0.;uy=0.;uz=0.;dens=0.;
				for (a=0;a<19;a++)
				{
					dens+=Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
					ux+= E0[a]*Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
					uy+= E1[a]*Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
					uz+= E2[a]*Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
				}
				if (k==0)
				  coef=0.5;
				else
				  coef=1.;
				momentum += coef*ux;
				denst += coef*dens;
			}
		}
	}
	
	momentum /= ((double)((xdf-1.)*(ydf-1.)*XDIM*YDIM*(NZs-2.)));
	denst /= ((double)((xdf-1.)*(ydf-1.)*XDIM*YDIM*(NZs-2.)));
	MPI_Allreduce(&momentum,&tmomentum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	momentum=tmomentum;

	MPI_Allreduce(&denst,&tdens,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	denst=tdens;
	
	V.force[2]=momentum;
	V.rho[2]=denst;
	
	V.force[0] = (V.rho[1]+V.rho[2])*ubulk - (V.force[1]+V.force[2]);
  
	return V;
}