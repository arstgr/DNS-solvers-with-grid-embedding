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

POINTER output(int xdc, int ydc, int zdc, int xdf, int ydf, int zdf, int GR, POINTER V, DP Fx)
{
	int i, j, k,a,q;
	FILE *ooutput;
	DP ux,uy,uz,dens;
	double *zps;
	int NZ,NZs;
	double h;
	
	NZs=2*zdf+GR*(zdc-4)+GR-1;
	NZ=2*zdf+(zdc-4);
	
	h=(NZs-2.)/2.;
	
	zps=(DP *)calloc(NZ,sizeof(DP));
	
	for (k=1;k<zdf;k++)
	  *(zps+k) = k-0.5;
	for (k=2;k<(zdc-2);k++)
	  *(zps+zdf-2+k) = GR*(k-1)+(zdf-1-0.5);
	for (k=0;k<(zdf-1);k++)
	  *(zps+zdf+zdc-4+k) = k + GR*(zdc-3) + zdf-1.-0.5;
	*(zps+2*zdf+zdc-5) = 2*zdf+GR*(zdc-4)+GR-3;

	ooutput=fopen("results.dat","w");
	fprintf(ooutput,"TITLE = \"3d velocity field\"\n");
	fprintf(ooutput,"VARIABLES = \"x/h\", \"y/h\", \"z/h\", \"rho\", \"u\", \"v\", \"w\" \n");
	fprintf(ooutput,"ZONE I=%d, J=%d, K=%d, F=POINT",xdf-1,ydf,zdf);
	for (k=0;k<zdf;k++)
	{
		for (j=0;j<ydf;j++)
		{
			for (i=0;i<(xdf-1);i++)
			{
				ux=0.;uy=0.;uz=0.;dens=0.;
				for (a=0;a<19;a++)
				{
					ux += E0[a]*Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
					uy += E1[a]*Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
					uz += E2[a]*Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
					dens += Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
				}
				fprintf(ooutput,"%f %f %f %f %f %f %f\n",(i+0.5)/h,(j+0.5)/h,(k-0.5)/h,dens,ux/ub,uy/ub,uz/ub);
			}
		}
	}

	fprintf(ooutput,"ZONE I=%d, J=%d, K=%d, F=POINT",xdc-1,ydc,zdc-4);
	for (k=2;k<(zdc-2);k++)
	{
		for (j=0;j<ydc;j++)
		{
			for (i=0;i<(xdc-1);i++)
			{
				ux=0.;uy=0.;uz=0.;dens=0.;
				for (a=0;a<19;a++)
				{
					ux += E0[a]*Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);
					uy += E1[a]*Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);
					uz += E2[a]*Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);
					dens += Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);
				}
				fprintf(ooutput,"%f %f %f %f %f %f %f\n",(i*GR+0.5)/h,(j*GR+0.5)/h,*(zps+zdf-2+k)/h,dens,ux/ub,uy/ub,uz/ub);
			}
		}
	}
	
	fprintf(ooutput,"ZONE I=%d, J=%d, K=%d, F=POINT",xdf-1,ydf,zdf);
	for (k=0;k<zdf;k++)
	{
		for (j=0;j<ydf;j++)
		{
			for (i=0;i<(xdf-1);i++)
			{
				ux=0.;uy=0.;uz=0.;dens=0.;
				for (a=0;a<19;a++)
				{
					ux += E0[a]*Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
					uy += E1[a]*Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
					uz += E2[a]*Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
					dens += Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
				}
				fprintf(ooutput,"%f %f %f %f %f %f %f\n",(i+0.5)/h,(j+0.5)/h,*(zps+zdf+zdc-4+k)/h,dens,ux/ub,uy/ub,uz/ub);
			}
		}
	}
	
	fprintf(ooutput,"\n\n");
	fclose(ooutput);
	
	free(zps);
return V;
}