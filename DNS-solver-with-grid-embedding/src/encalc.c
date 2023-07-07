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

/* Force calculaion for the constant flow rate formulation */
/* data is taken from all of the fine grids and non-overlapping parts of the coarse grid */
POINTER encalc(int xdc, int ydc , int zdc, int xdf, int ydf , int zdf, int GR, POINTER V, PDATA pd, int TM, long cntr,DP ubulk)
{
	int i,j,k,a;
	double coef,umax,dt;
	FILE *fh;
	double dens,ux,uy,uz;
	double momentum=0.,tmomentum=0.,tdens=0.,denst=0.;
	double cmomentum=0.,cdenst=0.,tcmomentum=0.,tcdenst=0.;
	char fn[30];
	int NZ=2*(zdf-GR-1)+(zdc);
	int NZs;
	
	NZs=2*zdf+GR*(zdc-4)+GR-1;
	umax=1.;
	dt=GR*1./((double)(0.5*(NZs-2.)/CFL));

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
	
	for (i=0;i<(xdc-1);i++)
	{
		for (j=0;j<(ydc-1);j++)
		{
			for (k=1;k<(zdc-1);k++)
			{
				ux=0.;uy=0.;uz=0.;dens=0.;
				for (a=0;a<19;a++)
				{
					dens+=Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);
					ux+= E0[a]*Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);
					uy+= E1[a]*Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);
					uz+= E2[a]*Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);
				}
				if ((k==1)||(k==(zdc-2)))
				  coef=0.5*GR*GR*GR;
				else
				  coef=1.*GR*GR*GR;
				momentum += coef*ux;
				cmomentum += coef*ux;
				
				denst += coef*dens;
				cdenst += coef*dens;
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
	cmomentum /= ((double)((xdf-1.)*(ydf-1.)*XDIM*YDIM*(NZs-2.)));
	
	denst /= ((double)((xdf-1.)*(ydf-1.)*XDIM*YDIM*(NZs-2.)));
	cdenst /= ((double)((xdf-1.)*(ydf-1.)*XDIM*YDIM*(NZs-2.)));

	MPI_Allreduce(&momentum,&tmomentum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	momentum=tmomentum;
	MPI_Allreduce(&cmomentum,&tcmomentum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	cmomentum=tcmomentum;

	MPI_Allreduce(&denst,&tdens,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	denst=tdens;
	MPI_Allreduce(&cdenst,&tcdenst,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	cdenst=tcdenst;

	V.rho[0]=denst;
	V.rho[1]=cdenst;

	V.force[0]=denst*ubulk-momentum;
	V.force[1]=cmomentum;
//	printf("V.rho =%f V.force =%f \n",V.rho, V.force);
//	V.force=denst*ubulk-momentum;
//	V.force=ubulk-momentum/denst;

	if (pd.myrank==0)
	{
			fh=fopen("flowrate.txt","a");
			fprintf(fh,"%f %.12f\n",((double)TM+dt*cntr),momentum*((ydf-1)*(NZs-2))/denst);
			fclose(fh);
			
			fh=fopen("density.txt","a");
			fprintf(fh,"%f %.12f\n",((double)TM+dt*cntr),denst);
			fclose(fh);
	}
	
	return V;
}
