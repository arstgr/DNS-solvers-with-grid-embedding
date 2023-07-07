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

/* 2d velocity profile for tecplot or gnuplot */
POINTER output_twoD(int xdf, int ydf, int zdf, int xdc, int ydc, int zdc, int GR, POINTER V, PDATA pd)
{
	int i,j,k,a;
	FILE *myfile;
	DP *y,*zps,*zpsf;
	DP *jj,dens,*rh0, *er,*tmp;
	int NZs,NZ=2*zdf+(zdc-4);
	DP ux;
	
	NZs=2*zdf+GR*(zdc-4)+GR-1;

	y=(DP *)calloc(NZs,sizeof(DP));
	zpsf=(DP *)calloc(NZs,sizeof(DP));
	jj=(DP *)calloc(NZ,sizeof(DP));
	rh0=(DP *)calloc(NZ,sizeof(DP));
	zps=(DP *)calloc(NZ,sizeof(DP));
	er = (DP *)calloc(NZ,sizeof(DP));
	tmp= (DP *)calloc(NZ,sizeof(DP));
	
	for (k=1;k<(NZs-1);k++)
	  *(zpsf+k) = k-0.5;
	*(zpsf+NZs-1)=NZs-2;
	
	for (k=0;k<NZs;k++)
	  *(zpsf+k) /= ((double)(NZs-2.));
	
	for (k=1;k<zdf;k++)
	  *(zps+k) = k-0.5;
	for (k=2;k<(zdc-2);k++)
	  *(zps+zdf-2+k) = GR*(k-1)+(zdf-1-0.5);
	for (k=0;k<(zdf-1);k++)
	  *(zps+zdf+zdc-4+k) = k + GR*(zdc-3) + zdf-1.-0.5;
	*(zps+2*zdf+zdc-5) = 2*zdf+GR*(zdc-4)+GR-3;
	
	for (k=0;k<NZ;k++)
	  *(zps+k) /= ((double)(NZs-2.));
	

	for (k=1;k<(NZs-1);k++)
	{
		*(y+k)=-6.*ub*((k-0.5)/(NZs-2.))*((k-0.5)/(NZs-2.)-1);
		//*(y+k)=-1.5*ub*(((k-0.5)/(NZs-2.))*((k-0.5)/(NZs-2.))-1.);	
	}

	for (i=0;i<(xdf-1);i++)
	{
		for (j=0;j<(ydf-1);j++)
		{
			for (k=1;k<zdf;k++)
			{
				dens=0.;ux=0.;
				for (a=0;a<19;a++)
				{
					dens += Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
					ux += E0[a]*Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
				}
				*(jj+k) += ux/dens;
				*(rh0+k) += dens;
			}
		}
	}
	for (k=1;k<zdf;k++)
	{
		jj[k] /= ((double)((xdf-1.)*XDIM*(ydf-1.)*YDIM));
		rh0[k] /= ((double)((xdf-1.)*XDIM*(ydf-1.)*YDIM));
	}
	
	for (i=0;i<(xdc-1);i++)
	{
		for (j=0;j<(ydc-1);j++)
		{
			for (k=2;k<(zdc-2);k++)
			{
				dens=0.;ux=0.;
				for (a=0;a<19;a++)
				{
					dens += Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);
					ux += E0[a]*Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);
				}
				*(jj+zdf-2+k) += ux/dens;
				*(rh0+zdf-2+k) += dens;
			}
		}
	}
	for (k=2;k<(zdc-2);k++)
	{
		jj[zdf-2+k] /= ((double)((xdc-1.)*XDIM*(ydc-1.)*YDIM));
		rh0[zdf-2+k] /= ((double)((xdc-1.)*XDIM*(ydc-1.)*YDIM));
	}
		
	for (i=0;i<(xdf-1);i++)
	{
		for (j=0;j<(ydf-1);j++)
		{
			for (k=0;k<(zdf-1);k++)
			{
				dens=0.;ux=0.;
				for (a=0;a<19;a++)
				{
					dens+=Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
					ux += E0[a]*Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
				}
				*(jj+zdf+zdc-4+k) += ux/dens;
				*(rh0+zdf+zdc-4+k) += dens;
			}
		}
	}
	for (k=0;k<(zdf-1);k++)
	{
		jj[zdf+zdc-4+k] /= (double)((xdf-1.)*XDIM*(ydf-1.)*YDIM);
		rh0[zdf+zdc-4+k] /= (double)((xdf-1.)*XDIM*(ydf-1.)*YDIM);
	}
	
	MPI_Allreduce(jj,tmp,NZ,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	memcpy(jj,tmp,NZ*sizeof(double));
	
	MPI_Allreduce(rh0,tmp,NZ,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	memcpy(rh0,tmp,NZ*sizeof(double));
	
	
// 	for (k=0;k<NZ;k++)
// 	{
// 		*(er+k) = fabs(*(jj+k));
// 	}
	
	if (pd.myrank==0)
	{
		myfile=fopen("2dresults-exact.dat","w");
		for (k=0;k<NZs;k++)
		  fprintf(myfile,"%.12f %.12f \n",*(zpsf+k),*(y+k));
		fclose(myfile);

		myfile=fopen("2dresults-lb.dat","w");
		for (k=0;k<NZ;k++)
		  fprintf(myfile,"%.12f %.12f %.12f\n",*(zps+k),*(jj+k),*(rh0+k));
		fclose(myfile);
	
	
		myfile=fopen("2dresults-lb.bdat","w");
		fwrite(jj,sizeof(DP),NZ,myfile);
		fclose(myfile);
	}
        
	free(y);
	free(rh0);
	free(zpsf);
	free(jj);
	free(zps);
	
	free(er);
	free(tmp);
	
	return V;
}
