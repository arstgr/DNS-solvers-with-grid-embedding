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

POINTER statistics(int xdf, int ydf, int zdf, int xdc, int ydc, int zdc, POINTER V, PDATA pd, long ts, int GR)
{
	int i,j,k;
	int a;
	DP u,v,w;
	DP UFLUC,VFLUC,WFLUC,DFLUC;
	DP *UAVG,*VAVG,*WAVG,*DUMMY,*DAVG;
	DP *UTIN,*VTIN,*WTIN,*REYS,*DTIN;
	DP dens, *zps, *y, *zpsf;
	FILE *sv;
	int exty,extx;
	const int NZ=2*zdf+(zdc-4);
	const int NZs = 2*zdf+GR*(zdc-4)+GR-1;
	char fn[20];

	UAVG=(DP *)calloc(NZ,sizeof(DP));
	VAVG=(DP *)calloc(NZ,sizeof(DP));
	WAVG=(DP *)calloc(NZ,sizeof(DP));
	DAVG=(DP *)calloc(NZ,sizeof(DP));
	DUMMY=(DP *)calloc(NZ,sizeof(DP));

	UTIN=(DP *)calloc(NZ,sizeof(DP));
	VTIN=(DP *)calloc(NZ,sizeof(DP));
	WTIN=(DP *)calloc(NZ,sizeof(DP));
	REYS=(DP *)calloc(NZ,sizeof(DP));
	DTIN=(DP *)calloc(NZ,sizeof(DP));
	
	y = (DP *)calloc(NZs,sizeof(DP));
	zpsf = (DP *)calloc(NZs,sizeof(DP));
	
	extx = xdf-1;
	exty = ydf-1;

	for (i=0;i<extx;i++)
	{
		for (j=0;j<exty;j++)
		{
			for (k=1;k<zdf;k++)
			{
					u=0.;v=0.;w=0.; dens=0.;
					for (a=0;a<19;a++)
					{
						dens += Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
						u += E0[a]*Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
						v += E1[a]*Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
						w += E2[a]*Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
					}

					*(UAVG+k) += (u/dens);
					*(VAVG+k) += (v/dens);
					*(WAVG+k) += (w/dens);
					*(DAVG+k) += dens;
			}
		}
	}

	for (k=0;k<zdf;k++)
	{
		*(UAVG+k) /= (double)(extx*exty*XDIM*YDIM);
		*(VAVG+k) /= (double)(extx*exty*XDIM*YDIM);
		*(WAVG+k) /= (double)(extx*exty*XDIM*YDIM);
		*(DAVG+k) /= (double)(extx*exty*XDIM*YDIM);
	}
	
	extx = xdc-1;
	exty = ydc-1;

	for (i=0;i<extx;i++)
	{
		for (j=0;j<exty;j++)
		{
			for (k=2;k<(zdc-2);k++)
			{
					u=0.;v=0.;w=0.; dens=0.;
					for (a=0;a<19;a++)
					{
						dens += Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);
						u += E0[a]*Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);
						v += E1[a]*Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);
						w += E2[a]*Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);
					}

					*(UAVG+zdf-2+k) += (u/dens);
					*(VAVG+zdf-2+k) += (v/dens);
					*(WAVG+zdf-2+k) += (w/dens);
					*(DAVG+zdf-2+k) += dens;
			}
		}
	}

	for (k=2;k<(zdc-2);k++)
	{
		*(UAVG+zdf-2+k) /= (double)(extx*exty*XDIM*YDIM);
		*(VAVG+zdf-2+k) /= (double)(extx*exty*XDIM*YDIM);
		*(WAVG+zdf-2+k) /= (double)(extx*exty*XDIM*YDIM);
		*(DAVG+zdf-2+k) /= (double)(extx*exty*XDIM*YDIM);
	}
	
	extx = xdf-1;
	exty = ydf-1;

	for (i=0;i<extx;i++)
	{
		for (j=0;j<exty;j++)
		{
			for (k=0;k<(zdf-1);k++)
			{
					u=0.;v=0.;w=0.; dens=0.;
					for (a=0;a<19;a++)
					{
						dens += Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
						u += E0[a]*Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
						v += E1[a]*Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
						w += E2[a]*Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
					}

					*(UAVG+zdf+zdc-4+k) += (u/dens);
					*(VAVG+zdf+zdc-4+k) += (v/dens);
					*(WAVG+zdf+zdc-4+k) += (w/dens);
					*(DAVG+zdf+zdc-4+k) += dens;
			}
		}
	}

	for (k=0;k<zdf;k++)
	{
		*(UAVG+zdf+zdc-4+k) /= (double)(extx*exty*XDIM*YDIM);
		*(VAVG+zdf+zdc-4+k) /= (double)(extx*exty*XDIM*YDIM);
		*(WAVG+zdf+zdc-4+k) /= (double)(extx*exty*XDIM*YDIM);
		*(DAVG+zdf+zdc-4+k) /= (double)(extx*exty*XDIM*YDIM);
	}

	MPI_Allreduce(UAVG,DUMMY,NZ,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	memcpy(UAVG,DUMMY,NZ*sizeof(DP));

	MPI_Allreduce(VAVG,DUMMY,NZ,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	memcpy(VAVG,DUMMY,NZ*sizeof(DP));

	MPI_Allreduce(WAVG,DUMMY,NZ,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	memcpy(WAVG,DUMMY,NZ*sizeof(DP));

	MPI_Allreduce(DAVG,DUMMY,NZ,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	memcpy(DAVG,DUMMY,NZ*sizeof(DP));
	
	extx = xdf-1;
	exty = ydf-1;

	for (i=0;i<extx;i++)
	{
		for (j=0;j<exty;j++)
		{
			for (k=1;k<zdf;k++)
			{
				u=0.;v=0.;w=0.; dens=0.;
				for (a=0;a<19;a++)
				{
					dens += Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
					u += E0[a]*Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
					v += E1[a]*Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
					w += E2[a]*Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
				}

				UFLUC = (u/dens)-(*(UAVG+k));
				VFLUC = (v/dens)-(*(VAVG+k));
				WFLUC = (w/dens)-(*(WAVG+k));
				DFLUC = dens-(*(DAVG+k));
				*(UTIN+k) += UFLUC*UFLUC;
				*(VTIN+k) += VFLUC*VFLUC;
				*(WTIN+k) += WFLUC*WFLUC;
				*(REYS+k) += UFLUC*WFLUC;
				*(DTIN+k) += DFLUC*DFLUC;
			}
		}
	}

	for (k=0;k<zdf;k++)
	{
		*(UTIN+k) /= (double)(extx*exty*XDIM*YDIM);
		*(VTIN+k) /= (double)(exty*extx*XDIM*YDIM);
		*(WTIN+k) /= (double)(exty*extx*XDIM*YDIM);
		*(REYS+k) /= (double)(extx*exty*XDIM*YDIM);
		*(DTIN+k) /= (double)(extx*exty*XDIM*YDIM);
	}
	
	extx = xdc-1;
	exty = ydc-1;

	for (i=0;i<extx;i++)
	{
		for (j=0;j<exty;j++)
		{
			for (k=2;k<(zdc-2);k++)
			{
				u=0.;v=0.;w=0.; dens=0.;
				for (a=0;a<19;a++)
				{
					dens += Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);
					u += E0[a]*Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);
					v += E1[a]*Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);
					w += E2[a]*Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);
				}

				UFLUC = (u/dens)-(*(UAVG+zdf-2+k));
				VFLUC = (v/dens)-(*(VAVG+zdf-2+k));
				WFLUC = (w/dens)-(*(WAVG+zdf-2+k));
				DFLUC = dens-(*(DAVG+zdf-2+k));
				*(UTIN+zdf-2+k) += UFLUC*UFLUC;
				*(VTIN+zdf-2+k) += VFLUC*VFLUC;
				*(WTIN+zdf-2+k) += WFLUC*WFLUC;
				*(REYS+zdf-2+k) += UFLUC*WFLUC;
				*(DTIN+zdf-2+k) += DFLUC*DFLUC;
			}
		}
	}

	for (k=2;k<(zdc-2);k++)
	{
		*(UTIN+zdf-2+k) /= (double)(extx*exty*XDIM*YDIM);
		*(VTIN+zdf-2+k) /= (double)(exty*extx*XDIM*YDIM);
		*(WTIN+zdf-2+k) /= (double)(exty*extx*XDIM*YDIM);
		*(REYS+zdf-2+k) /= (double)(extx*exty*XDIM*YDIM);
		*(DTIN+zdf-2+k) /= (double)(extx*exty*XDIM*YDIM);
	}
	
	extx = xdf-1;
	exty = ydf-1;

	for (i=0;i<extx;i++)
	{
		for (j=0;j<exty;j++)
		{
			for (k=0;k<(zdf-1);k++)
			{
				u=0.;v=0.;w=0.; dens=0.;
				for (a=0;a<19;a++)
				{
					dens += Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
					u += E0[a]*Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
					v += E1[a]*Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
					w += E2[a]*Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
				}

				UFLUC = (u/dens)-(*(UAVG+zdf+zdc-4+k));
				VFLUC = (v/dens)-(*(VAVG+zdf+zdc-4+k));
				WFLUC = (w/dens)-(*(WAVG+zdf+zdc-4+k));
				DFLUC = dens-(*(DAVG+zdf+zdc-4+k));
				*(UTIN+zdf+zdc-4+k) += UFLUC*UFLUC;
				*(VTIN+zdf+zdc-4+k) += VFLUC*VFLUC;
				*(WTIN+zdf+zdc-4+k) += WFLUC*WFLUC;
				*(REYS+zdf+zdc-4+k) += UFLUC*WFLUC;
				*(DTIN+zdf+zdc-4+k) += DFLUC*DFLUC;
			}
		}
	}

	for (k=0;k<zdf;k++)
	{
		*(UTIN+zdf+zdc-4+k) /= (double)(extx*exty*XDIM*YDIM);
		*(VTIN+zdf+zdc-4+k) /= (double)(exty*extx*XDIM*YDIM);
		*(WTIN+zdf+zdc-4+k) /= (double)(exty*extx*XDIM*YDIM);
		*(REYS+zdf+zdc-4+k) /= (double)(extx*exty*XDIM*YDIM);
		*(DTIN+zdf+zdc-4+k) /= (double)(extx*exty*XDIM*YDIM);
	}
	
	MPI_Allreduce(UTIN,DUMMY,NZ,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	memcpy(UTIN,DUMMY,NZ*sizeof(DP));

	MPI_Allreduce(VTIN,DUMMY,NZ,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	memcpy(VTIN,DUMMY,NZ*sizeof(DP));

	MPI_Allreduce(WTIN,DUMMY,NZ,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	memcpy(WTIN,DUMMY,NZ*sizeof(DP));

	MPI_Allreduce(REYS,DUMMY,NZ,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	memcpy(REYS,DUMMY,NZ*sizeof(DP));

	MPI_Allreduce(DTIN,DUMMY,NZ,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	memcpy(DTIN,DUMMY,NZ*sizeof(DP));

	for (k=0;k<NZ;k++)
	{
		*(UTIN+k) = sqrt((*(UTIN+k)));
		*(VTIN+k) = sqrt((*(VTIN+k)));
		*(WTIN+k) = sqrt((*(WTIN+k)));
		*(DTIN+k) = sqrt((*(DTIN+k)));
	}


	if (pd.myrank==0)
	{
		sprintf(fn,"turb-fld.%.3d.%.3ld",pd.myrank,ts);
		sv=fopen(fn,"wb");
		fwrite(UAVG,sizeof(DP),NZ,sv);
		fwrite(VAVG,sizeof(DP),NZ,sv);
		fwrite(WAVG,sizeof(DP),NZ,sv);
		fwrite(UTIN,sizeof(DP),NZ,sv);
		fwrite(VTIN,sizeof(DP),NZ,sv);
		fwrite(WTIN,sizeof(DP),NZ,sv);
		fwrite(REYS,sizeof(DP),NZ,sv);
		fwrite(DAVG,sizeof(DP),NZ,sv);
		fwrite(DTIN,sizeof(DP),NZ,sv);
		fclose(sv);
	
// Test print function
		zps = (double *)calloc(NZ,sizeof(double));
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
		
		sv=fopen("uavg-temp.dat","w");
		for (k=0;k<NZ;k++)
		  fprintf(sv,"%.12f %.12f\n",*(zps+k),*(UAVG+k));
		fclose(sv);
		
		for (k=1;k<(NZs-1);k++)
		  *(zpsf+k) = k-0.5;
		*(zpsf+NZs-1)=NZs-2;
	
		for (k=0;k<NZs;k++)
		  *(zpsf+k) /= ((double)(NZs-2.));
		
		sv=fopen("uavg-exact-lam.dat","w");
		for (k=0;k<NZs;k++)
		  fprintf(sv,"%.12f %.12f \n",*(zpsf+k),*(y+k));
		fclose(sv);
		
		
		
		free(zps);
// End of test print function
	}

	free(UAVG);
	free(VAVG);
	free(WAVG);
	free(DAVG);

	free(UTIN);
	free(VTIN);
	free(WTIN);
	free(REYS);
	free(DTIN);

	free(DUMMY);
	
	free(y);
	free(zpsf);
	
//	fprintf(stderr,"in function statistics, c error handler is %s\n",strerror(errno));

	return V;
}