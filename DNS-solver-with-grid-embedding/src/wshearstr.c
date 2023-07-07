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

POINTER wshearstr(int xdf, int ydf, int zdf, int zdc, int GR, POINTER V, PDATA pd, DP tauf, int TM, long cntr)
{
	int a,mm;
	int k=1,i,j;
	DP stress=0.,ux=0.,dens=0.;
	DP stup=0.,stdown=0.,dummy=0.;
	DP mu,davg=0.;
	FILE *fn;
	char fl[30];
	double dt,Re,umax;
	double ubulk=ub;
	DP ux1,ux2,ux3;
	DP pgrad,pgradt;
	const int L[5]={2,8,9,13,14};
	const int R[5]={1,7,10,11,12};
	const int NZs=2*zdf+GR*(zdc-4)+GR-1;

	umax=1.;//uts*(2.5*log(Re)+5.5);
	dt=GR*1./((double)(0.5*(NZs-2.)/CFL));

	mu=(1./6.)*(2./tauf-1);

	for (i=0;i<(xdf-1);i++)
	{
		for (j=0;j<(ydf-1);j++)
		{
			k=1;
			ux1=0.;dens=0.;
			for (a=0;a<19;a++)
			{
				dens+=Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
				ux1+= E0[a]*Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
			}
			ux1 = ux1/dens;
			davg += dens;
			k=2;
			ux2=0.;dens=0.;
			for (a=0;a<19;a++)
			{
				dens+=Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
				ux2+= E0[a]*Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
			}
			ux2 = ux2/dens;

			k=3;
			ux3=0.;dens=0.;
			for (a=0;a<19;a++)
			{
				dens+=Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
				ux3+= E0[a]*Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
			}
			ux3 = ux3/dens;

			stdown += (-2.*ux1+3.*ux2-ux3);
		}
	}

	davg /= ((double)((xdf-1.)*(ydf-1.)*YDIM*XDIM));
	stdown /= ((double)((xdf-1.)*(ydf-1.)*XDIM*YDIM));

	MPI_Allreduce(&stdown,&dummy,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	stdown=dummy;

	MPI_Allreduce(&davg,&dummy,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	davg=dummy;

	stdown *= mu*davg;

	k=zdf-2;
	davg=0.;
	for (i=0;i<(xdf-1);i++)
	{
		for (j=0;j<(ydf-1);j++)
		{
			k=zdf-2;
			ux1=0.;dens=0.;
			for (a=0;a<19;a++)
			{
				dens+=Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
				ux1+= E0[a]*Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
			}
			ux1 = ux1/dens;
			k=zdf-3;
			ux2=0.;dens=0.;
			for (a=0;a<19;a++)
			{
				dens+=Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
				ux2+= E0[a]*Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
			}
			ux2 = ux2/dens;

			k=zdf-4;
			ux3=0.;dens=0.;
			for (a=0;a<19;a++)
			{
				dens+=Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
				ux3+= E0[a]*Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
			}
			ux3 = ux3/dens;

			stup += -(2.*ux1-3.*ux2+ux3);
			davg += dens;
		}
	}

	davg /= ((double)((xdf-1.)*(ydf-1.)*XDIM*YDIM));
	stup /= ((double)((xdf-1.)*(ydf-1.)*XDIM*YDIM));
	dummy=0.;

	MPI_Allreduce(&stup,&dummy,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	stup=dummy;

	MPI_Allreduce(&davg,&dummy,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	davg=dummy;

	pgrad=V.force[0];

	pgrad = pgrad*(NZs-2.)/(2.*GR);

	stup *= mu*davg;
	
//	fprintf(stderr, "in function wshearstr, c error handler is %s\n",strerror(errno));

	stress=0.5*(stup+stdown);
	if (isnan(stress))
	{
		fprintf(stderr,"Error: NAN in wal shear stress\n");
		MPI_Abort(MPI_COMM_WORLD, mpi_errors);
	}

	if (pd.myrank==0)
	{
		sprintf(fl,"str-up.txt");
		fn=fopen(fl,"a");
		fprintf(fn,"%f %.12f\n",((double)TM+dt*cntr),stup);
		fclose(fn);

		sprintf(fl,"str-do.txt");
		fn=fopen(fl,"a");
		fprintf(fn,"%f %.12f\n",((double)TM+dt*cntr),stdown);
		fclose(fn);

		sprintf(fl,"w-str.txt");
		fn=fopen(fl,"a");
		fprintf(fn,"%f %.12f\n",((double)TM+dt*cntr),stress);
		fclose(fn);

		sprintf(fl,"p-grad.txt");
		fn=fopen(fl,"a");
		fprintf(fn,"%f %.12f\n",((double)TM+dt*cntr),pgrad);
                fclose(fn);
	}

	return V;
}