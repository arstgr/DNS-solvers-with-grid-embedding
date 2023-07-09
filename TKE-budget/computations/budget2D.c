/************************************************************************
 * By Amirreza Rastegari                                                *
 *   arstgri@gmail.com                                                  *
 *                                                                      *
 * Assumes a Multigrid data input                                       *
 * calculates rms, vorticity and mean velocities                        *
 * Calculates the energy budget2D                                       *
 * Skewness and Kurtosis                                                *
 * Assumes a 2D mean velocity i.e. U(y,z)                               *
 ************************************************************************/
#include "definitions.h"

POINTER budget2D(int xdc, int ydc, int zdc, int xdf, int ydf, int zdf, int GR, POINTER V, PDATA pd, long ts)
{
	int i,j,k,a,q,c;
	DP ufl,vfl,wfl,dfl,dens;
	int xd=0,yd=0,zd=0;
	int NZ=2*zdf+(zdc-4);
	int NZs=2*zdf+GR*(zdc-4)+GR-1;
	
	DP *UAVG,*VAVG,*WAVG,*DUMMY;
	DP *UTIN,*VTIN,*WTIN,*DTIN;
	double dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz;
	double *uu,*uv,*uw,*vv,*vw,*ww,*dd;
	DP *DUMMYXZ,*DAVGXZ,*DAVG;
	DP siu,siv,siw,sip;
	double *SKEWu,*SKEWv,*SKEWw,*SKEWp,*SKEWuw,*SKEWuv,*SKEWvw;
	double *KURTu,*KURTv,*KURTw,*KURTp,*KURTuw,*KURTuv,*KURTvw;
	double ox,oy,oz;
	double u,v,w;
	double *wx,*wy,*wz;
	double *UW;
	double dUdz[NZ];
	double *KE,KEt;
	double *Pii,*Puu,*Pvv,*Pww,*Puw,*Puv,*Pvw;
	double *Tiir,*Tiipai,*Tiis,*Epii;
	double *Paiuu,*Paivv,*Paiww,*Paiuw,*Paiuv,*Paivw,*Epuu,*Epvv,*Epww,*Epuw,*Epuv,*Epvw;
	double *Tiirt,*Tuur,*Tvvr,*Twwr,*Tuwr,*Tuvr,*Tvwr;
	double *Tuupai,*Tvvpai,*Twwpai,*Tuwpai,*Tuvpai,*Tvwpai;
	double *Tuus,*Tvvs,*Twws,*Tuws,*Tuvs,*Tvws;
	double tmp;
	double coef;
	
	DP pu[NZ],pv[NZ],pw[NZ],uuw[NZ],vvw[NZ],www[NZ],uww[NZ];
	DP usuw[NZ],vsvw[NZ],wsww[NZ],wsuw[NZ],usww[NZ];

	DP nu,tauc,tauf;
	FILE *sv;
	int jp,jn,in,ip;
	int jpp,jnn,inn,ipp;
	int num=0,temp=0,upn;
	int extx,exty;
	char fn[50],fh[100],fht[100];
	MPI_Status status;

	tauf=(3./Rebs)*ub*(NZs-2.)+0.5;
	tauc = (1./((double)GR))*(3./Rebs)*ub*(NZs-2.)+0.5;
	nu=(1./6.)*(2.*tauf-1.);
	coef = (1./((double)GR));
	
	
	zd=NZ;
	V.vel=NULL;
	V.den=NULL;
	
	UAVG=(DP *)calloc(zd,sizeof(DP));
	VAVG=(DP *)calloc(zd,sizeof(DP));
	WAVG=(DP *)calloc(zd,sizeof(DP));
	DAVG=(DP *)calloc(zd,sizeof(DP));
	DUMMY=(DP *)calloc(zd,sizeof(DP));
	UW=(DP *)calloc(zd,sizeof(DP));	

	UTIN=(DP *)calloc(zd,sizeof(DP));
	VTIN=(DP *)calloc(zd,sizeof(DP));
	WTIN=(DP *)calloc(zd,sizeof(DP));
	DTIN=(DP *)calloc(zd,sizeof(DP));
	
	uu=(DP *)calloc(zd,sizeof(DP));
	uv=(DP *)calloc(zd,sizeof(DP));
	uw=(DP *)calloc(zd,sizeof(DP));
	vv=(DP *)calloc(zd,sizeof(DP));
	vw=(DP *)calloc(zd,sizeof(DP));
	ww=(DP *)calloc(zd,sizeof(DP));
	dd=(DP *)calloc(zd,sizeof(DP));
	
	wx=(DP *)calloc(zd,sizeof(DP));
	wy=(DP *)calloc(zd,sizeof(DP));
	wz=(DP *)calloc(zd,sizeof(DP));
	
	SKEWu=(DP *)calloc(zd,sizeof(DP));
	SKEWv=(DP *)calloc(zd,sizeof(DP));
	SKEWw=(DP *)calloc(zd,sizeof(DP));
	SKEWp=(DP *)calloc(zd,sizeof(DP));
	SKEWuw=(DP *)calloc(zd,sizeof(DP));
	SKEWuv=(DP *)calloc(zd,sizeof(DP));
	SKEWvw=(DP *)calloc(zd,sizeof(DP));
	
	KURTu=(DP *)calloc(zd,sizeof(DP));
	KURTv=(DP *)calloc(zd,sizeof(DP));
	KURTw=(DP *)calloc(zd,sizeof(DP));
	KURTp=(DP *)calloc(zd,sizeof(DP));
	KURTuw=(DP *)calloc(zd,sizeof(DP));
	KURTvw=(DP *)calloc(zd,sizeof(DP));
	KURTuv=(DP *)calloc(zd,sizeof(DP));

	KE=(DP *)calloc(zd,sizeof(DP));

	Paiuu=(DP *)calloc(zd,sizeof(DP));
	Paivv=(DP *)calloc(zd,sizeof(DP));
	Paiww=(DP *)calloc(zd,sizeof(DP));
	Paiuw=(DP *)calloc(zd,sizeof(DP));
	Paiuv=(DP *)calloc(zd,sizeof(DP));
	Paivw=(DP *)calloc(zd,sizeof(DP));

	Tiirt=(DP *)calloc(zd,sizeof(DP));
	Tuur=(DP *)calloc(zd,sizeof(DP));
	Tvvr=(DP *)calloc(zd,sizeof(DP));
	Twwr=(DP *)calloc(zd,sizeof(DP));
	Tuwr=(DP *)calloc(zd,sizeof(DP));
	Tuvr=(DP *)calloc(zd,sizeof(DP));
	Tvwr=(DP *)calloc(zd,sizeof(DP));

	Tuupai=(DP *)calloc(zd,sizeof(DP));
	Tvvpai=(DP *)calloc(zd,sizeof(DP));
	Twwpai=(DP *)calloc(zd,sizeof(DP));
	Tuwpai=(DP *)calloc(zd,sizeof(DP));
	Tuvpai=(DP *)calloc(zd,sizeof(DP));
	Tvwpai=(DP *)calloc(zd,sizeof(DP));

	Tuus=(DP *)calloc(zd,sizeof(DP));
	Tvvs=(DP *)calloc(zd,sizeof(DP));
	Twws=(DP *)calloc(zd,sizeof(DP));
	Tuws=(DP *)calloc(zd,sizeof(DP));
	Tuvs=(DP *)calloc(zd,sizeof(DP));
	Tvws=(DP *)calloc(zd,sizeof(DP));

	Epuu=(DP *)calloc(zd,sizeof(DP));
	Epvv=(DP *)calloc(zd,sizeof(DP));
	Epww=(DP *)calloc(zd,sizeof(DP));
	Epuw=(DP *)calloc(zd,sizeof(DP));
	Epuv=(DP *)calloc(zd,sizeof(DP));
	Epvw=(DP *)calloc(zd,sizeof(DP));

	Pii=(DP *)calloc(zd,sizeof(DP));
	Puu=(DP *)calloc(zd,sizeof(DP));
	Pvv=(DP *)calloc(zd,sizeof(DP));
	Pww=(DP *)calloc(zd,sizeof(DP));
	Puw=(DP *)calloc(zd,sizeof(DP));
	Puv=(DP *)calloc(zd,sizeof(DP));
	Pvw=(DP *)calloc(zd,sizeof(DP));
	
	Tiir=(DP *)calloc(zd,sizeof(DP));
	Tiipai=(DP *)calloc(zd,sizeof(DP));
	Tiis=(DP *)calloc(zd,sizeof(DP));
	Epii=(DP *)calloc(zd,sizeof(DP));
	
	for (k=0;k<NZ;k++)
	{
		pu[k]=0.;
		pv[k]=0.;
		pw[k]=0.;
		uuw[k]=0.;
		vvw[k]=0.;
		www[k]=0.;
		uww[k]=0.;
		
		usuw[k]=0.;
		vsvw[k]=0.;
		wsww[k]=0.;
		wsuw[k]=0.;
		usww[k]=0.;
		
		dUdz[k]=0.;
	}
	
	extx = xdf-1;
	exty = ydf;
		
	for (j=0;j<exty;j++)
	{
		for (k=0;k<zdf;k++)
		{
			for (i=0;i<extx;i++)
			{
				*(UAVG+k) +=Fb(V.velfl,i,j,k,0,xdf,ydf,zdf,3);
				*(VAVG+k) +=Fb(V.velfl,i,j,k,1,xdf,ydf,zdf,3);
				*(WAVG+k) +=Fb(V.velfl,i,j,k,2,xdf,ydf,zdf,3);
				*(DAVG+k) +=(1./3.)*Fb(V.denfl,i,j,k,0,xdf,ydf,zdf,1);
			}
		}
	}
	
	for (k=0;k<zdf;k++)
	{
		*(UAVG+k) /= (double)(extx*exty*pd.numproc);
		*(VAVG+k) /= (double)(extx*exty*pd.numproc);
		*(WAVG+k) /= (double)(extx*exty*pd.numproc);
		*(DAVG+k) /= (double)(extx*exty*pd.numproc);
	}
	
	extx = xdc-1;
	exty = ydc;
		
	for (j=0;j<exty;j++)
	{
		for (k=2;k<(zdc-2);k++)
		{
			for (i=0;i<extx;i++)
			{
				UAVG[zdf-2+k] +=Fb(V.velc,i,j,k,0,xdc,ydc,zdc,3);
				VAVG[zdf-2+k] +=Fb(V.velc,i,j,k,1,xdc,ydc,zdc,3);
				WAVG[zdf-2+k] +=Fb(V.velc,i,j,k,2,xdc,ydc,zdc,3);
				DAVG[zdf-2+k] +=(1./3.)*Fb(V.denc,i,j,k,0,xdc,ydc,zdc,1);
			}
		}
	}
	for (k=2;k<(zdc-2);k++)
	{
		*(UAVG+zdf-2+k) /= (double)(extx*exty*pd.numproc);
		*(VAVG+zdf-2+k) /= (double)(extx*exty*pd.numproc);
		*(WAVG+zdf-2+k) /= (double)(extx*exty*pd.numproc);
		*(DAVG+zdf-2+k) /= (double)(extx*exty*pd.numproc);
	}
	
	extx = xdf-1;
	exty = ydf;
		
	for (j=0;j<exty;j++)
	{
		for (k=0;k<zdf;k++)
		{
			for (i=0;i<extx;i++)
			{
				UAVG[zdf+zdc-4+k]+=Fb(V.velfu,i,j,k,0,xdf,ydf,zdf,3);
				VAVG[zdf+zdc-4+k]+=Fb(V.velfu,i,j,k,1,xdf,ydf,zdf,3);
				WAVG[zdf+zdc-4+k]+=Fb(V.velfu,i,j,k,2,xdf,ydf,zdf,3);
				DAVG[zdf+zdc-4+k]+=(1./3.)*Fb(V.denfu,i,j,k,0,xdf,ydf,zdf,1);
			}
		}
	}
	for (k=0;k<zdf;k++)
	{
		*(UAVG+zdf+zdc-4+k) /= (double)(extx*exty*pd.numproc);
		*(VAVG+zdf+zdc-4+k) /= (double)(extx*exty*pd.numproc);
		*(WAVG+zdf+zdc-4+k) /= (double)(extx*exty*pd.numproc);
		*(DAVG+zdf+zdc-4+k) /= (double)(extx*exty*pd.numproc);
	}
	
	
	MPI_Reduce(UAVG,DUMMY,NZ,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,NZ,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(UAVG,DUMMY,NZ*sizeof(DP));	
	
	MPI_Reduce(VAVG,DUMMY,NZ,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,NZ,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(VAVG,DUMMY,NZ*sizeof(DP));
	
	MPI_Reduce(WAVG,DUMMY,NZ,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,NZ,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(WAVG,DUMMY,NZ*sizeof(DP));
	
	MPI_Reduce(DAVG,DUMMY,NZ,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,NZ,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(DAVG,DUMMY,NZ*sizeof(DP));
	
	extx=xdf-1;
	exty=ydf;
			
	for (i=-2;i<(extx+2);i++)
	{
		for (j=0;j<exty;j++)
		{
			for (k=0;k<zdf;k++)
			{
				Fb(V.velfl,i,j,k,0,xdf,ydf,zdf,3) -= UAVG[k];
				Fb(V.velfl,i,j,k,1,xdf,ydf,zdf,3) -= VAVG[k];
				Fb(V.velfl,i,j,k,2,xdf,ydf,zdf,3) -= WAVG[k];
			}
		}
	}
	
	extx=xdc-1;
	exty=ydc;
			
	for (i=-2;i<(extx+2);i++)
	{
		for (j=0;j<exty;j++)
		{
			for (k=2;k<(zdc-2);k++)
			{
				Fb(V.velc,i,j,k,0,xdc,ydc,zdc,3) -= UAVG[zdf-2+k];
				Fb(V.velc,i,j,k,1,xdc,ydc,zdc,3) -= VAVG[zdf-2+k];
				Fb(V.velc,i,j,k,2,xdc,ydc,zdc,3) -= WAVG[zdf-2+k];
			}
		}
	}
	
// 	for (i=-2;i<(xdf+1);i++)
// 	{
// 		for (j=0;j<ydf;j++)
// 		{
// 			Fb(V.velc,(i/GR),(j/GR),0,0,xdc,ydc,zdc,3) -= UAVG[zdf-1-GR];
// 			Fb(V.velc,(i/GR),(j/GR),0,1,xdc,ydc,zdc,3) -= VAVG[zdf-1-GR];
// 			Fb(V.velc,(i/GR),(j/GR),0,2,xdc,ydc,zdc,3) -= WAVG[zdf-1-GR];
// 			
// 			Fb(V.velc,(i/GR),(j/GR),1,0,xdc,ydc,zdc,3) -= UAVG[zdf-1];
// 			Fb(V.velc,(i/GR),(j/GR),1,1,xdc,ydc,zdc,3) -= VAVG[zdf-1];
// 			Fb(V.velc,(i/GR),(j/GR),1,2,xdc,ydc,zdc,3) -= WAVG[zdf-1];
// 			
// 			Fb(V.velc,(i/GR),(j/GR),(zdc-2),0,xdc,ydc,zdc,3) -= UAVG[zdf+zdc-4];
// 			Fb(V.velc,(i/GR),(j/GR),(zdc-2),1,xdc,ydc,zdc,3) -= VAVG[zdf+zdc-4];
// 			Fb(V.velc,(i/GR),(j/GR),(zdc-2),2,xdc,ydc,zdc,3) -= WAVG[zdf+zdc-4];
// 			
// 			Fb(V.velc,(i/GR),(j/GR),(zdc-1),0,xdc,ydc,zdc,3) -= UAVG[zdf+zdc-4+GR];
// 			Fb(V.velc,(i/GR),(j/GR),(zdc-1),1,xdc,ydc,zdc,3) -= VAVG[zdf+zdc-4+GR];
// 			Fb(V.velc,(i/GR),(j/GR),(zdc-1),2,xdc,ydc,zdc,3) -= WAVG[zdf+zdc-4+GR];
// 		}
// 	}
	
	extx=xdf-1;
	exty=ydf;
			
	for (i=-2;i<(extx+2);i++)
	{
		for (j=0;j<exty;j++)
		{
			for (k=0;k<zdf;k++)
			{
				Fb(V.velfu,i,j,k,0,xdf,ydf,zdf,3) -= UAVG[zdf+zdc-4+k];
				Fb(V.velfu,i,j,k,1,xdf,ydf,zdf,3) -= VAVG[zdf+zdc-4+k];
				Fb(V.velfu,i,j,k,2,xdf,ydf,zdf,3) -= WAVG[zdf+zdc-4+k];
			}
		}
	}
	
	extx=xdf-1;
	exty=ydf;
	xd=xdf;yd=ydf;zd=zdf;
	
	
	for (k=1;k<zd;k++)
	{
		
		if (k==1)	
		{
			tmp = (2./3.)*(UAVG[k+3]-UAVG[k+1])-(1./12.)*(UAVG[k+4]-UAVG[k]);
			dUdz[k] = -(7./3.)*UAVG[k]+6.*UAVG[k+1]-3.*UAVG[k+2]-(2./3.)*UAVG[k+3]+3.*tmp;
		}
		else if (k==(zd-1))
		{
// 			tmp = (2./3.)*(UAVG[k-1]-UAVG[k-3])-(1./12.)*(UAVG[k]-UAVG[k-4]);
// 			dUdz[k] = (7./3.)*UAVG[k]-6.*UAVG[k-1]+3.*UAVG[k-2]+(2./3.)*UAVG[k-3]+3.*tmp;
			dUdz[k] = (coef*(UAVG[k+1]-UAVG[k])+GR*(UAVG[k]-UAVG[k-1]))*1./(1.+GR);
		}
		else if (k==2)
		{
// 			tmp = (2./3.)*(UAVG[k+2]-UAVG[k])-(1./12.)*(UAVG[k+3]-UAVG[k-1]);
// 			dUdz[k] = -(1./6.)*UAVG[k-1]-(3./2.)*UAVG[k]+(3./2.)*UAVG[k+1]+(1./6.)*UAVG[k+2]-tmp;
			dUdz[k] = 0.5*(UAVG[k+1]-UAVG[k-1]);
		}
		else if (k==(zd-2))
		{
// 			tmp = (2./3.)*(UAVG[k]-UAVG[k-2])-(1./12.)*(UAVG[k+1]-UAVG[k-3]);
// 			dUdz[k] = (1./6.)*UAVG[k+1]+(3./2.)*UAVG[k]-(3./2.)*UAVG[k-1]-(1./6.)*UAVG[k-2]-tmp;
			dUdz[k] = 0.5*(UAVG[k+1]-UAVG[k-1]);
		}
		else
		{
			dUdz[k] = (2./3.)*(UAVG[k+1]-UAVG[k-1])-(1./12.)*(UAVG[k+2]-UAVG[k-2]);
		}
	}
	
	extx=xdc-1;
	exty=ydc;
	xd=xdc;yd=ydc;zd=zdc;
	
	for (k=2;k<(zd-2);k++)
	{	
		if (k==2)
		{
// 			tmp = (2./3.)*(UAVG[zdf-2+k+3]-UAVG[zdf-2+k+1])-(1./12.)*(UAVG[zdf-2+k+4]-UAVG[zdf-2+k]);
// 			dUdz[zdf-2+k] = coef*(-(7./3.)*UAVG[zdf-2+k]+6.*UAVG[zdf-2+k+1]-3.*UAVG[zdf-2+k+2]-(2./3.)*UAVG[zdf-2+k+3]+3.*tmp);
			dUdz[zdf-2+k] = 0.5*coef*(UAVG[zdf-2+k+1]-UAVG[zdf-2+k-1]);
		}
		else if (k==3)
		{
// 			tmp = (2./3.)*(UAVG[zdf-2+k+2]-UAVG[zdf-2+k])-(1./12.)*(UAVG[zdf-2+k+3]-UAVG[zdf-2+k-1]);
// 			dUdz[zdf-2+k] = coef*(-(1./6.)*UAVG[zdf-2+k-1]-(3./2.)*UAVG[zdf-2+k]+(3./2.)*UAVG[zdf-2+k+1]+(1./6.)*UAVG[zdf-2+k+2]-tmp);
			dUdz[zdf-2+k] = 0.5*coef*(UAVG[zdf-2+k+1]-UAVG[zdf-2+k-1]);
		}
		else if (k==(zdc-4))
		{
// 			tmp = (2./3.)*(UAVG[zdf-2+k]-UAVG[zdf-2+k-2])-(1./12.)*(UAVG[zdf-2+k+1]-UAVG[zdf-2+k-3]);
// 			dUdz[zdf-2+k] = coef*((1./6.)*UAVG[zdf-2+k+1]+(3./2.)*UAVG[zdf-2+k]-(3./2.)*UAVG[zdf-2+k-1]-(1./6.)*UAVG[zdf-2+k-2]-tmp);
			dUdz[zdf-2+k] = 0.5*coef*(UAVG[zdf-2+k+1]-UAVG[zdf-2+k-1]);
		}
		else if (k==(zdc-3))
		{
// 			tmp = (2./3.)*(UAVG[zdf-2+k-1]-UAVG[zdf-2+k-3])-(1./12.)*(UAVG[zdf-2+k]-UAVG[zdf-2+k-4]);
// 			dUdz[zdf-2+k] = coef*((7./3.)*UAVG[zdf-2+k]-6.*UAVG[zdf-2+k-1]+3.*UAVG[zdf-2+k-2]+(2./3.)*UAVG[zdf-2+k-3]+3.*tmp);
			dUdz[zdf-2+k] = 0.5*coef*(UAVG[zdf-2+k+1]-UAVG[zdf-2+k-1]);
		}
		else
		{
			dUdz[zdf-2+k] = coef*((2./3.)*(UAVG[zdf-2+k+1]-UAVG[zdf-2+k-1])-(1./12.)*(UAVG[zdf-2+k+2]-UAVG[zdf-2+k-2]));
		}
	}
	
	extx=xdf-1;
	exty=ydf;
	xd=xdf;yd=ydf;zd=zdf;
	
	for (k=0;k<(zd-1);k++)
	{
		
		if (k==0)	
		{
// 			tmp = (2./3.)*(UAVG[zdf+zdc-4+k+3]-UAVG[zdf+zdc-4+k+1])-(1./12.)*(UAVG[zdf+zdc-4+k+4]-UAVG[zdf+zdc-4+k]);
// 			dUdz[zdf+zdc-4+k] = -(7./3.)*UAVG[zdf+zdc-4+k]+6.*UAVG[zdf+zdc-4+k+1]-3.*UAVG[zdf+zdc-4+k+2]-(2./3.)*UAVG[zdf+zdc-4+k+3]+3.*tmp;
			dUdz[zdf+zdc-4+k] = (GR*(UAVG[zdf+zdc-4+k+1]-UAVG[zdf+zdc-4+k])+coef*(UAVG[zdf+zdc-4+k]-UAVG[zdf+zdc-4+k-1]))*1./(1.+GR);
		}
		else if (k==(zd-2))
		{
			tmp = (2./3.)*(UAVG[zdf+zdc-4+k-1]-UAVG[zdf+zdc-4+k-3])-(1./12.)*(UAVG[zdf+zdc-4+k]-UAVG[zdf+zdc-4+k-4]);
			dUdz[zdf+zdc-4+k] = (7./3.)*UAVG[zdf+zdc-4+k]-6.*UAVG[zdf+zdc-4+k-1]+3.*UAVG[zdf+zdc-4+k-2]+(2./3.)*UAVG[zdf+zdc-4+k-3]+3.*tmp;
		}
		else if (k==1)
		{
// 			tmp = (2./3.)*(UAVG[zdf+zdc-4+k+2]-UAVG[zdf+zdc-4+k])-(1./12.)*(UAVG[zdf+zdc-4+k+3]-UAVG[zdf+zdc-4+k-1]);
// 			dUdz[zdf+zdc-4+k] = -(1./6.)*UAVG[zdf+zdc-4+k-1]-(3./2.)*UAVG[zdf+zdc-4+k]+(3./2.)*UAVG[zdf+zdc-4+k+1]+(1./6.)*UAVG[zdf+zdc-4+k+2]-tmp;
			dUdz[zdf+zdc-4+k] = 0.5*(UAVG[zdf+zdc-4+k+1]-UAVG[zdf+zdc-4+k-1]);
		}
		else if (k==(zd-3))
		{
// 			tmp = (2./3.)*(UAVG[zdf+zdc-4+k]-UAVG[zdf+zdc-4+k-2])-(1./12.)*(UAVG[zdf+zdc-4+k+1]-UAVG[zdf+zdc-4+k-3]);
// 			dUdz[zdf+zdc-4+k] = (1./6.)*UAVG[zdf+zdc-4+k+1]+(3./2.)*UAVG[zdf+zdc-4+k]-(3./2.)*UAVG[zdf+zdc-4+k-1]-(1./6.)*UAVG[zdf+zdc-4+k-2]-tmp;
			dUdz[zdf+zdc-4+k] = 0.5*(UAVG[zdf+zdc-4+k+1]-UAVG[zdf+zdc-4+k-1]);
		}
		else
		{
			dUdz[zdf+zdc-4+k] = (2./3.)*(UAVG[zdf+zdc-4+k+1]-UAVG[zdf+zdc-4+k-1])-(1./12.)*(UAVG[zdf+zdc-4+k+2]-UAVG[zdf+zdc-4+k-2]);
		}
	}
	
	extx=xdf-1;
	exty=ydf;
	xd=xdf;yd=ydf;zd=zdf;
	V.vel=V.velfl;
	V.den=V.denfl;
	
	for (i=0;i<extx;i++)
	{
		for (j=0;j<exty;j++)
		{
			for (k=1;k<zd;k++)
			{
				
				u=Fb(V.vel,i,j,k,0,xd,yd,zd,3);
				v=Fb(V.vel,i,j,k,1,xd,yd,zd,3);
				w=Fb(V.vel,i,j,k,2,xd,yd,zd,3);

				dens = Fb(V.den,i,j,k,0,xd,yd,zd,1);
				dfl = (1./3.)*dens-DAVG[k];

				*(uu+k) += u*u;
				*(uv+k) += u*v;
				*(uw+k) += u*w;
				*(vv+k) += v*v;
				*(vw+k) += v*w;
				*(ww+k) += w*w;
				
				*(SKEWu+k) += u*u*u;
				*(SKEWv+k) += v*v*v;
				*(SKEWw+k) += w*w*w;
				*(SKEWp+k) += dfl*dfl*dfl;
				*(SKEWuw+k) += u*u*u*w*w*w;
				*(SKEWuv+k) += u*u*u*v*v*v;
				*(SKEWvw+k) += v*v*v*w*w*w;
				
				*(KURTu+k) += u*u*u*u;
				*(KURTv+k) += v*v*v*v;
				*(KURTw+k) += w*w*w*w;
				*(KURTp+k) += dfl*dfl*dfl*dfl;
				*(KURTuw+k) += u*u*u*u*w*w*w*w;
				*(KURTvw+k) += v*v*v*v*w*w*w*w;
				*(KURTuv+k) += u*u*u*u*v*v*v*v;
	
				*(dd+k) += dfl*dfl;
				
				*(uuw+k) += u*u*w;
				*(vvw+k) += v*v*w;
				*(www+k) += w*w*w;
				*(uww+k) += u*w*w;

				*(pu+k) += u*dfl;
				*(pv+k) += v*dfl;
				*(pw+k) += w*dfl;
			
				if (j==0)
				{
					jp=yd-1;
					jpp=yd-2;
					jn=1;
					jnn=2;
				}
				else if (j==1)
				{
					jp=0;
					jpp=yd-1;
					jn=2;
					jnn=3;
				}
				else if (j==(yd-2))
				{
					jp=j-1;
					jpp=j-2;
					jn=j+1;
					jnn=0;
				}
				else if (j==(yd-1))
				{
					jp=j-1;
					jpp=j-2;
					jn=0;
					jnn=1;
				}
				else
				{
					jp=j-1;
					jpp=j-2;
					jn=j+1;
					jnn=j+2;
				}
				ip=i-1;
				ipp=i-2;
				in=i+1;
				inn=i+2;
				
				dudx=((2./3.)*(Fb(V.vel,in,j,k,0,xd,yd,zd,3)-Fb(V.vel,ip,j,k,0,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,inn,j,k,0,xd,yd,zd,3)-Fb(V.vel,ipp,j,k,0,xd,yd,zd,3)));	
				dvdx=((2./3.)*(Fb(V.vel,in,j,k,1,xd,yd,zd,3)-Fb(V.vel,ip,j,k,1,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,inn,j,k,1,xd,yd,zd,3)-Fb(V.vel,ipp,j,k,1,xd,yd,zd,3)));	
				dwdx=((2./3.)*(Fb(V.vel,in,j,k,2,xd,yd,zd,3)-Fb(V.vel,ip,j,k,2,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,inn,j,k,2,xd,yd,zd,3)-Fb(V.vel,ipp,j,k,2,xd,yd,zd,3)));
				
				dudy=((2./3.)*(Fb(V.vel,i,jn,k,0,xd,yd,zd,3)-Fb(V.vel,i,jp,k,0,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,jnn,k,0,xd,yd,zd,3)-Fb(V.vel,i,jpp,k,0,xd,yd,zd,3)));	
				dvdy=((2./3.)*(Fb(V.vel,i,jn,k,1,xd,yd,zd,3)-Fb(V.vel,i,jp,k,1,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,jnn,k,1,xd,yd,zd,3)-Fb(V.vel,i,jpp,k,1,xd,yd,zd,3)));	
				dwdy=((2./3.)*(Fb(V.vel,i,jn,k,2,xd,yd,zd,3)-Fb(V.vel,i,jp,k,2,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,jnn,k,2,xd,yd,zd,3)-Fb(V.vel,i,jpp,k,2,xd,yd,zd,3)));
				
				if (k==1)	
				{
					tmp=(2./3.)*(Fb(V.vel,i,j,(k+3),2,xd,yd,zd,3)-Fb(V.vel,i,j,(k+1),2,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+4),2,xd,yd,zd,3)-Fb(V.vel,i,j,k,2,xd,yd,zd,3));
					dwdz=-(7./3.)*Fb(V.vel,i,j,k,2,xd,yd,zd,3)+6.*Fb(V.vel,i,j,(k+1),2,xd,yd,zd,3)-3.*Fb(V.vel,i,j,(k+2),2,xd,yd,zd,3)-(2./3.)*Fb(V.vel,i,j,(k+3),2,xd,yd,zd,3)+3.*tmp;
					tmp=(2./3.)*(Fb(V.vel,i,j,(k+3),1,xd,yd,zd,3)-Fb(V.vel,i,j,(k+1),1,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+4),1,xd,yd,zd,3)-Fb(V.vel,i,j,k,1,xd,yd,zd,3));
					dvdz=-(7./3.)*Fb(V.vel,i,j,k,1,xd,yd,zd,3)+6.*Fb(V.vel,i,j,(k+1),1,xd,yd,zd,3)-3.*Fb(V.vel,i,j,(k+2),1,xd,yd,zd,3)-(2./3.)*Fb(V.vel,i,j,(k+3),1,xd,yd,zd,3)+3.*tmp;
					tmp=(2./3.)*(Fb(V.vel,i,j,(k+3),0,xd,yd,zd,3)-Fb(V.vel,i,j,(k+1),0,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+4),0,xd,yd,zd,3)-Fb(V.vel,i,j,k,0,xd,yd,zd,3));
					dudz=-(7./3.)*Fb(V.vel,i,j,k,0,xd,yd,zd,3)+6.*Fb(V.vel,i,j,(k+1),0,xd,yd,zd,3)-3.*Fb(V.vel,i,j,(k+2),0,xd,yd,zd,3)-(2./3.)*Fb(V.vel,i,j,(k+3),0,xd,yd,zd,3)+3.*tmp;
				}
				else if (k==(zd-1))
				{
// 					tmp=(2./3.)*(Fb(V.vel,i,j,(k-1),2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-3),2,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,k,2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-4),2,xd,yd,zd,3));
// 					dwdz=(7./3.)*Fb(V.vel,i,j,k,2,xd,yd,zd,3)-6.*Fb(V.vel,i,j,(k-1),2,xd,yd,zd,3)+3.*Fb(V.vel,i,j,(k-2),2,xd,yd,zd,3)+(2./3.)*Fb(V.vel,i,j,(k-3),2,xd,yd,zd,3)+3.*tmp;
// 					tmp=(2./3.)*(Fb(V.vel,i,j,(k-1),1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-3),1,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,k,1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-4),1,xd,yd,zd,3));
// 					dvdz=(7./3.)*Fb(V.vel,i,j,k,1,xd,yd,zd,3)-6.*Fb(V.vel,i,j,(k-1),1,xd,yd,zd,3)+3.*Fb(V.vel,i,j,(k-2),1,xd,yd,zd,3)+(2./3.)*Fb(V.vel,i,j,(k-3),1,xd,yd,zd,3)+3.*tmp;
// 					tmp=(2./3.)*(Fb(V.vel,i,j,(k-1),0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-3),0,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,k,0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-4),0,xd,yd,zd,3));
// 					dudz=(7./3.)*Fb(V.vel,i,j,k,0,xd,yd,zd,3)-6.*Fb(V.vel,i,j,(k-1),0,xd,yd,zd,3)+3.*Fb(V.vel,i,j,(k-2),0,xd,yd,zd,3)+(2./3.)*Fb(V.vel,i,j,(k-3),0,xd,yd,zd,3)+3.*tmp;
					dudz = (coef*(Fb(V.velc,(i/GR),(j/GR),2,0,xdc,ydc,zdc,3)-Fb(V.vel,i,j,k,0,xd,yd,zd,3))+GR*(Fb(V.vel,i,j,k,0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),0,xd,yd,zd,3)))*1./(1.+GR);
					dvdz = (coef*(Fb(V.velc,(i/GR),(j/GR),2,1,xdc,ydc,zdc,3)-Fb(V.vel,i,j,k,1,xd,yd,zd,3))+GR*(Fb(V.vel,i,j,k,1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),1,xd,yd,zd,3)))*1./(1.+GR);
					dwdz = (coef*(Fb(V.velc,(i/GR),(j/GR),2,2,xdc,ydc,zdc,3)-Fb(V.vel,i,j,k,2,xd,yd,zd,3))+GR*(Fb(V.vel,i,j,k,2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),2,xd,yd,zd,3)))*1./(1.+GR);
				}
				else if (k==2)
				{
// 					tmp=(2./3.)*(Fb(V.vel,i,j,(k+2),2,xd,yd,zd,3)-Fb(V.vel,i,j,k,2,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+3),2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),2,xd,yd,zd,3));
// 					dwdz=-(1./6.)*Fb(V.vel,i,j,(k-1),2,xd,yd,zd,3)-(3./2.)*Fb(V.vel,i,j,k,2,xd,yd,zd,3)+(3./2.)*Fb(V.vel,i,j,(k+1),2,xd,yd,zd,3)+(1./6.)*Fb(V.vel,i,j,(k+2),2,xd,yd,zd,3)-tmp;
// 					tmp=(2./3.)*(Fb(V.vel,i,j,(k+2),1,xd,yd,zd,3)-Fb(V.vel,i,j,k,1,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+3),1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),1,xd,yd,zd,3));
// 					dvdz=-(1./6.)*Fb(V.vel,i,j,(k-1),1,xd,yd,zd,3)-(3./2.)*Fb(V.vel,i,j,k,1,xd,yd,zd,3)+(3./2.)*Fb(V.vel,i,j,(k+1),1,xd,yd,zd,3)+(1./6.)*Fb(V.vel,i,j,(k+2),1,xd,yd,zd,3)-tmp;
// 					tmp=(2./3.)*(Fb(V.vel,i,j,(k+2),0,xd,yd,zd,3)-Fb(V.vel,i,j,k,0,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+3),0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),0,xd,yd,zd,3));
// 					dudz=-(1./6.)*Fb(V.vel,i,j,(k-1),0,xd,yd,zd,3)-(3./2.)*Fb(V.vel,i,j,k,0,xd,yd,zd,3)+(3./2.)*Fb(V.vel,i,j,(k+1),0,xd,yd,zd,3)+(1./6.)*Fb(V.vel,i,j,(k+2),0,xd,yd,zd,3)-tmp;
					dudz = 0.5*(Fb(V.vel,i,j,(k+1),0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),0,xd,yd,zd,3));
					dvdz = 0.5*(Fb(V.vel,i,j,(k+1),1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),1,xd,yd,zd,3));
					dwdz = 0.5*(Fb(V.vel,i,j,(k+1),2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),2,xd,yd,zd,3));
				}
				else if (k==(zd-2))
				{
// 					tmp=(2./3.)*(Fb(V.vel,i,j,k,2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-2),2,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+1),2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-3),2,xd,yd,zd,3));
// 					dwdz=(1./6.)*Fb(V.vel,i,j,(k+1),2,xd,yd,zd,3)+(3./2.)*Fb(V.vel,i,j,k,2,xd,yd,zd,3)-(3./2.)*Fb(V.vel,i,j,(k-1),2,xd,yd,zd,3)-(1./6.)*Fb(V.vel,i,j,(k-2),2,xd,yd,zd,3)-tmp;
// 					tmp=(2./3.)*(Fb(V.vel,i,j,k,1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-2),1,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+1),1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-3),1,xd,yd,zd,3));
// 					dvdz=(1./6.)*Fb(V.vel,i,j,(k+1),1,xd,yd,zd,3)+(3./2.)*Fb(V.vel,i,j,k,1,xd,yd,zd,3)-(3./2.)*Fb(V.vel,i,j,(k-1),1,xd,yd,zd,3)-(1./6.)*Fb(V.vel,i,j,(k-2),1,xd,yd,zd,3)-tmp;
// 					tmp=(2./3.)*(Fb(V.vel,i,j,k,0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-2),0,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+1),0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-3),0,xd,yd,zd,3));
// 					dudz=(1./6.)*Fb(V.vel,i,j,(k+1),0,xd,yd,zd,3)+(3./2.)*Fb(V.vel,i,j,k,0,xd,yd,zd,3)-(3./2.)*Fb(V.vel,i,j,(k-1),0,xd,yd,zd,3)-(1./6.)*Fb(V.vel,i,j,(k-2),0,xd,yd,zd,3)-tmp;
					dudz = 0.5*(Fb(V.vel,i,j,(k+1),0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),0,xd,yd,zd,3));
					dvdz = 0.5*(Fb(V.vel,i,j,(k+1),1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),1,xd,yd,zd,3));
					dwdz = 0.5*(Fb(V.vel,i,j,(k+1),2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),2,xd,yd,zd,3));
				}
				else
				{
//					dwdz=(2./3.)*(Fb(V.vel,i,j,(k+1),2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),2,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+2),2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-2),2,xd,yd,zd,3));
//					dvdz=(2./3.)*(Fb(V.vel,i,j,(k+1),1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),1,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+2),1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-2),1,xd,yd,zd,3));
//					dudz=(2./3.)*(Fb(V.vel,i,j,(k+1),0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),0,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+2),0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-2),0,xd,yd,zd,3));
					dudz = 0.5*(Fb(V.vel,i,j,(k+1),0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),0,xd,yd,zd,3));
                                        dvdz = 0.5*(Fb(V.vel,i,j,(k+1),1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),1,xd,yd,zd,3));
                                        dwdz = 0.5*(Fb(V.vel,i,j,(k+1),2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),2,xd,yd,zd,3));
				}
				
				ox=(dwdy-dvdz);
				oy=(dudz-dwdx);
				oz=(dvdx-dudy);

				*(wx+k) += ox*ox;
             			*(wy+k) += oy*oy;
           			*(wz+k) += oz*oz;		
				
				*(Paiuu+k) += dfl*dudx;
				*(Paivv+k) += dfl*dvdy;
				*(Paiww+k) += dfl*dwdz;
				*(Paiuw+k) += dfl*(dudz+dwdx);
				
				*(usuw+k) += 0.5*u*(dudz+dwdx);
				*(vsvw+k) += 0.5*v*(dvdz+dwdy);
				*(wsww+k) += 0.5*w*(dwdz+dwdz);
				*(wsuw+k) += 0.5*w*(dudz+dwdx);
				*(usww+k) += 0.5*u*(dwdz+dwdz);
				
				*(Epuu+k) += dudx*(dudx+dudx)+dudy*(dudy+dvdx)+dudz*(dudz+dwdx);
				*(Epvv+k) += dvdx*(dvdx+dudy)+dvdy*(dvdy+dvdy)+dvdz*(dvdz+dwdy);
				*(Epww+k) += dwdx*(dwdx+dudz)+dwdy*(dwdy+dvdz)+dwdz*(dwdz+dwdz);
				*(Epuw+k) += dudx*(dwdx+dudz)+dudy*(dwdy+dvdz)+dudz*(dwdz+dwdz)+dwdx*(dudx+dudx)+dwdy*(dudy+dvdx)+dwdz*(dudz+dwdx);
			}
		}
	}
	
	extx=xdc-1;
	exty=ydc;
	xd=xdc;yd=ydc;zd=zdc;
	V.vel=V.velc;
	V.den=V.denc;
	
	for (i=0;i<extx;i++)
	{
		for (j=0;j<exty;j++)
		{
			for (k=2;k<(zd-2);k++)
			{
				
				u=Fb(V.vel,i,j,k,0,xd,yd,zd,3);
				v=Fb(V.vel,i,j,k,1,xd,yd,zd,3);
				w=Fb(V.vel,i,j,k,2,xd,yd,zd,3);

				dens = Fb(V.den,i,j,k,0,xd,yd,zd,1);
				dfl = (1./3.)*dens-DAVG[zdf-2+k];
				
				*(uu+zdf-2+k) += u*u;
				*(uv+zdf-2+k) += u*v;
				*(uw+zdf-2+k) += u*w;
				*(vv+zdf-2+k) += v*v;
				*(vw+zdf-2+k) += v*w;
				*(ww+zdf-2+k) += w*w;
				
				*(SKEWu+zdf-2+k) += u*u*u;
				*(SKEWv+zdf-2+k) += v*v*v;
				*(SKEWw+zdf-2+k) += w*w*w;
				*(SKEWp+zdf-2+k) += dfl*dfl*dfl;
				*(SKEWuw+zdf-2+k) += u*u*u*w*w*w;
				*(SKEWuv+zdf-2+k) += u*u*u*v*v*v;
				*(SKEWvw+zdf-2+k) += v*v*v*w*w*w;
				
				*(KURTu+zdf-2+k) += u*u*u*u;
				*(KURTv+zdf-2+k) += v*v*v*v;
				*(KURTw+zdf-2+k) += w*w*w*w;
				*(KURTp+zdf-2+k) += dfl*dfl*dfl*dfl;
				*(KURTuw+zdf-2+k) += u*u*u*u*w*w*w*w;
				*(KURTvw+zdf-2+k) += v*v*v*v*w*w*w*w;
				*(KURTuv+zdf-2+k) += u*u*u*u*v*v*v*v;
	
				*(dd+zdf-2+k) += dfl*dfl;
				
				*(uuw+zdf-2+k) += u*u*w;
				*(vvw+zdf-2+k) += v*v*w;
				*(www+zdf-2+k) += w*w*w;
				*(uww+zdf-2+k) += u*w*w;

				*(pu+zdf-2+k) += u*dfl;
				*(pv+zdf-2+k) += v*dfl;
				*(pw+zdf-2+k) += w*dfl;
			
				if (j==0)
				{
					jp=yd-1;
					jpp=yd-2;
					jn=1;
					jnn=2;
				}
				else if (j==1)
				{
					jp=0;
					jpp=yd-1;
					jn=2;
					jnn=3;
				}
				else if (j==(yd-2))
				{
					jp=j-1;
					jpp=j-2;
					jn=j+1;
					jnn=0;
				}
				else if (j==(yd-1))
				{
					jp=j-1;
					jpp=j-2;
					jn=0;
					jnn=1;
				}
				else
				{
					jp=j-1;
					jpp=j-2;
					jn=j+1;
					jnn=j+2;
				}
				ip=i-1;
				ipp=i-2;
				in=i+1;
				inn=i+2;
				
				dudx=coef*((2./3.)*(Fb(V.vel,in,j,k,0,xd,yd,zd,3)-Fb(V.vel,ip,j,k,0,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,inn,j,k,0,xd,yd,zd,3)-Fb(V.vel,ipp,j,k,0,xd,yd,zd,3)));	
				dvdx=coef*((2./3.)*(Fb(V.vel,in,j,k,1,xd,yd,zd,3)-Fb(V.vel,ip,j,k,1,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,inn,j,k,1,xd,yd,zd,3)-Fb(V.vel,ipp,j,k,1,xd,yd,zd,3)));	
				dwdx=coef*((2./3.)*(Fb(V.vel,in,j,k,2,xd,yd,zd,3)-Fb(V.vel,ip,j,k,2,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,inn,j,k,2,xd,yd,zd,3)-Fb(V.vel,ipp,j,k,2,xd,yd,zd,3)));
				
				dudy=coef*((2./3.)*(Fb(V.vel,i,jn,k,0,xd,yd,zd,3)-Fb(V.vel,i,jp,k,0,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,jnn,k,0,xd,yd,zd,3)-Fb(V.vel,i,jpp,k,0,xd,yd,zd,3)));	
				dvdy=coef*((2./3.)*(Fb(V.vel,i,jn,k,1,xd,yd,zd,3)-Fb(V.vel,i,jp,k,1,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,jnn,k,1,xd,yd,zd,3)-Fb(V.vel,i,jpp,k,1,xd,yd,zd,3)));	
				dwdy=coef*((2./3.)*(Fb(V.vel,i,jn,k,2,xd,yd,zd,3)-Fb(V.vel,i,jp,k,2,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,jnn,k,2,xd,yd,zd,3)-Fb(V.vel,i,jpp,k,2,xd,yd,zd,3)));
				
				if (k==2)
				{
// 					tmp=(2./3.)*(Fb(V.vel,i,j,(k+3),2,xd,yd,zd,3)-Fb(V.vel,i,j,(k+1),2,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+4),2,xd,yd,zd,3)-Fb(V.vel,i,j,k,2,xd,yd,zd,3));
// 					dwdz=coef*(-(7./3.)*Fb(V.vel,i,j,k,2,xd,yd,zd,3)+6.*Fb(V.vel,i,j,(k+1),2,xd,yd,zd,3)-3.*Fb(V.vel,i,j,(k+2),2,xd,yd,zd,3)-(2./3.)*Fb(V.vel,i,j,(k+3),2,xd,yd,zd,3)+3.*tmp);
// 					tmp=(2./3.)*(Fb(V.vel,i,j,(k+3),1,xd,yd,zd,3)-Fb(V.vel,i,j,(k+1),1,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+4),1,xd,yd,zd,3)-Fb(V.vel,i,j,k,1,xd,yd,zd,3));
// 					dvdz=coef*(-(7./3.)*Fb(V.vel,i,j,k,1,xd,yd,zd,3)+6.*Fb(V.vel,i,j,(k+1),1,xd,yd,zd,3)-3.*Fb(V.vel,i,j,(k+2),1,xd,yd,zd,3)-(2./3.)*Fb(V.vel,i,j,(k+3),1,xd,yd,zd,3)+3.*tmp);
// 					tmp=(2./3.)*(Fb(V.vel,i,j,(k+3),0,xd,yd,zd,3)-Fb(V.vel,i,j,(k+1),0,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+4),0,xd,yd,zd,3)-Fb(V.vel,i,j,k,0,xd,yd,zd,3));
// 					dudz=coef*(-(7./3.)*Fb(V.vel,i,j,k,0,xd,yd,zd,3)+6.*Fb(V.vel,i,j,(k+1),0,xd,yd,zd,3)-3.*Fb(V.vel,i,j,(k+2),0,xd,yd,zd,3)-(2./3.)*Fb(V.vel,i,j,(k+3),0,xd,yd,zd,3)+3.*tmp);
					dudz = 0.5*coef*(Fb(V.vel,i,j,(k+1),0,xd,yd,zd,3)-Fb(V.velfl,(i*GR),(j*GR),(zdf-1),0,xdf,ydf,zdf,3));
					dvdz = 0.5*coef*(Fb(V.vel,i,j,(k+1),1,xd,yd,zd,3)-Fb(V.velfl,(i*GR),(j*GR),(zdf-1),1,xdf,ydf,zdf,3));
					dwdz = 0.5*coef*(Fb(V.vel,i,j,(k+1),2,xd,yd,zd,3)-Fb(V.velfl,(i*GR),(j*GR),(zdf-1),2,xdf,ydf,zdf,3));
				}
				else if (k==3)
				{
// 					tmp=(2./3.)*(Fb(V.vel,i,j,(k+2),2,xd,yd,zd,3)-Fb(V.vel,i,j,k,2,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+3),2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),2,xd,yd,zd,3));
// 					dwdz=coef*(-(1./6.)*Fb(V.vel,i,j,(k-1),2,xd,yd,zd,3)-(3./2.)*Fb(V.vel,i,j,k,2,xd,yd,zd,3)+(3./2.)*Fb(V.vel,i,j,(k+1),2,xd,yd,zd,3)+(1./6.)*Fb(V.vel,i,j,(k+2),2,xd,yd,zd,3)-tmp);
// 					tmp=(2./3.)*(Fb(V.vel,i,j,(k+2),1,xd,yd,zd,3)-Fb(V.vel,i,j,k,1,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+3),1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),1,xd,yd,zd,3));
// 					dvdz=coef*(-(1./6.)*Fb(V.vel,i,j,(k-1),1,xd,yd,zd,3)-(3./2.)*Fb(V.vel,i,j,k,1,xd,yd,zd,3)+(3./2.)*Fb(V.vel,i,j,(k+1),1,xd,yd,zd,3)+(1./6.)*Fb(V.vel,i,j,(k+2),1,xd,yd,zd,3)-tmp);
// 					tmp=(2./3.)*(Fb(V.vel,i,j,(k+2),0,xd,yd,zd,3)-Fb(V.vel,i,j,k,0,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+3),0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),0,xd,yd,zd,3));
// 					dudz=coef*(-(1./6.)*Fb(V.vel,i,j,(k-1),0,xd,yd,zd,3)-(3./2.)*Fb(V.vel,i,j,k,0,xd,yd,zd,3)+(3./2.)*Fb(V.vel,i,j,(k+1),0,xd,yd,zd,3)+(1./6.)*Fb(V.vel,i,j,(k+2),0,xd,yd,zd,3)-tmp);
					dudz = 0.5*coef*(Fb(V.vel,i,j,(k+1),0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),0,xd,yd,zd,3));
					dvdz = 0.5*coef*(Fb(V.vel,i,j,(k+1),1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),1,xd,yd,zd,3));
					dwdz = 0.5*coef*(Fb(V.vel,i,j,(k+1),2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),2,xd,yd,zd,3));
				}
				else if (k==(zdc-4))
				{
// 					tmp=(2./3.)*(Fb(V.vel,i,j,k,2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-2),2,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+1),2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-3),2,xd,yd,zd,3));
// 					dwdz=coef*((1./6.)*Fb(V.vel,i,j,(k+1),2,xd,yd,zd,3)+(3./2.)*Fb(V.vel,i,j,k,2,xd,yd,zd,3)-(3./2.)*Fb(V.vel,i,j,(k-1),2,xd,yd,zd,3)-(1./6.)*Fb(V.vel,i,j,(k-2),2,xd,yd,zd,3)-tmp);
// 					tmp=(2./3.)*(Fb(V.vel,i,j,k,1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-2),1,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+1),1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-3),1,xd,yd,zd,3));
// 					dvdz=coef*((1./6.)*Fb(V.vel,i,j,(k+1),1,xd,yd,zd,3)+(3./2.)*Fb(V.vel,i,j,k,1,xd,yd,zd,3)-(3./2.)*Fb(V.vel,i,j,(k-1),1,xd,yd,zd,3)-(1./6.)*Fb(V.vel,i,j,(k-2),1,xd,yd,zd,3)-tmp);
// 					tmp=(2./3.)*(Fb(V.vel,i,j,k,0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-2),0,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+1),0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-3),0,xd,yd,zd,3));
// 					dudz=coef*((1./6.)*Fb(V.vel,i,j,(k+1),0,xd,yd,zd,3)+(3./2.)*Fb(V.vel,i,j,k,0,xd,yd,zd,3)-(3./2.)*Fb(V.vel,i,j,(k-1),0,xd,yd,zd,3)-(1./6.)*Fb(V.vel,i,j,(k-2),0,xd,yd,zd,3)-tmp);
					dudz = 0.5*coef*(Fb(V.vel,i,j,(k+1),0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),0,xd,yd,zd,3));
					dvdz = 0.5*coef*(Fb(V.vel,i,j,(k+1),1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),1,xd,yd,zd,3));
					dwdz = 0.5*coef*(Fb(V.vel,i,j,(k+1),2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),2,xd,yd,zd,3));
				}
				else if (k==(zdc-3))
				{
// 					tmp=(2./3.)*(Fb(V.vel,i,j,(k-1),2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-3),2,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,k,2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-4),2,xd,yd,zd,3));
// 					dwdz=coef*((7./3.)*Fb(V.vel,i,j,k,2,xd,yd,zd,3)-6.*Fb(V.vel,i,j,(k-1),2,xd,yd,zd,3)+3.*Fb(V.vel,i,j,(k-2),2,xd,yd,zd,3)+(2./3.)*Fb(V.vel,i,j,(k-3),2,xd,yd,zd,3)+3.*tmp);
// 					tmp=(2./3.)*(Fb(V.vel,i,j,(k-1),1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-3),1,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,k,1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-4),1,xd,yd,zd,3));
// 					dvdz=coef*((7./3.)*Fb(V.vel,i,j,k,1,xd,yd,zd,3)-6.*Fb(V.vel,i,j,(k-1),1,xd,yd,zd,3)+3.*Fb(V.vel,i,j,(k-2),1,xd,yd,zd,3)+(2./3.)*Fb(V.vel,i,j,(k-3),1,xd,yd,zd,3)+3.*tmp);
// 					tmp=(2./3.)*(Fb(V.vel,i,j,(k-1),0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-3),0,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,k,0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-4),0,xd,yd,zd,3));
// 					dudz=coef*((7./3.)*Fb(V.vel,i,j,k,0,xd,yd,zd,3)-6.*Fb(V.vel,i,j,(k-1),0,xd,yd,zd,3)+3.*Fb(V.vel,i,j,(k-2),0,xd,yd,zd,3)+(2./3.)*Fb(V.vel,i,j,(k-3),0,xd,yd,zd,3)+3.*tmp);
					dudz = 0.5*coef*(Fb(V.velfu,(i*GR),(j*GR),0,0,xdf,ydf,zdf,3)-Fb(V.vel,i,j,(k-1),0,xd,yd,zd,3));
					dvdz = 0.5*coef*(Fb(V.velfu,(i*GR),(j*GR),0,1,xdf,ydf,zdf,3)-Fb(V.vel,i,j,(k-1),1,xd,yd,zd,3));
					dwdz = 0.5*coef*(Fb(V.velfu,(i*GR),(j*GR),0,2,xdf,ydf,zdf,3)-Fb(V.vel,i,j,(k-1),2,xd,yd,zd,3));
				}
				else
				{
//					dwdz=coef*((2./3.)*(Fb(V.vel,i,j,(k+1),2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),2,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+2),2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-2),2,xd,yd,zd,3)));
//					dvdz=coef*((2./3.)*(Fb(V.vel,i,j,(k+1),1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),1,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+2),1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-2),1,xd,yd,zd,3)));
//					dudz=coef*((2./3.)*(Fb(V.vel,i,j,(k+1),0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),0,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+2),0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-2),0,xd,yd,zd,3)));
					dudz = 0.5*coef*(Fb(V.vel,i,j,(k+1),0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),0,xd,yd,zd,3));
                                        dvdz = 0.5*coef*(Fb(V.vel,i,j,(k+1),1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),1,xd,yd,zd,3));
                                        dwdz = 0.5*coef*(Fb(V.vel,i,j,(k+1),2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),2,xd,yd,zd,3));
				}
				
				ox=(dwdy-dvdz);
				oy=(dudz-dwdx);
				oz=(dvdx-dudy);

				*(wx+zdf-2+k) += ox*ox;
             			*(wy+zdf-2+k) += oy*oy;
           			*(wz+zdf-2+k) += oz*oz;	
				
				*(Paiuu+zdf-2+k) += dfl*dudx;
				*(Paivv+zdf-2+k) += dfl*dvdy;
				*(Paiww+zdf-2+k) += dfl*dwdz;
				*(Paiuw+zdf-2+k) += dfl*(dudz+dwdx);
				
				*(usuw+zdf-2+k) += 0.5*u*(dudz+dwdx);
				*(vsvw+zdf-2+k) += 0.5*v*(dvdz+dwdy);
				*(wsww+zdf-2+k) += 0.5*w*(dwdz+dwdz);
				*(wsuw+zdf-2+k) += 0.5*w*(dudz+dwdx);
				*(usww+zdf-2+k) += 0.5*u*(dwdz+dwdz);
				
				*(Epuu+zdf-2+k) += dudx*(dudx+dudx)+dudy*(dudy+dvdx)+dudz*(dudz+dwdx);
				*(Epvv+zdf-2+k) += dvdx*(dvdx+dudy)+dvdy*(dvdy+dvdy)+dvdz*(dvdz+dwdy);
				*(Epww+zdf-2+k) += dwdx*(dwdx+dudz)+dwdy*(dwdy+dvdz)+dwdz*(dwdz+dwdz);
				*(Epuw+zdf-2+k) += dudx*(dwdx+dudz)+dudy*(dwdy+dvdz)+dudz*(dwdz+dwdz)+dwdx*(dudx+dudx)+dwdy*(dudy+dvdx)+dwdz*(dudz+dwdx);
			}
		}
	}
	
	extx=xdf-1;
	exty=ydf;
	xd=xdf;yd=ydf;zd=zdf;
	V.vel=V.velfu;
	V.den=V.denfu;
	
	for (i=0;i<extx;i++)
	{
		for (j=0;j<exty;j++)
		{
			for (k=0;k<(zd-1);k++)
			{
				
				u=Fb(V.vel,i,j,k,0,xd,yd,zd,3);
				v=Fb(V.vel,i,j,k,1,xd,yd,zd,3);
				w=Fb(V.vel,i,j,k,2,xd,yd,zd,3);

				dens = Fb(V.den,i,j,k,0,xd,yd,zd,1);
				dfl = (1./3.)*dens-DAVG[zdf+zdc-4+k];

				*(uu+zdf+zdc-4+k) += u*u;
				*(uv+zdf+zdc-4+k) += u*v;
				*(uw+zdf+zdc-4+k) += u*w;
				*(vv+zdf+zdc-4+k) += v*v;
				*(vw+zdf+zdc-4+k) += v*w;
				*(ww+zdf+zdc-4+k) += w*w;
				
				*(SKEWu+zdf+zdc-4+k) += u*u*u;
				*(SKEWv+zdf+zdc-4+k) += v*v*v;
				*(SKEWw+zdf+zdc-4+k) += w*w*w;
				*(SKEWp+zdf+zdc-4+k) += dfl*dfl*dfl;
				*(SKEWuw+zdf+zdc-4+k) += u*u*u*w*w*w;
				*(SKEWuv+zdf+zdc-4+k) += u*u*u*v*v*v;
				*(SKEWvw+zdf+zdc-4+k) += v*v*v*w*w*w;
				
				*(KURTu+zdf+zdc-4+k) += u*u*u*u;
				*(KURTv+zdf+zdc-4+k) += v*v*v*v;
				*(KURTw+zdf+zdc-4+k) += w*w*w*w;
				*(KURTp+zdf+zdc-4+k) += dfl*dfl*dfl*dfl;
				*(KURTuw+zdf+zdc-4+k) += u*u*u*u*w*w*w*w;
				*(KURTvw+zdf+zdc-4+k) += v*v*v*v*w*w*w*w;
				*(KURTuv+zdf+zdc-4+k) += u*u*u*u*v*v*v*v;
	
				*(dd+zdf+zdc-4+k) += dfl*dfl;
				
				*(uuw+zdf+zdc-4+k) += u*u*w;
				*(vvw+zdf+zdc-4+k) += v*v*w;
				*(www+zdf+zdc-4+k) += w*w*w;
				*(uww+zdf+zdc-4+k) += u*w*w;

				*(pu+zdf+zdc-4+k) += u*dfl;
				*(pv+zdf+zdc-4+k) += v*dfl;
				*(pw+zdf+zdc-4+k) += w*dfl;
			
				if (j==0)
				{
					jp=yd-1;
					jpp=yd-2;
					jn=1;
					jnn=2;
				}
				else if (j==1)
				{
					jp=0;
					jpp=yd-1;
					jn=2;
					jnn=3;
				}
				else if (j==(yd-2))
				{
					jp=j-1;
					jpp=j-2;
					jn=j+1;
					jnn=0;
				}
				else if (j==(yd-1))
				{
					jp=j-1;
					jpp=j-2;
					jn=0;
					jnn=1;
				}
				else
				{
					jp=j-1;
					jpp=j-2;
					jn=j+1;
					jnn=j+2;
				}
				ip=i-1;
				ipp=i-2;
				in=i+1;
				inn=i+2;
				
				dudx=((2./3.)*(Fb(V.vel,in,j,k,0,xd,yd,zd,3)-Fb(V.vel,ip,j,k,0,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,inn,j,k,0,xd,yd,zd,3)-Fb(V.vel,ipp,j,k,0,xd,yd,zd,3)));	
				dvdx=((2./3.)*(Fb(V.vel,in,j,k,1,xd,yd,zd,3)-Fb(V.vel,ip,j,k,1,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,inn,j,k,1,xd,yd,zd,3)-Fb(V.vel,ipp,j,k,1,xd,yd,zd,3)));	
				dwdx=((2./3.)*(Fb(V.vel,in,j,k,2,xd,yd,zd,3)-Fb(V.vel,ip,j,k,2,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,inn,j,k,2,xd,yd,zd,3)-Fb(V.vel,ipp,j,k,2,xd,yd,zd,3)));
				
				dudy=((2./3.)*(Fb(V.vel,i,jn,k,0,xd,yd,zd,3)-Fb(V.vel,i,jp,k,0,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,jnn,k,0,xd,yd,zd,3)-Fb(V.vel,i,jpp,k,0,xd,yd,zd,3)));	
				dvdy=((2./3.)*(Fb(V.vel,i,jn,k,1,xd,yd,zd,3)-Fb(V.vel,i,jp,k,1,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,jnn,k,1,xd,yd,zd,3)-Fb(V.vel,i,jpp,k,1,xd,yd,zd,3)));	
				dwdy=((2./3.)*(Fb(V.vel,i,jn,k,2,xd,yd,zd,3)-Fb(V.vel,i,jp,k,2,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,jnn,k,2,xd,yd,zd,3)-Fb(V.vel,i,jpp,k,2,xd,yd,zd,3)));
				
				if (k==0)	
				{
// 					tmp=(2./3.)*(Fb(V.vel,i,j,(k+3),2,xd,yd,zd,3)-Fb(V.vel,i,j,(k+1),2,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+4),2,xd,yd,zd,3)-Fb(V.vel,i,j,k,2,xd,yd,zd,3));
// 					dwdz=-(7./3.)*Fb(V.vel,i,j,k,2,xd,yd,zd,3)+6.*Fb(V.vel,i,j,(k+1),2,xd,yd,zd,3)-3.*Fb(V.vel,i,j,(k+2),2,xd,yd,zd,3)-(2./3.)*Fb(V.vel,i,j,(k+3),2,xd,yd,zd,3)+3.*tmp;
// 					tmp=(2./3.)*(Fb(V.vel,i,j,(k+3),1,xd,yd,zd,3)-Fb(V.vel,i,j,(k+1),1,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+4),1,xd,yd,zd,3)-Fb(V.vel,i,j,k,1,xd,yd,zd,3));
// 					dvdz=-(7./3.)*Fb(V.vel,i,j,k,1,xd,yd,zd,3)+6.*Fb(V.vel,i,j,(k+1),1,xd,yd,zd,3)-3.*Fb(V.vel,i,j,(k+2),1,xd,yd,zd,3)-(2./3.)*Fb(V.vel,i,j,(k+3),1,xd,yd,zd,3)+3.*tmp;
// 					tmp=(2./3.)*(Fb(V.vel,i,j,(k+3),0,xd,yd,zd,3)-Fb(V.vel,i,j,(k+1),0,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+4),0,xd,yd,zd,3)-Fb(V.vel,i,j,k,0,xd,yd,zd,3));
// 					dudz=-(7./3.)*Fb(V.vel,i,j,k,0,xd,yd,zd,3)+6.*Fb(V.vel,i,j,(k+1),0,xd,yd,zd,3)-3.*Fb(V.vel,i,j,(k+2),0,xd,yd,zd,3)-(2./3.)*Fb(V.vel,i,j,(k+3),0,xd,yd,zd,3)+3.*tmp;
					dudz = (GR*(Fb(V.vel,i,j,(k+1),0,xd,yd,zd,3)-Fb(V.vel,i,j,k,0,xd,yd,zd,3))+coef*(Fb(V.vel,i,j,k,0,xd,yd,zd,3)-Fb(V.velc,(i/GR),(j/GR),(zdc-3),0,xdc,ydc,zdc,3)))*1./(1.+GR);
					dvdz = (GR*(Fb(V.vel,i,j,(k+1),1,xd,yd,zd,3)-Fb(V.vel,i,j,k,1,xd,yd,zd,3))+coef*(Fb(V.vel,i,j,k,1,xd,yd,zd,3)-Fb(V.velc,(i/GR),(j/GR),(zdc-3),1,xdc,ydc,zdc,3)))*1./(1.+GR);
					dwdz = (GR*(Fb(V.vel,i,j,(k+1),2,xd,yd,zd,3)-Fb(V.vel,i,j,k,2,xd,yd,zd,3))+coef*(Fb(V.vel,i,j,k,2,xd,yd,zd,3)-Fb(V.velc,(i/GR),(j/GR),(zdc-3),2,xdc,ydc,zdc,3)))*1./(1.+GR);
				}
				else if (k==(zd-2))
				{
					tmp=(2./3.)*(Fb(V.vel,i,j,(k-1),2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-3),2,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,k,2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-4),2,xd,yd,zd,3));
					dwdz=(7./3.)*Fb(V.vel,i,j,k,2,xd,yd,zd,3)-6.*Fb(V.vel,i,j,(k-1),2,xd,yd,zd,3)+3.*Fb(V.vel,i,j,(k-2),2,xd,yd,zd,3)+(2./3.)*Fb(V.vel,i,j,(k-3),2,xd,yd,zd,3)+3.*tmp;
					tmp=(2./3.)*(Fb(V.vel,i,j,(k-1),1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-3),1,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,k,1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-4),1,xd,yd,zd,3));
					dvdz=(7./3.)*Fb(V.vel,i,j,k,1,xd,yd,zd,3)-6.*Fb(V.vel,i,j,(k-1),1,xd,yd,zd,3)+3.*Fb(V.vel,i,j,(k-2),1,xd,yd,zd,3)+(2./3.)*Fb(V.vel,i,j,(k-3),1,xd,yd,zd,3)+3.*tmp;
					tmp=(2./3.)*(Fb(V.vel,i,j,(k-1),0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-3),0,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,k,0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-4),0,xd,yd,zd,3));
					dudz=(7./3.)*Fb(V.vel,i,j,k,0,xd,yd,zd,3)-6.*Fb(V.vel,i,j,(k-1),0,xd,yd,zd,3)+3.*Fb(V.vel,i,j,(k-2),0,xd,yd,zd,3)+(2./3.)*Fb(V.vel,i,j,(k-3),0,xd,yd,zd,3)+3.*tmp;
				}
				else if (k==1)
				{
// 					tmp=(2./3.)*(Fb(V.vel,i,j,(k+2),2,xd,yd,zd,3)-Fb(V.vel,i,j,k,2,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+3),2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),2,xd,yd,zd,3));
// 					dwdz=-(1./6.)*Fb(V.vel,i,j,(k-1),2,xd,yd,zd,3)-(3./2.)*Fb(V.vel,i,j,k,2,xd,yd,zd,3)+(3./2.)*Fb(V.vel,i,j,(k+1),2,xd,yd,zd,3)+(1./6.)*Fb(V.vel,i,j,(k+2),2,xd,yd,zd,3)-tmp;
// 					tmp=(2./3.)*(Fb(V.vel,i,j,(k+2),1,xd,yd,zd,3)-Fb(V.vel,i,j,k,1,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+3),1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),1,xd,yd,zd,3));
// 					dvdz=-(1./6.)*Fb(V.vel,i,j,(k-1),1,xd,yd,zd,3)-(3./2.)*Fb(V.vel,i,j,k,1,xd,yd,zd,3)+(3./2.)*Fb(V.vel,i,j,(k+1),1,xd,yd,zd,3)+(1./6.)*Fb(V.vel,i,j,(k+2),1,xd,yd,zd,3)-tmp;
// 					tmp=(2./3.)*(Fb(V.vel,i,j,(k+2),0,xd,yd,zd,3)-Fb(V.vel,i,j,k,0,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+3),0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),0,xd,yd,zd,3));
// 					dudz=-(1./6.)*Fb(V.vel,i,j,(k-1),0,xd,yd,zd,3)-(3./2.)*Fb(V.vel,i,j,k,0,xd,yd,zd,3)+(3./2.)*Fb(V.vel,i,j,(k+1),0,xd,yd,zd,3)+(1./6.)*Fb(V.vel,i,j,(k+2),0,xd,yd,zd,3)-tmp;
					dudz = 0.5*(Fb(V.vel,i,j,(k+1),0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),0,xd,yd,zd,3));
					dvdz = 0.5*(Fb(V.vel,i,j,(k+1),1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),1,xd,yd,zd,3));
					dwdz = 0.5*(Fb(V.vel,i,j,(k+1),2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),2,xd,yd,zd,3));
				}
				else if (k==(zd-3))
				{
// 					tmp=(2./3.)*(Fb(V.vel,i,j,k,2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-2),2,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+1),2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-3),2,xd,yd,zd,3));
// 					dwdz=(1./6.)*Fb(V.vel,i,j,(k+1),2,xd,yd,zd,3)+(3./2.)*Fb(V.vel,i,j,k,2,xd,yd,zd,3)-(3./2.)*Fb(V.vel,i,j,(k-1),2,xd,yd,zd,3)-(1./6.)*Fb(V.vel,i,j,(k-2),2,xd,yd,zd,3)-tmp;
// 					tmp=(2./3.)*(Fb(V.vel,i,j,k,1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-2),1,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+1),1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-3),1,xd,yd,zd,3));
// 					dvdz=(1./6.)*Fb(V.vel,i,j,(k+1),1,xd,yd,zd,3)+(3./2.)*Fb(V.vel,i,j,k,1,xd,yd,zd,3)-(3./2.)*Fb(V.vel,i,j,(k-1),1,xd,yd,zd,3)-(1./6.)*Fb(V.vel,i,j,(k-2),1,xd,yd,zd,3)-tmp;
// 					tmp=(2./3.)*(Fb(V.vel,i,j,k,0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-2),0,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+1),0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-3),0,xd,yd,zd,3));
// 					dudz=(1./6.)*Fb(V.vel,i,j,(k+1),0,xd,yd,zd,3)+(3./2.)*Fb(V.vel,i,j,k,0,xd,yd,zd,3)-(3./2.)*Fb(V.vel,i,j,(k-1),0,xd,yd,zd,3)-(1./6.)*Fb(V.vel,i,j,(k-2),0,xd,yd,zd,3)-tmp;
					dudz = 0.5*(Fb(V.vel,i,j,(k+1),0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),0,xd,yd,zd,3));
					dvdz = 0.5*(Fb(V.vel,i,j,(k+1),1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),1,xd,yd,zd,3));
					dwdz = 0.5*(Fb(V.vel,i,j,(k+1),2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),2,xd,yd,zd,3));
				}
				else
				{
//					dwdz=(2./3.)*(Fb(V.vel,i,j,(k+1),2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),2,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+2),2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-2),2,xd,yd,zd,3));
//					dvdz=(2./3.)*(Fb(V.vel,i,j,(k+1),1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),1,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+2),1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-2),1,xd,yd,zd,3));
//					dudz=(2./3.)*(Fb(V.vel,i,j,(k+1),0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),0,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+2),0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-2),0,xd,yd,zd,3));
					dudz = 0.5*(Fb(V.vel,i,j,(k+1),0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),0,xd,yd,zd,3));
                                        dvdz = 0.5*(Fb(V.vel,i,j,(k+1),1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),1,xd,yd,zd,3));
                                        dwdz = 0.5*(Fb(V.vel,i,j,(k+1),2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),2,xd,yd,zd,3));
				}
				
				ox=(dwdy-dvdz);
				oy=(dudz-dwdx);
				oz=(dvdx-dudy);

				*(wx+zdf+zdc-4+k) += ox*ox;
             			*(wy+zdf+zdc-4+k) += oy*oy;
           			*(wz+zdf+zdc-4+k) += oz*oz;
				
				*(Paiuu+zdf+zdc-4+k) += dfl*dudx;
				*(Paivv+zdf+zdc-4+k) += dfl*dvdy;
				*(Paiww+zdf+zdc-4+k) += dfl*dwdz;
				*(Paiuw+zdf+zdc-4+k) += dfl*(dudz+dwdx);
				
				*(usuw+zdf+zdc-4+k) += 0.5*u*(dudz+dwdx);
				*(vsvw+zdf+zdc-4+k) += 0.5*v*(dvdz+dwdy);
				*(wsww+zdf+zdc-4+k) += 0.5*w*(dwdz+dwdz);
				*(wsuw+zdf+zdc-4+k) += 0.5*w*(dudz+dwdx);
				*(usww+zdf+zdc-4+k) += 0.5*u*(dwdz+dwdz);
				
				*(Epuu+zdf+zdc-4+k) += dudx*(dudx+dudx)+dudy*(dudy+dvdx)+dudz*(dudz+dwdx);
				*(Epvv+zdf+zdc-4+k) += dvdx*(dvdx+dudy)+dvdy*(dvdy+dvdy)+dvdz*(dvdz+dwdy);
				*(Epww+zdf+zdc-4+k) += dwdx*(dwdx+dudz)+dwdy*(dwdy+dvdz)+dwdz*(dwdz+dwdz);
				*(Epuw+zdf+zdc-4+k) += dudx*(dwdx+dudz)+dudy*(dwdy+dvdz)+dudz*(dwdz+dwdz)+dwdx*(dudx+dudx)+dwdy*(dudy+dvdx)+dwdz*(dudz+dwdx);
			}
		}
	}
	
//	num=(xdf-1)*ydf*pd.numproc;
	V.vel=NULL;
	V.den=NULL;
	xd=yd=zd=0;
	for (k=0;k<NZ;k++)
	{
		if (k<zdf)
		  num=(xdf-1)*ydf*pd.numproc;
		else if ((k<(zdf+zdc-4))&&(k>=zdf))
		  num=(xdc-1)*ydc*pd.numproc;
		else if (k>(zdf+zdc-5))
		  num=(xdf-1)*ydf*pd.numproc;
		
		*(uu+k) /= (double)(num);
		*(uv+k) /= (double)(num);
		*(uw+k) /= (double)(num);
		*(vv+k) /= (double)(num);
		*(vw+k) /= (double)(num);
		*(ww+k) /= (double)(num);
		*(dd+k) /= (double)(num);
		
		*(SKEWu+k) /= (double)(num);
		*(SKEWv+k) /= (double)(num);
		*(SKEWw+k) /= (double)(num);
		*(SKEWp+k) /= (double)(num);
		*(SKEWuw+k) /= (double)(num);
		*(SKEWuv+k) /= (double)(num);
		*(SKEWvw+k) /= (double)(num);
		
		*(KURTu+k) /= (double)(num);
		*(KURTv+k) /= (double)(num);
		*(KURTw+k) /= (double)(num);
		*(KURTp+k) /= (double)(num);
		*(KURTuw+k) /= (double)(num);
		*(KURTuv+k) /= (double)(num);
		*(KURTvw+k) /= (double)(num);

		*(wx+k) /= (double)(num);
     		*(wy+k) /= (double)(num);
     		*(wz+k) /= (double)(num);
		
		*(uuw+k) /= (double)(num);
		*(vvw+k) /= (double)(num);
		*(www+k) /= (double)(num);
		*(uww+k) /= (double)(num);

		*(pu+k) /= (double)(num);
		*(pv+k) /= (double)(num);
		*(pw+k) /= (double)(num);

		*(Paiuu+k) /= (double)(num);
		*(Paivv+k) /= (double)(num);
		*(Paiww+k) /= (double)(num);
		*(Paiuw+k) /= (double)(num);

		*(Epuu+k) /= (double)(num);
		*(Epvv+k) /= (double)(num);
		*(Epww+k) /= (double)(num);
		*(Epuw+k) /= (double)(num);
		
		*(usuw+k) /= (double)(num);
		*(vsvw+k) /= (double)(num);
		*(wsww+k) /= (double)(num);
		*(wsuw+k) /= (double)(num);
		*(usww+k) /= (double)(num);
	}
	
// ***********************************************************************************
	zd=NZ;
	
	MPI_Reduce(uu,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(uu,DUMMY,zd*sizeof(DP));

	MPI_Reduce(uv,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(uv,DUMMY,zd*sizeof(DP));

	MPI_Reduce(uw,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(uw,DUMMY,zd*sizeof(DP));

	MPI_Reduce(vv,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(vv,DUMMY,zd*sizeof(DP));

	MPI_Reduce(vw,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(vw,DUMMY,zd*sizeof(DP));

	MPI_Reduce(ww,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(ww,DUMMY,zd*sizeof(DP));

	MPI_Reduce(dd,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(dd,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(SKEWu,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(SKEWu,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(SKEWv,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(SKEWv,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(SKEWw,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(SKEWw,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(SKEWp,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(SKEWp,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(SKEWuw,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(SKEWuw,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(SKEWuv,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(SKEWuv,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(SKEWvw,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(SKEWvw,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(KURTu,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(KURTu,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(KURTv,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(KURTv,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(KURTw,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(KURTw,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(KURTp,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(KURTp,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(KURTuw,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(KURTuw,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(KURTvw,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(KURTvw,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(KURTuv,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(KURTuv,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(wx,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(wx,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(wy,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(wy,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(wz,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(wz,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(uuw,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(uuw,DUMMY,zd*sizeof(DP));

	MPI_Reduce(vvw,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(vvw,DUMMY,zd*sizeof(DP));

	MPI_Reduce(www,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(www,DUMMY,zd*sizeof(DP));

	MPI_Reduce(uww,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(uww,DUMMY,zd*sizeof(DP));

	MPI_Reduce(pu,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(pu,DUMMY,zd*sizeof(DP));

	MPI_Reduce(pv,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(pv,DUMMY,zd*sizeof(DP));

	MPI_Reduce(pw,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(pw,DUMMY,zd*sizeof(DP));

	MPI_Reduce(Paiuu,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(Paiuu,DUMMY,zd*sizeof(DP));

	MPI_Reduce(Paivv,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(Paivv,DUMMY,zd*sizeof(DP));

	MPI_Reduce(Paiww,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(Paiww,DUMMY,zd*sizeof(DP));

	MPI_Reduce(Paiuw,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(Paiuw,DUMMY,zd*sizeof(DP));

	MPI_Reduce(Epuu,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(Epuu,DUMMY,zd*sizeof(DP));

	MPI_Reduce(Epvv,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(Epvv,DUMMY,zd*sizeof(DP));

	MPI_Reduce(Epww,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(Epww,DUMMY,zd*sizeof(DP));

	MPI_Reduce(Epuw,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(Epuw,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(usuw,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(usuw,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(vsvw,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(vsvw,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(wsww,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(wsww,DUMMY,zd*sizeof(DP));

	MPI_Reduce(wsuw,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(wsuw,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(usww,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(usww,DUMMY,zd*sizeof(DP));

	zd=0;
//*****************************************************************************************************************

	for (k=1;k<(NZ-1);k++)
	{
		siu = *(UTIN+k) = sqrt((*(uu+k)));
		siv = *(VTIN+k) = sqrt((*(vv+k)));
		siw = *(WTIN+k) = sqrt((*(ww+k)));
		sip = *(DTIN+k) = sqrt((*(dd+k)));
		
		*(KE+k) += ((*(uu+k))+(*(vv+k))+(*(ww+k)));
		
		*(wx+k) = sqrt((*(wx+k)));
		*(wy+k) = sqrt((*(wy+k)));
		*(wz+k) = sqrt((*(wz+k)));
		
		*(SKEWu+k) /= (siu*siu*siu);
		*(SKEWv+k) /= (siv*siv*siv);
		*(SKEWw+k) /= (siw*siw*siw);
		*(SKEWp+k) /= (sip*sip*sip);
		*(SKEWuw+k) /= (siu*siu*siu*siw*siw*siw);
		*(SKEWvw+k) /= (siv*siv*siv*siw*siw*siw);
		*(SKEWuv+k) /= (siu*siu*siu*siv*siv*siv);
		
		*(KURTu+k) /= (siu*siu*siu*siu);
		*(KURTv+k) /= (siv*siv*siv*siv);
		*(KURTw+k) /= (siw*siw*siw*siw);
		*(KURTp+k) /= (sip*sip*sip*sip);
		*(KURTuw+k) /= (siu*siu*siu*siu*siw*siw*siw*siw);
		*(KURTuv+k) /= (siu*siu*siu*siu*siv*siv*siv*siv);
		*(KURTvw+k) /= (siv*siv*siv*siv*siw*siw*siw*siw);
		
		*(Epuu+k) *= nu;
		*(Epvv+k) *= nu;
		*(Epww+k) *= nu;
		*(Epuw+k) *= (2.*nu);
		*(Epii+k) = (*(Epuu+k))+(*(Epvv+k))+(*(Epww+k));
	}	
	
	
	for (k=1;k<zdf;k++)
	{
		*(Pii+k) = -(*(uw+k))*(*(dUdz+k));
		*(Puu+k) = 2.*(*(Pii+k));
		*(Pvv+k) = 0.;
		*(Pww+k) = 0.;
		*(Puw+k) = -(*(ww+k))*(*(dUdz+k));

		*(Tuupai+k) = 0.;
		*(Tvvpai+k) = 0.;
		
		if (k==1)	
		{
			*(Tuur+k) = -DDZ1F(uuw,k);
			*(Tvvr+k) = -DDZ1F(vvw,k);
			*(Twwr+k) = -DDZ1F(www,k);
			*(Tuwr+k) = -DDZ1F(uww,k);
			
			*(Tiipai+k) = -DDZ1F(pw,k);

			*(Twwpai+k) = -2.*DDZ1F(pw,k);
			*(Tuwpai+k) = -DDZ1F(pu,k);

			*(Tuus+k) = 2.*nu*DDZ1F(usuw,k);
			*(Tvvs+k) = 2.*nu*DDZ1F(vsvw,k);
			*(Twws+k) = 2.*nu*DDZ1F(wsww,k);
			*(Tuws+k) = 2.*nu*DDZ1F(wsuw,k) + 2.*nu*DDZ1F(wsww,k);
		}
		else if (k==(zdf-1))
		{
// 			*(Tuur+k) = -DDZ1B(uuw,k);
// 			*(Tvvr+k) = -DDZ1B(vvw,k);
// 			*(Twwr+k) = -DDZ1B(www,k);
// 			*(Tuwr+k) = -DDZ1B(uww,k);
// 			
// 			*(Tiipai+k) = -DDZ1B(pw,k);
// 
// 			*(Twwpai+k) = -2.*DDZ1B(pw,k);
// 			*(Tuwpai+k) = -DDZ1B(pu,k);
// 
// 			*(Tuus+k) = 2.*nu*DDZ1B(usuw,k);
// 			*(Tvvs+k) = 2.*nu*DDZ1B(vsvw,k);
// 			*(Twws+k) = 2.*nu*DDZ1B(wsww,k);
// 			*(Tuws+k) = 2.*nu*DDZ1B(wsuw,k) + 2.*nu*DDZ1B(wsww,k);
			*(Tuur+k) = -DDZ2CNUl(uuw,k,GR);
			*(Tvvr+k) = -DDZ2CNUl(vvw,k,GR);
			*(Twwr+k) = -DDZ2CNUl(www,k,GR);
			*(Tuwr+k) = -DDZ2CNUl(uww,k,GR);
			
			*(Tiipai+k) = -DDZ2CNUl(pw,k,GR);

			*(Twwpai+k) = -2.*DDZ2CNUl(pw,k,GR);
			*(Tuwpai+k) = -DDZ2CNUl(pu,k,GR);

			*(Tuus+k) = 2.*nu*DDZ2CNUl(usuw,k,GR);
			*(Tvvs+k) = 2.*nu*DDZ2CNUl(vsvw,k,GR);
			*(Twws+k) = 2.*nu*DDZ2CNUl(wsww,k,GR);
			*(Tuws+k) = 2.*nu*DDZ2CNUl(wsuw,k,GR) + 2.*nu*DDZ2CNUl(wsww,k,GR);
		}
		else if ((k==2)||(k==(zdf-2)))
		{
			*(Tuur+k) = -DDZ2C(uuw,k);
			*(Tvvr+k) = -DDZ2C(vvw,k);
			*(Twwr+k) = -DDZ2C(www,k);
			*(Tuwr+k) = -DDZ2C(uww,k);
			
			*(Tiipai+k) = -DDZ2C(pw,k);

			*(Twwpai+k) = -2.*DDZ2C(pw,k);
			*(Tuwpai+k) = -DDZ2C(pu,k);
			
			*(Tuus+k) = 2.*nu*DDZ2C(usuw,k);
			*(Tvvs+k) = 2.*nu*DDZ2C(vsvw,k);
			*(Twws+k) = 2.*nu*DDZ2C(wsww,k);
			*(Tuws+k) = 2.*nu*DDZ2C(wsuw,k) + 2.*nu*DDZ2C(wsww,k);
		}
		else
		{
			*(Tuur+k) = -0.5*DDZ2C(uuw,k);
			*(Tvvr+k) = -0.5*DDZ2C(vvw,k);
			*(Twwr+k) = -0.5*DDZ2C(www,k);
			*(Tuwr+k) = -DDZ2C(uww,k);
			
			*(Tiipai+k) = -DDZ2C(pw,k);

			*(Twwpai+k) = -2.*DDZ2C(pw,k);
			*(Tuwpai+k) = -DDZ2C(pu,k);

			*(Tuus+k) = 2.*nu*DDZ2C(usuw,k);
			*(Tvvs+k) = 2.*nu*DDZ2C(vsvw,k);
			*(Twws+k) = 2.*nu*DDZ2C(wsww,k);
			*(Tuws+k) = 2.*nu*DDZ2C(wsuw,k) + 2.*nu*DDZ2C(wsww,k);
		}
		Tiis[k] = 0.5*(Tuus[k]+Tvvs[k]+Twws[k]);
		Tiir[k] = 0.5*(Tuur[k]+Tvvr[k]+Twwr[k]);
	}
	
	for (k=2;k<(zdc-2);k++)
	{
		*(Pii+zdf-2+k) = -(*(uw+zdf-2+k))*(*(dUdz+zdf-2+k));
		*(Puu+zdf-2+k) = 2.*(*(Pii+zdf-2+k));
		*(Pvv+zdf-2+k) = 0.;
		*(Pww+zdf-2+k) = 0.;
		*(Puw+zdf-2+k) = -(*(ww+zdf-2+k))*(*(dUdz+zdf-2+k));

		*(Tuupai+zdf-2+k) = 0.;
		*(Tvvpai+zdf-2+k) = 0.;
		
		if (k==2)
		{
			*(Tuur+zdf-2+k) = -coef*DDZ2C(uuw,(zdf-2+k));
			*(Tvvr+zdf-2+k) = -coef*DDZ2C(vvw,(zdf-2+k));
			*(Twwr+zdf-2+k) = -coef*DDZ2C(www,(zdf-2+k));
			*(Tuwr+zdf-2+k) = -coef*DDZ2C(uww,(zdf-2+k));
			
			*(Tiipai+zdf-2+k) = -coef*DDZ2C(pw,(zdf-2+k));

			*(Twwpai+zdf-2+k) = -2.*coef*DDZ2C(pw,(zdf-2+k));
			*(Tuwpai+zdf-2+k) = -coef*DDZ2C(pu,(zdf-2+k));

			*(Tuus+zdf-2+k) = 2.*nu*coef*DDZ2C(usuw,(zdf-2+k));
			*(Tvvs+zdf-2+k) = 2.*nu*coef*DDZ2C(vsvw,(zdf-2+k));
			*(Twws+zdf-2+k) = 2.*nu*coef*DDZ2C(wsww,(zdf-2+k));
			*(Tuws+zdf-2+k) = 2.*nu*coef*DDZ2C(wsuw,(zdf-2+k)) + 2.*nu*coef*DDZ2C(wsww,(zdf-2+k));
		}
		else if ((k==3)||(k==(zdc-4)))
		{
			*(Tuur+zdf-2+k) = -coef*DDZ2C(uuw,(zdf-2+k));
			*(Tvvr+zdf-2+k) = -coef*DDZ2C(vvw,(zdf-2+k));
			*(Twwr+zdf-2+k) = -coef*DDZ2C(www,(zdf-2+k));
			*(Tuwr+zdf-2+k) = -coef*DDZ2C(uww,(zdf-2+k));
			
			*(Tiipai+zdf-2+k) = -coef*DDZ2C(pw,(zdf-2+k));

			*(Twwpai+zdf-2+k) = -2.*coef*DDZ2C(pw,(zdf-2+k));
			*(Tuwpai+zdf-2+k) = -coef*DDZ2C(pu,(zdf-2+k));
			
			*(Tuus+zdf-2+k) = 2.*nu*coef*DDZ2C(usuw,(zdf-2+k));
			*(Tvvs+zdf-2+k) = 2.*nu*coef*DDZ2C(vsvw,(zdf-2+k));
			*(Twws+zdf-2+k) = 2.*nu*coef*DDZ2C(wsww,(zdf-2+k));
			*(Tuws+zdf-2+k) = 2.*nu*coef*DDZ2C(wsuw,(zdf-2+k)) + 2.*nu*coef*DDZ2C(wsww,(zdf-2+k));
		}
		else if (k==(zdc-3))
		{
			*(Tuur+zdf-2+k) = -coef*DDZ2C(uuw,(zdf-2+k));
			*(Tvvr+zdf-2+k) = -coef*DDZ2C(vvw,(zdf-2+k));
			*(Twwr+zdf-2+k) = -coef*DDZ2C(www,(zdf-2+k));
			*(Tuwr+zdf-2+k) = -coef*DDZ2C(uww,(zdf-2+k));
			
			*(Tiipai+zdf-2+k) = -coef*DDZ2C(pw,(zdf-2+k));

			*(Twwpai+zdf-2+k) = -2.*coef*DDZ2C(pw,(zdf-2+k));
			*(Tuwpai+zdf-2+k) = -coef*DDZ2C(pu,(zdf-2+k));

			*(Tuus+zdf-2+k) = 2.*nu*coef*DDZ2C(usuw,(zdf-2+k));
			*(Tvvs+zdf-2+k) = 2.*nu*coef*DDZ2C(vsvw,(zdf-2+k));
			*(Twws+zdf-2+k) = 2.*nu*coef*DDZ2C(wsww,(zdf-2+k));
			*(Tuws+zdf-2+k) = 2.*nu*coef*DDZ2C(wsuw,(zdf-2+k)) + 2.*nu*coef*DDZ2C(wsww,(zdf-2+k));
		}
		else
		{
			*(Tuur+zdf-2+k) = -coef*0.5*DDZ2C(uuw,(zdf-2+k));
			*(Tvvr+zdf-2+k) = -coef*0.5*DDZ2C(vvw,(zdf-2+k));
			*(Twwr+zdf-2+k) = -coef*0.5*DDZ2C(www,(zdf-2+k));
			*(Tuwr+zdf-2+k) = -coef*DDZ2C(uww,(zdf-2+k));
			
			*(Tiipai+zdf-2+k) = -coef*DDZ2C(pw,(zdf-2+k));

			*(Twwpai+zdf-2+k) = -2.*coef*DDZ2C(pw,(zdf-2+k));
			*(Tuwpai+zdf-2+k) = -coef*DDZ2C(pu,(zdf-2+k));

			*(Tuus+zdf-2+k) = coef*2.*nu*DDZ2C(usuw,(zdf-2+k));
			*(Tvvs+zdf-2+k) = coef*2.*nu*DDZ2C(vsvw,(zdf-2+k));
			*(Twws+zdf-2+k) = coef*2.*nu*DDZ2C(wsww,(zdf-2+k));
			*(Tuws+zdf-2+k) = coef*2.*nu*DDZ2C(wsuw,(zdf-2+k)) + coef*2.*nu*DDZ2C(wsww,(zdf-2+k));
		}

		Tiis[zdf-2+k] = 0.5*(Tuus[(zdf-2+k)]+Tvvs[(zdf-2+k)]+Twws[(zdf-2+k)]);
		Tiir[zdf-2+k] = 0.5*(Tuur[(zdf-2+k)]+Tvvr[(zdf-2+k)]+Twwr[(zdf-2+k)]);
	}
	
	for (k=0;k<(zdf-1);k++)
	{
		*(Pii+zdf+zdc-4+k) = -(*(uw+zdf+zdc-4+k))*(*(dUdz+zdf+zdc-4+k));
		*(Puu+zdf+zdc-4+k) = 2.*(*(Pii+zdf+zdc-4+k));
		*(Pvv+zdf+zdc-4+k) = 0.;
		*(Pww+zdf+zdc-4+k) = 0.;
		*(Puw+zdf+zdc-4+k) = -(*(ww+zdf+zdc-4+k))*(*(dUdz+zdf+zdc-4+k));

		*(Tuupai+zdf+zdc-4+k) = 0.;
		*(Tvvpai+zdf+zdc-4+k) = 0.;
		
		if (k==0)	
		{
			*(Tuur+zdf+zdc-4+k) = -DDZ2CNUu(uuw,(zdf+zdc-4+k),GR);
			*(Tvvr+zdf+zdc-4+k) = -DDZ2CNUu(vvw,(zdf+zdc-4+k),GR);
			*(Twwr+zdf+zdc-4+k) = -DDZ2CNUu(www,(zdf+zdc-4+k),GR);
			*(Tuwr+zdf+zdc-4+k) = -DDZ2CNUu(uww,(zdf+zdc-4+k),GR);
			
			*(Tiipai+zdf+zdc-4+k) = -DDZ2CNUu(pw,(zdf+zdc-4+k),GR);

			*(Twwpai+zdf+zdc-4+k) = -2.*DDZ2CNUu(pw,(zdf+zdc-4+k),GR);
			*(Tuwpai+zdf+zdc-4+k) = -DDZ2CNUu(pu,(zdf+zdc-4+k),GR);

			*(Tuus+zdf+zdc-4+k) = 2.*nu*DDZ2CNUu(usuw,(zdf+zdc-4+k),GR);
			*(Tvvs+zdf+zdc-4+k) = 2.*nu*DDZ2CNUu(vsvw,(zdf+zdc-4+k),GR);
			*(Twws+zdf+zdc-4+k) = 2.*nu*DDZ2CNUu(wsww,(zdf+zdc-4+k),GR);
			*(Tuws+zdf+zdc-4+k) = 2.*nu*DDZ2CNUu(wsuw,(zdf+zdc-4+k),GR) + 2.*nu*DDZ2CNUu(wsww,(zdf+zdc-4+k),GR);
		}
		else if (k==(zdf-2))
		{
			*(Tuur+zdf+zdc-4+k) = -DDZ1B(uuw,(zdf+zdc-4+k));
			*(Tvvr+zdf+zdc-4+k) = -DDZ1B(vvw,(zdf+zdc-4+k));
			*(Twwr+zdf+zdc-4+k) = -DDZ1B(www,(zdf+zdc-4+k));
			*(Tuwr+zdf+zdc-4+k) = -DDZ1B(uww,(zdf+zdc-4+k));
			
			*(Tiipai+zdf+zdc-4+k) = -DDZ1B(pw,(zdf+zdc-4+k));

			*(Twwpai+zdf+zdc-4+k) = -2.*DDZ1B(pw,(zdf+zdc-4+k));
			*(Tuwpai+zdf+zdc-4+k) = -DDZ1B(pu,(zdf+zdc-4+k));

			*(Tuus+zdf+zdc-4+k) = 2.*nu*DDZ1B(usuw,(zdf+zdc-4+k));
			*(Tvvs+zdf+zdc-4+k) = 2.*nu*DDZ1B(vsvw,(zdf+zdc-4+k));
			*(Twws+zdf+zdc-4+k) = 2.*nu*DDZ1B(wsww,(zdf+zdc-4+k));
			*(Tuws+zdf+zdc-4+k) = 2.*nu*DDZ1B(wsuw,(zdf+zdc-4+k)) + 2.*nu*DDZ1B(wsww,(zdf+zdc-4+k));
		}
		else if ((k==2)||(k==(zdf-3)))
		{
			*(Tuur+zdf+zdc-4+k) = -DDZ2C(uuw,(zdf+zdc-4+k));
			*(Tvvr+zdf+zdc-4+k) = -DDZ2C(vvw,(zdf+zdc-4+k));
			*(Twwr+zdf+zdc-4+k) = -DDZ2C(www,(zdf+zdc-4+k));
			*(Tuwr+zdf+zdc-4+k) = -DDZ2C(uww,(zdf+zdc-4+k));
			
			*(Tiipai+zdf+zdc-4+k) = -DDZ2C(pw,(zdf+zdc-4+k));

			*(Twwpai+zdf+zdc-4+k) = -2.*DDZ2C(pw,(zdf+zdc-4+k));
			*(Tuwpai+zdf+zdc-4+k) = -DDZ2C(pu,(zdf+zdc-4+k));
			
			*(Tuus+zdf+zdc-4+k) = 2.*nu*DDZ2C(usuw,(zdf+zdc-4+k));
			*(Tvvs+zdf+zdc-4+k) = 2.*nu*DDZ2C(vsvw,(zdf+zdc-4+k));
			*(Twws+zdf+zdc-4+k) = 2.*nu*DDZ2C(wsww,(zdf+zdc-4+k));
			*(Tuws+zdf+zdc-4+k) = 2.*nu*DDZ2C(wsuw,(zdf+zdc-4+k)) + 2.*nu*DDZ2C(wsww,(zdf+zdc-4+k));
		}
		else
		{
			*(Tuur+zdf+zdc-4+k) = -0.5*DDZ2C(uuw,(zdf+zdc-4+k));
			*(Tvvr+zdf+zdc-4+k) = -0.5*DDZ2C(vvw,(zdf+zdc-4+k));
			*(Twwr+zdf+zdc-4+k) = -0.5*DDZ2C(www,(zdf+zdc-4+k));
			*(Tuwr+zdf+zdc-4+k) = -DDZ2C(uww,(zdf+zdc-4+k));
			
			*(Tiipai+zdf+zdc-4+k) = -DDZ2C(pw,(zdf+zdc-4+k));

			*(Twwpai+zdf+zdc-4+k) = -2.*DDZ2C(pw,(zdf+zdc-4+k));
			*(Tuwpai+zdf+zdc-4+k) = -DDZ2C(pu,(zdf+zdc-4+k));

			*(Tuus+zdf+zdc-4+k) = 2.*nu*DDZ2C(usuw,(zdf+zdc-4+k));
			*(Tvvs+zdf+zdc-4+k) = 2.*nu*DDZ2C(vsvw,(zdf+zdc-4+k));
			*(Twws+zdf+zdc-4+k) = 2.*nu*DDZ2C(wsww,(zdf+zdc-4+k));
			*(Tuws+zdf+zdc-4+k) = 2.*nu*DDZ2C(wsuw,(zdf+zdc-4+k)) + 2.*nu*DDZ2C(wsww,(zdf+zdc-4+k));
		}
		Tiis[zdf+zdc-4+k] = 0.5*(Tuus[zdf+zdc-4+k]+Tvvs[zdf+zdc-4+k]+Twws[zdf+zdc-4+k]);
		Tiir[zdf+zdc-4+k] = 0.5*(Tuur[zdf+zdc-4+k]+Tvvr[zdf+zdc-4+k]+Twwr[zdf+zdc-4+k]);
	}

	zd=NZ;
	if (pd.myrank==0)
	{
		sprintf(fn,"turb-fld-2d.%.3d.%.4d",pd.myrank,ts);
		sv=fopen(fn,"wb");
		fwrite(UAVG,sizeof(DP),zd,sv);
		fwrite(VAVG,sizeof(DP),zd,sv);
		fwrite(WAVG,sizeof(DP),zd,sv);
		fwrite(UTIN,sizeof(DP),zd,sv);
		fwrite(VTIN,sizeof(DP),zd,sv);
		fwrite(WTIN,sizeof(DP),zd,sv);
//		fwrite(uv,sizeof(DP),zd,sv);
		fwrite(uw,sizeof(DP),zd,sv);
//		fwrite(vw,sizeof(DP),zd,sv);
		fwrite(DAVG,sizeof(DP),zd,sv);
		fwrite(DTIN,sizeof(DP),zd,sv);
		fwrite(SKEWu,sizeof(DP),zd,sv);
		fwrite(SKEWv,sizeof(DP),zd,sv);
		fwrite(SKEWw,sizeof(DP),zd,sv);
//		fwrite(SKEWuv,sizeof(DP),zd,sv);
		fwrite(SKEWuw,sizeof(DP),zd,sv);
//		fwrite(SKEWvw,sizeof(DP),zd,sv);
		fwrite(SKEWp,sizeof(DP),zd,sv);
		fwrite(KURTu,sizeof(DP),zd,sv);
		fwrite(KURTv,sizeof(DP),zd,sv);
		fwrite(KURTw,sizeof(DP),zd,sv);
//		fwrite(KURTuv,sizeof(DP),zd,sv);
		fwrite(KURTuw,sizeof(DP),zd,sv);
//		fwrite(KURTvw,sizeof(DP),zd,sv);
		fwrite(KURTp,sizeof(DP),zd,sv);
//		fwrite(UW,sizeof(DP),zd,sv);
		fclose(sv);
		
		sprintf(fn,"vort-fluc-2d.%.3d.%.4d",pd.myrank,ts);	
		sv=fopen(fn,"wb");
		fwrite(wx,sizeof(DP),zd,sv);
		fwrite(wy,sizeof(DP),zd,sv);
		fwrite(wz,sizeof(DP),zd,sv);
		fclose(sv);
                
                sprintf(fn,"KE-Budget-2d.%.3d.%.4d",pd.myrank,ts);
                sv=fopen(fn,"wb");
		fwrite(Pii,sizeof(DP),zd,sv);
		fwrite(Puu,sizeof(DP),zd,sv);
		fwrite(Pvv,sizeof(DP),zd,sv);
		fwrite(Pww,sizeof(DP),zd,sv);
//		fwrite(Puv,sizeof(DP),zd,sv);
		fwrite(Puw,sizeof(DP),zd,sv);
//		fwrite(Pvw,sizeof(DP),zd,sv);
		
		fwrite(Tiir,sizeof(DP),zd,sv);
		fwrite(Tuur,sizeof(DP),zd,sv);
		fwrite(Tvvr,sizeof(DP),zd,sv);
		fwrite(Twwr,sizeof(DP),zd,sv);
//		fwrite(Tuvr,sizeof(DP),zd,sv);
		fwrite(Tuwr,sizeof(DP),zd,sv);
//		fwrite(Tvwr,sizeof(DP),zd,sv);
		
		fwrite(Tiipai,sizeof(DP),zd,sv);
		fwrite(Tuupai,sizeof(DP),zd,sv);
		fwrite(Tvvpai,sizeof(DP),zd,sv);
		fwrite(Twwpai,sizeof(DP),zd,sv);
//		fwrite(Tuvpai,sizeof(DP),zd,sv);
		fwrite(Tuwpai,sizeof(DP),zd,sv);
//		fwrite(Tvwpai,sizeof(DP),zd,sv);
		
		fwrite(Paiuu,sizeof(DP),zd,sv);
		fwrite(Paivv,sizeof(DP),zd,sv);
		fwrite(Paiww,sizeof(DP),zd,sv);
//		fwrite(Paiuv,sizeof(DP),zd,sv);
		fwrite(Paiuw,sizeof(DP),zd,sv);
//		fwrite(Paivw,sizeof(DP),zd,sv);
		
		fwrite(Tiis,sizeof(DP),zd,sv);
		fwrite(Tuus,sizeof(DP),zd,sv);
		fwrite(Tvvs,sizeof(DP),zd,sv);
		fwrite(Twws,sizeof(DP),zd,sv);
//		fwrite(Tuvs,sizeof(DP),zd,sv);
		fwrite(Tuws,sizeof(DP),zd,sv);
//		fwrite(Tvws,sizeof(DP),zd,sv);
		
		fwrite(Epii,sizeof(DP),zd,sv);
		fwrite(Epuu,sizeof(DP),zd,sv);
		fwrite(Epvv,sizeof(DP),zd,sv);
		fwrite(Epww,sizeof(DP),zd,sv);
//		fwrite(Epuv,sizeof(DP),zd,sv);
		fwrite(Epuw,sizeof(DP),zd,sv);
//		fwrite(Epvw,sizeof(DP),zd,sv);
		
		fwrite(KE,sizeof(DP),zd,sv);
		fclose(sv);
	}
	
	
	free(UAVG);
	free(VAVG);
	free(WAVG);
	free(DAVG);
	free(DUMMY);
	free(UW);
	
	free(UTIN);
	free(VTIN);
	free(WTIN);
	free(DTIN);
	
	free(uu);
	free(vv);
	free(ww);
	free(uw);
	free(vw);
	free(uv);
	free(dd);
	
	free(wx);
	free(wy);
	free(wz);
	
	free(SKEWu);
	free(SKEWv);
	free(SKEWw);
	free(SKEWp);
	free(SKEWuw);
	free(SKEWuv);
	free(SKEWvw);
	
	free(KURTu);
	free(KURTv);
	free(KURTw);
	free(KURTp);
	free(KURTuw);
	free(KURTuv);
	free(KURTvw);
	
	free(KE);
	
	free(Pii);
	free(Puu);
	free(Pvv);
	free(Pww);
	free(Puw);
	free(Puv);
	free(Pvw);
	
	free(Tiir);
	free(Tiipai);
	free(Tiis);
	free(Epii);

	free(Paiuu);
	free(Paivv);
	free(Paiww);
	free(Paiuw);
	free(Paiuv);
	free(Paivw);

	free(Tiirt);
	free(Tuur);
	free(Tvvr);
	free(Twwr);
	free(Tuwr);
	free(Tuvr);
	free(Tvwr);

	free(Tuus);
	free(Tvvs);
	free(Twws);
	free(Tuws);
	free(Tuvs);
	free(Tvws);

	free(Tuupai);
	free(Tvvpai);
	free(Twwpai);
	free(Tuwpai);
	free(Tuvpai);
	free(Tvwpai);

	free(Epuu);
	free(Epvv);
	free(Epww);
	free(Epuw);
	free(Epuv);
	free(Epvw);
	
return V;
}
