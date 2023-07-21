/********************************************************************************
 * Copyright (c) 2012								  *
 * Amirreza rastegari								  *
 * amirreza@umich.edu								  *
 * Creates the required files to plot 2 point velocity correlations		  *
 * It should be used with correlations.c					  *
 ********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define Reb 7200.0000
#define ubulk 0.0866034
#define TSTR 916
#define TEND 1195
#define ZPTGT 5    /* Z+ for the requested plane */
//if defined has priority
//#define ZOH 0.8   /* Z/H for the requested plane */
#define DT 1
#define TEXT "1DAVG"
// SLPF: 0 both, -1 no-slip, -2 slip
#define SLPF -3
#define SLPTXT "1DAVG"


typedef double DP;

#define xdv 512
#define ydv 256
#define zdv 197

/**************************************************************************
 * fine grid's parameters                                                 *
 * gratio : dx^c/dx^f ratio of the coarse to fine grid's resolution       *
 * ***********************************************************************/
#define gratio 4

// Number of z grid point on fine grid 
// a +1 is added inside the code to z for the wall points
#define zds 56

#define Fy(a,j,k,b,yd,zd,d) (*(a+(k*yd*d)+(j*d)+b))
#define Fx(a,i,k,b,xd,zd,d) (*(a+(k*xd*d)+(i*d)+b))

void correl_res(int xdc, int ydc, int zdc, int xdf, int ydf, int zdf, int GR, int st, int en)
{
	int i,j,k,a;
	FILE *sv;
	FILE *tcp,*wstr;
	char fn[100],fn1[100],fn2[100];
	int num=(en-st)+1;
	DP *z,dplus,Ret;
	DP utauavg=0.,utvt,time,nu;
	int cnt=0,ct=0,rd;
	int up;
	DP ubavg=0.,ubt=0.,twall,val,tau,ub=0.;
	char stri[60];
	DP uttemp,uttemp1,ttemp;
	double zn,zp;
	int nz;
	double *ycorlc,*xcorlc,*ycorltc,*xcorltc;
	double *ycorlfl,*xcorlfl,*ycorltfl,*xcorltfl;
	double *ycorlfu,*xcorlfu,*ycorltfu,*xcorltfu;
	int extx,exty;
	int NZ=2*(zdf-GR-1)+(zdc);
	int NZs=2*zdf+GR*(zdc-4)+GR-1;
	
	xcorlc = (DP *)calloc((xdv*zdc*3),sizeof(DP));
	ycorlc = (DP *)calloc((ydc*zdc*3),sizeof(DP));
	xcorltc = (DP *)calloc((xdv*zdc*3),sizeof(DP));
	ycorltc = (DP *)calloc((ydc*zdc*3),sizeof(DP));
	
	xcorlfl = (DP *)calloc((xdv*GR*zdf*3),sizeof(DP));
	ycorlfl = (DP *)calloc((ydf*zdf*3),sizeof(DP));
	xcorltfl = (DP *)calloc((xdv*GR*zdf*3),sizeof(DP));
	ycorltfl = (DP *)calloc((ydf*zdf*3),sizeof(DP));
	
	xcorlfu = (DP *)calloc((xdv*GR*zdf*3),sizeof(DP));
	ycorlfu = (DP *)calloc((ydf*zdf*3),sizeof(DP));
	xcorltfu = (DP *)calloc((xdv*GR*zdf*3),sizeof(DP));
	ycorltfu = (DP *)calloc((ydf*zdf*3),sizeof(DP));
	
	num/=DT;
	z=(DP *)calloc(NZ,sizeof(DP));

	wstr=fopen("p-grad.txt","r");
	ttemp=TSTR;
	uttemp=0.;
	while (fscanf(wstr,"%lf %lf\n",&ttemp,&uttemp1)!=EOF)
	{
		if ((ttemp>=TSTR)&&(ttemp<=TEND))
		{
			uttemp+=uttemp1;
			ct++;
		}
	}
	fclose(wstr);

	utauavg=sqrt(uttemp/ct);
	tau=(3./Reb)*ubulk*(NZs-2.)+0.5;
	nu=(1./6.)*(2.*tau-1.);
	Ret=utauavg*(NZs-2.)*0.5/nu;
	dplus=Ret/(0.5*(NZs-2));
	printf("ut=%f nu=%f dplus=%f Ret=%f\n",utauavg,nu,dplus,Ret);

	for (k=1;k<(zdf-GR-1);k++)
	  *(z+k) = k-0.5;
	for (k=0;k<(zdc);k++)
	  *(z+zdf-GR-1+k) = GR*(k)+(zdf-GR-1-0.5);
	for (k=GR+1;k<(zdf-1);k++)
	  *(z+zdf+zdc-GR-2-GR+k) = k-GR + GR*(zdc-1) + zdf-GR-1.-0.5;

 	*(z+NZ-1) = 2*zdf+GR*(zdc-4)+GR-3;

        for (k=0;k<NZ;k++)
          *(z+k) /= ((double)(NZs-2.));

        for (k=0;k<NZ;k++)
          *(z+k) *= 2.;

        for (k=0;k<NZ;k++)
          *(z+k) -= 1.;

	zn = -1;
	k=0;
	while (zn<ZPTGT)
	{
		zp = zn;
		zn = ((*(z+k))+1)*Ret;
		k++;
	}
	nz = k-2;
	if (fabs(zn-ZPTGT)<fabs(zp-ZPTGT))
	{
		zp = zn;
		nz++;
		printf("found it\n");
	}
	printf("nz=%d z+=%f\n",nz,zp); 

#ifdef ZOH
	zn=-1;
	k=0;
	while (zn<ZOH)
	{
		zp =zn;
		zn = ((*(z+k))+1);
		k++;
	}
	nz = k-2;
	if (fabs(zn-ZOH)<fabs(zp-ZOH))
	{
		zp = zn;
		nz++;
	}
	printf("forget about the previous one, nz=%d z/h=%f\n",nz,zp); 
#endif
		
	for (k=st;k<(en+1);k+=DT)
	{
		extx = zdf*xdf*3;
		exty = zdf*ydf*3;
		
		sprintf(fn1,"xcorl-fl.%.4d",k);
//		printf("h: %s\n",fn1);
		sv=fopen(fn1,"rb");
		fread(xcorltfl,sizeof(DP),extx,sv);
		fclose(sv);
		
		sprintf(fn2,"ycorl-fl.%.4d",k);
//		printf("h: %s\n",fn2);
		sv=fopen(fn2,"rb");
		fread(ycorltfl,sizeof(DP),exty,sv);
		fclose(sv);

		for (a=0;a<extx;a++)
		{
			*(xcorlfl+a) += (*(xcorltfl+a));
		}
		for (a=0;a<exty;a++)
		{
			*(ycorlfl+a) += (*(ycorltfl+a));
		}
//		printf("a1\n");	
		sprintf(fn1,"xcorl-fu.%.4d",k);
//		printf("h: %s\n",fn1);
		sv=fopen(fn1,"rb");
		fread(xcorltfu,sizeof(DP),extx,sv);
		fclose(sv);
		
		sprintf(fn2,"ycorl-fu.%.4d",k);
//		printf("h: %s\n",fn2);
		sv=fopen(fn2,"rb");
		fread(ycorltfu,sizeof(DP),exty,sv);
		fclose(sv);

		for (a=0;a<extx;a++)
		{
			*(xcorlfu+a) += (*(xcorltfu+a));
		}
		for (a=0;a<exty;a++)
		{
			*(ycorlfu+a) += (*(ycorltfu+a));
		}
//		printf("a2\n");
		
		extx = zdc*xdc*3;
		exty = zdc*ydc*3;
		
		sprintf(fn1,"xcorl-c.%.4d",k);
//		printf("h: %s\n",fn1);
		sv=fopen(fn1,"rb");
		fread(xcorltc,sizeof(DP),extx,sv);
		fclose(sv);
		
		sprintf(fn2,"ycorl-c.%.4d",k);
//		printf("h: %s\n",fn2);
		sv=fopen(fn2,"rb");
		fread(ycorltc,sizeof(DP),exty,sv);
		fclose(sv);

		for (a=0;a<extx;a++)
		{
			*(xcorlc+a) += (*(xcorltc+a));
		}
		for (a=0;a<exty;a++)
		{
			*(ycorlc+a) += (*(ycorltc+a));
		}
	}
//	printf("a3\n");
	extx = zdf*xdf*3;
	exty = zdf*ydf*3;
	for (a=0;a<extx;a++)
	{
		*(xcorlfl+a) /= ((double)num);
	}
	for (a=0;a<exty;a++)
	{
		*(ycorlfl+a) /= ((double)num);
	}
	
	for (a=0;a<extx;a++)
	{
		*(xcorlfu+a) /= ((double)num);
	}
	for (a=0;a<exty;a++)
	{
		*(ycorlfu+a) /= ((double)num);
	}

	extx = zdc*xdc*3;
	exty = zdc*ydc*3;
	for (a=0;a<extx;a++)
	{
		*(xcorlc+a) /= ((double)num);
	}
	for (a=0;a<exty;a++)
	{
		*(ycorlc+a) /= ((double)num);
	}
	
	for (k=0;k<zdf;k++)
	{
		for (j=0;j<ydf;j++)
		{
			Fy(ycorlfl,j,k,0,ydf,zdf,3) = (Fy(ycorlfl,j,k,0,ydf,zdf,3)+Fy(ycorlfu,j,(zdf-1-k),0,ydf,zdf,3))*0.5;
			Fy(ycorlfl,j,k,1,ydf,zdf,3) = (Fy(ycorlfl,j,k,1,ydf,zdf,3)+Fy(ycorlfu,j,(zdf-1-k),1,ydf,zdf,3))*0.5;
			Fy(ycorlfl,j,k,2,ydf,zdf,3) = (Fy(ycorlfl,j,k,2,ydf,zdf,3)+Fy(ycorlfu,j,(zdf-1-k),2,ydf,zdf,3))*0.5;
		}
	}
	for (k=0;k<zdf;k++)
	{
		for (i=0;i<GR*xdv;i++)
		{
			Fx(xcorlfl,i,k,0,xdv*GR,zdf,3) = (Fx(xcorlfl,i,k,0,xdv*GR,zdf,3)+Fx(xcorlfu,i,(zdf-1-k),0,xdv*GR,zdf,3))*0.5;
			Fx(xcorlfl,i,k,1,xdv*GR,zdf,3) = (Fx(xcorlfl,i,k,1,xdv*GR,zdf,3)+Fx(xcorlfu,i,(zdf-1-k),1,xdv*GR,zdf,3))*0.5;
			Fx(xcorlfl,i,k,2,xdv*GR,zdf,3) = (Fx(xcorlfl,i,k,2,xdv*GR,zdf,3)+Fx(xcorlfu,i,(zdf-1-k),2,xdv*GR,zdf,3))*0.5;
		}
	}
	for (k=0;k<zdf;k++)
	{
		for (j=0;j<(ydf+1)/2;j++)
		{
			Fy(ycorlfl,j,k,0,ydf,zdf,3) = (Fy(ycorlfl,j,k,0,ydf,zdf,3)+Fy(ycorlfl,(ydf-1-j),k,0,ydf,zdf,3))*0.5;
			Fy(ycorlfl,j,k,1,ydf,zdf,3) = (Fy(ycorlfl,j,k,1,ydf,zdf,3)+Fy(ycorlfl,(ydf-1-j),k,1,ydf,zdf,3))*0.5;
			Fy(ycorlfl,j,k,2,ydf,zdf,3) = (Fy(ycorlfl,j,k,2,ydf,zdf,3)+Fy(ycorlfl,(ydf-1-j),k,2,ydf,zdf,3))*0.5;
		}
	}
	for (k=0;k<zdf;k++)
	{
		for (i=0;i<(GR*xdv+1)/2;i++)
		{
			Fx(xcorlfl,i,k,0,xdv*GR,zdf,3) = (Fx(xcorlfl,i,k,0,xdv*GR,zdf,3)+Fx(xcorlfl,(xdv*GR-1-i),k,0,xdv*GR,zdf,3))*0.5;
			Fx(xcorlfl,i,k,1,xdv*GR,zdf,3) = (Fx(xcorlfl,i,k,1,xdv*GR,zdf,3)+Fx(xcorlfl,(xdv*GR-1-i),k,1,xdv*GR,zdf,3))*0.5;
			Fx(xcorlfl,i,k,2,xdv*GR,zdf,3) = (Fx(xcorlfl,i,k,2,xdv*GR,zdf,3)+Fx(xcorlfl,(xdv*GR-1-i),k,2,xdv*GR,zdf,3))*0.5;
		}
	}
	
	
	for (k=0;k<(zdc+1)/2;k++)
	{
		for (j=0;j<ydc;j++)
		{
			Fy(ycorlc,j,k,0,ydc,zdc,3) = (Fy(ycorlc,j,k,0,ydc,zdc,3)+Fy(ycorlc,j,(zdc-1-k),0,ydc,zdc,3))*0.5;
			Fy(ycorlc,j,k,1,ydc,zdc,3) = (Fy(ycorlc,j,k,1,ydc,zdc,3)+Fy(ycorlc,j,(zdc-1-k),1,ydc,zdc,3))*0.5;
			Fy(ycorlc,j,k,2,ydc,zdc,3) = (Fy(ycorlc,j,k,2,ydc,zdc,3)+Fy(ycorlc,j,(zdc-1-k),2,ydc,zdc,3))*0.5;
		}
	}
	for (k=0;k<(zdc+1)/2;k++)
	{
		for (i=0;i<xdv;i++)
		{
			Fx(xcorlc,i,k,0,xdv,zdc,3) = (Fx(xcorlc,i,k,0,xdv,zdc,3)+Fx(xcorlc,i,(zdc-1-k),0,xdv,zdc,3))*0.5;
			Fx(xcorlc,i,k,1,xdv,zdc,3) = (Fx(xcorlc,i,k,1,xdv,zdc,3)+Fx(xcorlc,i,(zdc-1-k),1,xdv,zdc,3))*0.5;
			Fx(xcorlc,i,k,2,xdv,zdc,3) = (Fx(xcorlc,i,k,2,xdv,zdc,3)+Fx(xcorlc,i,(zdc-1-k),2,xdv,zdc,3))*0.5;
		}
	}
	for (k=0;k<zdc;k++)
	{
		for (j=0;j<(ydc+1)/2;j++)
		{
			Fy(ycorlc,j,k,0,ydc,zdc,3) = (Fy(ycorlc,j,k,0,ydc,zdc,3)+Fy(ycorlc,(ydc-1-j),k,0,ydc,zdc,3))*0.5;
			Fy(ycorlc,j,k,1,ydc,zdc,3) = (Fy(ycorlc,j,k,1,ydc,zdc,3)+Fy(ycorlc,(ydc-1-j),k,1,ydc,zdc,3))*0.5;
			Fy(ycorlc,j,k,2,ydc,zdc,3) = (Fy(ycorlc,j,k,2,ydc,zdc,3)+Fy(ycorlc,(ydc-1-j),k,2,ydc,zdc,3))*0.5;
		}
	}
	for (k=0;k<zdc;k++)
	{
		for (i=0;i<(xdv+1)/2;i++)
		{
			Fx(xcorlc,i,k,0,xdv,zdc,3) = (Fx(xcorlc,i,k,0,xdv,zdc,3)+Fx(xcorlc,(xdv-1-i),k,0,xdv,zdc,3))*0.5;
			Fx(xcorlc,i,k,1,xdv,zdc,3) = (Fx(xcorlc,i,k,1,xdv,zdc,3)+Fx(xcorlc,(xdv-1-i),k,1,xdv,zdc,3))*0.5;
			Fx(xcorlc,i,k,2,xdv,zdc,3) = (Fx(xcorlc,i,k,2,xdv,zdc,3)+Fx(xcorlc,(xdv-1-i),k,2,xdv,zdc,3))*0.5;
		}
	}
	
	if (nz < (zdf-GR-1) )
	{
			k=nz;
			sprintf(fn,"RX_z%d.dat",((int)ZPTGT));
#ifdef ZOH
			sprintf(fn,"RX_z%.2f.dat",(ZOH));
#endif
			sv=fopen(fn,"w");
			for (i=0;i<(xdv*GR+1)/2;i++)
			{
			  fprintf(sv,"%f ",(i+0.5)*utauavg/nu);
			  fprintf(sv,"%f ",Fx(xcorlfl,i,k,0,xdv*GR,zdf,3)/Fx(xcorlfl,0,k,0,xdv*GR,zdf,3));
			  fprintf(sv,"%f ",Fx(xcorlfl,i,k,1,xdv*GR,zdf,3)/Fx(xcorlfl,0,k,1,xdv*GR,zdf,3));
			  fprintf(sv,"%f\n",Fx(xcorlfl,i,k,2,xdv*GR,zdf,3)/Fx(xcorlfl,0,k,2,xdv*GR,zdf,3));
			}
			fclose(sv);
			
			sprintf(fn,"RY_z%d.dat",((int)ZPTGT));
#ifdef ZOH
			sprintf(fn,"RY_z%.2f.dat",(ZOH));
#endif
			sv=fopen(fn,"w");
			for (j=0;j<(ydf+1)/2;j++)
			{
			  fprintf(sv,"%f ",(j+0.5)*utauavg/nu);
			  fprintf(sv,"%f ",Fy(ycorlfl,j,k,0,ydf,zdf,3)/Fy(ycorlfl,0,k,0,ydf,zdf,3));
			  fprintf(sv,"%f ",Fy(ycorlfl,j,k,1,ydf,zdf,3)/Fy(ycorlfl,0,k,1,ydf,zdf,3));
			  fprintf(sv,"%f\n",Fy(ycorlfl,j,k,2,ydf,zdf,3)/Fy(ycorlfl,0,k,2,ydf,zdf,3));
			}
			fclose(sv);
	}
	else
	{
			k=nz-zdf+GR+1;
			sprintf(fn,"RX_z%d.dat",((int)ZPTGT));
#ifdef ZOH
			sprintf(fn,"RX_z%.2f.dat",(ZOH));
#endif
			sv=fopen(fn,"w");
			for (i=0;i<(xdc+1)/2;i++)
			{
			  fprintf(sv,"%f ",(i+0.5)*utauavg/nu);
			  fprintf(sv,"%f ",Fx(xcorlc,i,k,0,xdv,zdc,3)/Fx(xcorlc,0,k,0,xdv,zdc,3));
			  fprintf(sv,"%f ",Fx(xcorlc,i,k,1,xdv,zdc,3)/Fx(xcorlc,0,k,1,xdv,zdc,3));
			  fprintf(sv,"%f\n",Fx(xcorlc,i,k,2,xdv,zdc,3)/Fx(xcorlc,0,k,2,xdv,zdc,3));
			}
			fclose(sv);
			
			sprintf(fn,"RY_z%d.dat",((int)ZPTGT));
#ifdef ZOH
			sprintf(fn,"RY_z%.2f.dat",(ZOH));
#endif
			sv=fopen(fn,"w");
			for (j=0;j<(ydc+1)/2;j++)
			{
			  fprintf(sv,"%f ",(j+0.5)*utauavg/nu);
			  fprintf(sv,"%f ",Fy(ycorlc,j,k,0,ydc,zdc,3)/Fy(ycorlc,0,k,0,ydc,zdc,3));
			  fprintf(sv,"%f ",Fy(ycorlc,j,k,1,ydc,zdc,3)/Fy(ycorlc,0,k,1,ydc,zdc,3));
			  fprintf(sv,"%f\n",Fy(ycorlc,j,k,2,ydc,zdc,3)/Fy(ycorlc,0,k,2,ydc,zdc,3));
			}
			fclose(sv);
	}
	
	free(ycorlc);
	free(xcorlc);
	free(ycorlfl);
	free(xcorlfl);
	free(ycorlfu);
	free(xcorlfu);

	free(ycorltc);
	free(xcorltc);
	free(ycorltfl);
	free(xcorltfl);
	free(ycorltfu);
	free(xcorltfu);
}	

int main()
{
	int xdc,ydc,zdc;
	int GR, xdf,ydf,zdf;

	zdc=zdv;
	xdc=xdv;
	ydc=ydv;
	
	GR = gratio;
	
	xdf = GR*xdc;
	ydf = GR*ydc;
	zdf = zds+1;
	
	correl_res(xdc,ydc,zdc,xdf,ydf,zdf,GR,TSTR,TEND);
	
	return(0);
}


	
	
	
	
	
	
	
	
	
	
