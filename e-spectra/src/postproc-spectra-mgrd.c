#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* postporcessing the results of energy spectra code */
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

#define Reb 7200.0000
#define ubulk 0.0866034
#define TSTR 916
#define TEND 1216
#define DT 1
#define ZPTGT 110.

#define ALPHA 1.35
#define BETA 2.5

typedef double DP;
#define PI (4.*atan(1.))

void e_spact(int xdc, int ydc, int zdc, int xdf, int ydf, int zdf, int GR, int st, int en)
{
	int i,j,k,a,b;
	
	double Euuxc[xdc][zdc],Evvxc[xdc][zdc],Ewwxc[xdc][zdc];
	double Euuxtc[xdc][zdc],Evvxtc[xdc][zdc],Ewwxtc[xdc][zdc];
	double Euuxec[xdc/2],Evvxec[xdc/2],Ewwxec[xdc/2];
	
	double Euuxfl[xdf][zdf],Evvxfl[xdf][zdf],Ewwxfl[xdf][zdf];
	double Euuxtfl[xdf][zdf],Evvxtfl[xdf][zdf],Ewwxtfl[xdf][zdf];
	
	double Euuxfu[xdf][zdf],Evvxfu[xdf][zdf],Ewwxfu[xdf][zdf];
	double Euuxtfu[xdf][zdf],Evvxtfu[xdf][zdf],Ewwxtfu[xdf][zdf];
	double Euuxef[xdf/2],Evvxef[xdf/2],Ewwxef[xdf/2];
	
	double Euuyc[ydc/2 +1][zdc],Evvyc[ydc/2 +1][zdc],Ewwyc[ydc/2 +1][zdc];
	double Euuytc[ydc/2 +1][zdc],Evvytc[ydc/2 +1][zdc],Ewwytc[ydc/2 +1][zdc];
	double Euuyec[ydc/2],Evvyec[ydc/2],Ewwyec[ydc/2];
	
	double Euuyfl[ydf/2 +1][zdf],Evvyfl[ydf/2 +1][zdf],Ewwyfl[ydf/2 +1][zdf];
	double Euuytfl[ydf/2 +1][zdf],Evvytfl[ydf/2 +1][zdf],Ewwytfl[ydf/2 +1][zdf];
	
	double Euuyfu[ydf/2 +1][zdf],Evvyfu[ydf/2 +1][zdf],Ewwyfu[ydf/2 +1][zdf];
	double Euuytfu[ydf/2 +1][zdf],Evvytfu[ydf/2 +1][zdf],Ewwytfu[ydf/2 +1][zdf];
	double Euuyef[ydf/2],Evvyef[ydf/2],Ewwyef[ydf/2];
	
	double Kxc[xdc/2],Kyc[ydc/2];
	double Kxf[xdf/2],Kyf[ydf/2];
	int NZ=2*(zdf-GR-1)+(zdc);
	int NZs=2*zdf+GR*(zdc-4)+GR-1;
	int nz;
	double z[NZ],zp,zn;
	
	double fac;
	
	int num=(en-st)+1;
	int up=1+(NZs-1)/2;
	
	double nu,tau,utauavg,uttemp,uttemp1,ttemp,dplus,Ret;
	int ct=0;
	FILE *wstr,*fh;
	char fn[80];
	int nxs,nys;
	
	up=1+(NZs-1)/2;
	num /= DT;
	
	for (i=0;i<xdc/2;i++)
	  Kxc[i] = i*ALPHA*GR/(0.5*(NZs-2.));
	
	for (j=0;j<ydc/2;j++)
	  Kyc[j] = j*BETA*GR/(0.5*(NZs-2.));
	
	
	for (i=0;i<xdf/2;i++)
	  Kxf[i] = i*ALPHA/(0.5*(NZs-2.));
	
	for (j=0;j<ydf/2;j++)
	  Kyf[j] = j*BETA/(0.5*(NZs-2.));
	
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
	printf("ut=%f nu=%f  Ret=%f\n",utauavg,nu,Ret);
	
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
	
	/* Initialization */
	printf("initialization\n");
	for (a=0;a<xdc;a++)
	{
		for (b=0;b<zdc;b++)
		{
			Euuxc[a][b] = 0.;
			Evvxc[a][b] = 0.;
			Ewwxc[a][b] = 0.;
		}
	}
	for (a=0;a<(ydc/2 +1);a++)
	{
		for (b=0;b<zdc;b++)
		{
			Euuyc[a][b] = 0.;
			Evvyc[a][b] = 0.;
			Ewwyc[a][b] = 0.;
		}
	}
	for (a=0;a<xdc/2;a++)
	{
		Euuxec[a] = 0.;
		Evvxec[a] = 0.;
		Ewwxec[a] = 0.;
	}
	for (a=0;a<ydc/2;a++)
	{
		Euuyec[a] = 0.;
		Evvyec[a] = 0.;
		Ewwyec[a] = 0.;
	}
	
	for (a=0;a<xdf;a++)
	{
		for (b=0;b<zdf;b++)
		{
			Euuxfl[a][b] = 0.;
			Evvxfl[a][b] = 0.;
			Ewwxfl[a][b] = 0.;
		}
	}
	for (a=0;a<(ydf/2 +1);a++)
	{
		for (b=0;b<zdf;b++)
		{
			Euuyfl[a][b] = 0.;
			Evvyfl[a][b] = 0.;
			Ewwyfl[a][b] = 0.;
		}
	}
	for (a=0;a<xdf/2;a++)
	{
		Euuxef[a] = 0.;
		Evvxef[a] = 0.;
		Ewwxef[a] = 0.;
	}
	for (a=0;a<ydf/2;a++)
	{
		Euuyef[a] = 0.;
		Evvyef[a] = 0.;
		Ewwyef[a] = 0.;
	}
	
	for (a=0;a<xdf;a++)
	{
		for (b=0;b<zdf;b++)
		{
			Euuxfu[a][b] = 0.;
			Evvxfu[a][b] = 0.;
			Ewwxfu[a][b] = 0.;
		}
	}
	for (a=0;a<(ydf/2 +1);a++)
	{
		for (b=0;b<zdf;b++)
		{
			Euuyfu[a][b] = 0.;
			Evvyfu[a][b] = 0.;
			Ewwyfu[a][b] = 0.;
		}
	}
	
	printf("reading\n");
	/* nz is the z plane of interest */
	for (k=st;k<(en+1);k+=DT)
	{
		sprintf(fn,"Euux-c.%.4d",k);
		fh=fopen(fn,"rb");
		fread(&Euuxtc[0][0],sizeof(DP),(zdc*xdc),fh);
		fclose(fh);
		
		sprintf(fn,"Evvx-c.%.4d",k);
		fh=fopen(fn,"rb");
		fread(&Evvxtc[0][0],sizeof(DP),(zdc*xdc),fh);
		fclose(fh);
		
		sprintf(fn,"Ewwx-c.%.4d",k);
		fh=fopen(fn,"rb");
		fread(&Ewwxtc[0][0],sizeof(DP),(zdc*xdc),fh);
		fclose(fh);
		
		sprintf(fn,"Euuy-c.%.4d",k);
		fh=fopen(fn,"rb");
		fread(&Euuytc[0][0],sizeof(DP),(zdc*(ydc/2 +1)),fh);
		fclose(fh);
		
		sprintf(fn,"Evvy-c.%.4d",k);
		fh=fopen(fn,"rb");
		fread(&Evvytc[0][0],sizeof(DP),(zdc*(ydc/2 +1)),fh);
		fclose(fh);
		
		sprintf(fn,"Ewwy-c.%.4d",k);
		fh=fopen(fn,"rb");
		fread(&Ewwytc[0][0],sizeof(DP),(zdc*(ydc/2 +1)),fh);
		fclose(fh);
		/********************************/
		sprintf(fn,"Euux-fl.%.4d",k);
		fh=fopen(fn,"rb");
		fread(&Euuxtfl[0][0],sizeof(DP),(zdf*xdf),fh);
		fclose(fh);
		
		sprintf(fn,"Evvx-fl.%.4d",k);
		fh=fopen(fn,"rb");
		fread(&Evvxtfl[0][0],sizeof(DP),(zdf*xdf),fh);
		fclose(fh);
		
		sprintf(fn,"Ewwx-fl.%.4d",k);
		fh=fopen(fn,"rb");
		fread(&Ewwxtfl[0][0],sizeof(DP),(zdf*xdf),fh);
		fclose(fh);
		
		sprintf(fn,"Euuy-fl.%.4d",k);
		fh=fopen(fn,"rb");
		fread(&Euuytfl[0][0],sizeof(DP),(zdf*(ydf/2 +1)),fh);
		fclose(fh);
		
		sprintf(fn,"Evvy-fl.%.4d",k);
		fh=fopen(fn,"rb");
		fread(&Evvytfl[0][0],sizeof(DP),(zdf*(ydf/2 +1)),fh);
		fclose(fh);
		
		sprintf(fn,"Ewwy-fl.%.4d",k);
		fh=fopen(fn,"rb");
		fread(&Ewwytfl[0][0],sizeof(DP),(zdf*(ydf/2 +1)),fh);
		fclose(fh);
		/*******************************/
		sprintf(fn,"Euux-fu.%.4d",k);
		fh=fopen(fn,"rb");
		fread(&Euuxtfu[0][0],sizeof(DP),(zdf*xdf),fh);
		fclose(fh);
		
		sprintf(fn,"Evvx-fu.%.4d",k);
		fh=fopen(fn,"rb");
		fread(&Evvxtfu[0][0],sizeof(DP),(zdf*xdf),fh);
		fclose(fh);
		
		sprintf(fn,"Ewwx-fu.%.4d",k);
		fh=fopen(fn,"rb");
		fread(&Ewwxtfu[0][0],sizeof(DP),(zdf*xdf),fh);
		fclose(fh);
		
		sprintf(fn,"Euuy-fu.%.4d",k);
		fh=fopen(fn,"rb");
		fread(&Euuytfu[0][0],sizeof(DP),(zdf*(ydf/2 +1)),fh);
		fclose(fh);
		
		sprintf(fn,"Evvy-fu.%.4d",k);
		fh=fopen(fn,"rb");
		fread(&Evvytfu[0][0],sizeof(DP),(zdf*(ydf/2 +1)),fh);
		fclose(fh);
		
		sprintf(fn,"Ewwy-fu.%.4d",k);
		fh=fopen(fn,"rb");
		fread(&Ewwytfu[0][0],sizeof(DP),(zdf*(ydf/2 +1)),fh);
		fclose(fh);
//		printf("adding\n");
		for (a=0;a<xdc;a++)
		{
			for (b=0;b<zdc;b++)
			{
				Euuxc[a][b] += Euuxtc[a][b];
				Evvxc[a][b] += Evvxtc[a][b];
				Ewwxc[a][b] += Ewwxtc[a][b];
			}
		}
		for (a=0;a<(ydc/2 +1);a++)
		{
			for (b=0;b<zdc;b++)
			{
				Euuyc[a][b] += Euuytc[a][b];
				Evvyc[a][b] += Evvytc[a][b];
				Ewwyc[a][b] += Ewwytc[a][b];
			}
		}
		
		for (a=0;a<xdf;a++)
		{
			for (b=0;b<zdf;b++)
			{
				Euuxfl[a][b] += Euuxtfl[a][b];
				Evvxfl[a][b] += Evvxtfl[a][b];
				Ewwxfl[a][b] += Ewwxtfl[a][b];
			}
		}
		for (a=0;a<(ydf/2 +1);a++)
		{
			for (b=0;b<zdf;b++)
			{
				Euuyfl[a][b] += Euuytfl[a][b];
				Evvyfl[a][b] += Evvytfl[a][b];
				Ewwyfl[a][b] += Ewwytfl[a][b];
			}
		}
		
		for (a=0;a<xdf;a++)
		{
			for (b=0;b<zdf;b++)
			{
				Euuxfu[a][b] += Euuxtfu[a][b];
				Evvxfu[a][b] += Evvxtfu[a][b];
				Ewwxfu[a][b] += Ewwxtfu[a][b];
			}
		}
		for (a=0;a<(ydf/2 +1);a++)
		{
			for (b=0;b<zdf;b++)
			{
				Euuyfu[a][b] += Euuytfu[a][b];
				Evvyfu[a][b] += Evvytfu[a][b];
				Ewwyfu[a][b] += Ewwytfu[a][b];
			}
		}
	}
	printf("end of reading \n");
	for (a=0;a<xdc;a++)
	{
		for (b=0;b<zdc;b++)
		{
			Euuxc[a][b] /= num;
			Evvxc[a][b] /= num;
			Ewwxc[a][b] /= num;
		}
	}
	for (a=0;a<(ydc/2 +1);a++)
	{
		for (b=0;b<zdc;b++)
		{
			Euuyc[a][b] /= num;
			Evvyc[a][b] /= num;
			Ewwyc[a][b] /= num;
		}
	}
	
	for (a=0;a<xdf;a++)
	{
		for (b=0;b<zdf;b++)
		{
			Euuxfl[a][b] /= num;
			Evvxfl[a][b] /= num;
			Ewwxfl[a][b] /= num;
		}
	}
	for (a=0;a<(ydf/2 +1);a++)
	{
		for (b=0;b<zdf;b++)
		{
			Euuyfl[a][b] /= num;
			Evvyfl[a][b] /= num;
			Ewwyfl[a][b] /= num;
		}
	}
	
	for (a=0;a<xdf;a++)
	{
		for (b=0;b<zdf;b++)
		{
			Euuxfu[a][b] /= num;
			Evvxfu[a][b] /= num;
			Ewwxfu[a][b] /= num;
		}
	}
	for (a=0;a<(ydf/2 +1);a++)
	{
		for (b=0;b<zdf;b++)
		{
			Euuyfu[a][b] /= num;
			Evvyfu[a][b] /= num;
			Ewwyfu[a][b] /= num;
		}
	}
	
	/* X spectra */
	if (nz < (zdf-GR-1) )
	{
		k=nz;
		a=0;
		Euuxef[a] = 0.5*(Euuxfl[a][k] + Euuxfu[a][zdf-1-k]);
		Evvxef[a] = 0.5*(Evvxfl[a][k] + Evvxfu[a][zdf-1-k]);
		Ewwxef[a] = 0.5*(Ewwxfl[a][k] + Ewwxfu[a][zdf-1-k]);

		for (a=1;a<xdf/2;a++)
		{
			Euuxef[a] = 0.5*(Euuxfl[a][k]+Euuxfu[a][zdf-1-k]) + 0.5*(Euuxfl[xdf-a][k]+Euuxfu[xdf-a][zdf-1-k]);
			Evvxef[a] = 0.5*(Evvxfl[a][k]+Evvxfu[a][zdf-1-k]) + 0.5*(Evvxfl[xdf-a][k]+Evvxfu[xdf-a][zdf-1-k]);
			Ewwxef[a] = 0.5*(Ewwxfl[a][k]+Ewwxfu[a][zdf-1-k]) + 0.5*(Ewwxfl[xdf-a][k]+Ewwxfu[xdf-a][zdf-1-k]);
		}
	  
	}
	else if ((nz<(zdf-GR+zdc-1))&&(nz>=(zdf-GR-1)))
	{
		k=nz-zdf+GR+1;
		a=0;
		Euuxec[a] = 0.5*(Euuxc[a][k] + Euuxc[a][zdc-1-k]);
		Evvxec[a] = 0.5*(Evvxc[a][k] + Evvxc[a][zdc-1-k]);
		Ewwxec[a] = 0.5*(Ewwxc[a][k] + Ewwxc[a][zdc-1-k]);

		for (a=1;a<xdc/2;a++)
		{
			Euuxec[a] = 0.5*(Euuxc[a][k]+Euuxc[a][zdc-1-k]) + 0.5*(Euuxc[xdc-a][k]+Euuxc[xdc-a][zdc-1-k]);
			Evvxec[a] = 0.5*(Evvxc[a][k]+Evvxc[a][zdc-1-k]) + 0.5*(Evvxc[xdc-a][k]+Evvxc[xdc-a][zdc-1-k]);
			Ewwxec[a] = 0.5*(Ewwxc[a][k]+Ewwxc[a][zdc-1-k]) + 0.5*(Ewwxc[xdc-a][k]+Ewwxc[xdc-a][zdc-1-k]);
		}
	}
	else if (nz>(zdf-GR-1+zdc-1))
	{
		k=nz-zdf-zdc+2*GR+2;
		a=0;
		Euuxef[a] = 0.5*(Euuxfu[a][k] + Euuxfl[a][zdf-1-k]);
		Evvxef[a] = 0.5*(Evvxfu[a][k] + Evvxfl[a][zdf-1-k]);
		Ewwxef[a] = 0.5*(Ewwxfu[a][k] + Ewwxfl[a][zdf-1-k]);

		for (a=1;a<xdf/2;a++)
		{
			Euuxef[a] = 0.5*(Euuxfu[a][k]+Euuxfl[a][zdf-1-k]) + 0.5*(Euuxfu[xdf-a][k]+Euuxfl[xdf-a][zdf-1-k]);
			Evvxef[a] = 0.5*(Evvxfu[a][k]+Evvxfl[a][zdf-1-k]) + 0.5*(Evvxfu[xdf-a][k]+Evvxfl[xdf-a][zdf-1-k]);
			Ewwxef[a] = 0.5*(Ewwxfu[a][k]+Ewwxfl[a][zdf-1-k]) + 0.5*(Ewwxfu[xdf-a][k]+Ewwxfl[xdf-a][zdf-1-k]);
		}
	}
	
	/* Y spectra */
	if (nz < (zdf-GR-1) )
	{
		k=nz;
		for (a=0;a<ydf/2;a++)
		{
			Euuyef[a] = 0.5*(Euuyfl[a][k]+Euuyfu[a][zdf-1-k]);
			Evvyef[a] = 0.5*(Evvyfl[a][k]+Evvyfu[a][zdf-1-k]);
			Ewwyef[a] = 0.5*(Ewwyfl[a][k]+Ewwyfu[a][zdf-1-k]);
		}
	}
	else if ((nz<(zdf-GR+zdc-1))&&(nz>=(zdf-GR-1)))
	{
		k=nz-zdf+GR+1;
		for (a=0;a<ydc/2;a++)
		{
			Euuyec[a] = 0.5*(Euuyc[a][k]+Euuyc[a][zdc-1-k]);
			Evvyec[a] = 0.5*(Evvyc[a][k]+Evvyc[a][zdc-1-k]);
			Ewwyec[a] = 0.5*(Ewwyc[a][k]+Ewwyc[a][zdc-1-k]);
		}
	}
	else if (nz>(zdf-GR-1+zdc-1))
	{
		k=nz-zdf-zdc+2*GR+2;
		for (a=0;a<ydf/2;a++)
		{
			Euuyef[a] = 0.5*(Euuyfu[a][k]+Euuyfl[a][zdf-1-k]);
			Evvyef[a] = 0.5*(Evvyfu[a][k]+Evvyfl[a][zdf-1-k]);
			Ewwyef[a] = 0.5*(Ewwyfu[a][k]+Ewwyfl[a][zdf-1-k]);
		}
	}

	/* Dealiasing */
	nxs = 2*(xdc/2)/3;
	nys = 2*((ydc/2)/3);
	for (a=1;a<xdc/2;a++)
	{
		if (a > nxs)
		{
			Euuxec[a] = 0.;
			Evvxec[a] = 0.;
			Ewwxec[a] = 0.;
		}
	}
	for (a=0;a<ydc/2;a++)
	{
		if (a > nys)
		{
			Euuyec[a] = 0.;
			Evvyec[a] = 0.;
			Ewwyec[a] = 0.;
		}
	}
	nxs = 2*(xdf/2)/3;
	nys = 2*((ydf/2)/3);
	for (a=1;a<xdf/2;a++)
	{
		if (a > nxs)
		{
			Euuxef[a] = 0.;
			Evvxef[a] = 0.;
			Ewwxef[a] = 0.;
		}
	}
	for (a=0;a<ydf/2;a++)
	{
		if (a > nys)
		{
			Euuyef[a] = 0.;
			Evvyef[a] = 0.;
			Ewwyef[a] = 0.;
		}
	}
	
	if (nz < (zdf-GR-1) )
	{
		k=nz;
		sprintf(fn,"EX-zpl%d.dat",(int)(ZPTGT));
		fh=fopen(fn,"w");
		fprintf(fh,"# Kxplus Euuxplus Evvxpl Ewwxpl\n");
		for (a=0;a<xdf/2;a++)
			fprintf(fh,"%.14f %.14f %.14f %.14f\n",Kxf[a]*nu/utauavg,Euuxef[a]*0.5*(NZs-2.)/(ALPHA*nu*utauavg),Evvxef[a]*0.5*(NZs-2.)/(ALPHA*nu*utauavg),Ewwxef[a]*0.5*(NZs-2.)/(ALPHA*nu*utauavg));
		fclose(fh);
		
		sprintf(fn,"EY-zpl%d.dat",(int)(ZPTGT));
		fh=fopen(fn,"w");
		fprintf(fh,"# Kyplus Euuyplus Evvypl Ewwypl\n");
		for (a=0;a<ydf/2;a++)
			fprintf(fh,"%.14f %.14f %.14f %.14f\n",Kyf[a]*nu/utauavg,Euuyef[a]*0.5*(NZs-2.)/(BETA*nu*utauavg),Evvyef[a]*0.5*(NZs-2.)/(BETA*nu*utauavg),Ewwyef[a]*0.5*(NZs-2.)/(BETA*nu*utauavg));
		fclose(fh);
	}
	else if ((nz<(zdf-GR+zdc-1))&&(nz>=(zdf-GR-1)))
	{
		k=nz-zdf+GR+1;
		sprintf(fn,"EX-zpl%d.dat",(int)(ZPTGT));
		fh=fopen(fn,"w");
		fprintf(fh,"# Kxplus Euuxplus Evvxpl Ewwxpl\n");
		for (a=0;a<xdc/2;a++)
			fprintf(fh,"%.14f %.14f %.14f %.14f\n",Kxc[a]*nu/utauavg,Euuxec[a]*0.5*(NZs-2.)/(ALPHA*nu*utauavg),Evvxec[a]*0.5*(NZs-2.)/(ALPHA*nu*utauavg),Ewwxec[a]*0.5*(NZs-2.)/(ALPHA*nu*utauavg));
		fclose(fh);
		
		sprintf(fn,"EY-zpl%d.dat",(int)(ZPTGT));
		fh=fopen(fn,"w");
		fprintf(fh,"# Kyplus Euuyplus Evvypl Ewwypl\n");
		for (a=0;a<ydc/2;a++)
			fprintf(fh,"%.14f %.14f %.14f %.14f\n",Kyc[a]*nu/utauavg,Euuyec[a]*0.5*(NZs-2.)/(BETA*nu*utauavg),Evvyec[a]*0.5*(NZs-2.)/(BETA*nu*utauavg),Ewwyec[a]*0.5*(NZs-2.)/(BETA*nu*utauavg));
		fclose(fh);
	}
	else if (nz>(zdf-GR-1+zdc-1))
	{
		k=nz-zdf-zdc+2*GR+2;
		sprintf(fn,"EX-zpl%d.dat",(int)(ZPTGT));
		fh=fopen(fn,"w");
		fprintf(fh,"# Kxplus Euuxplus Evvxpl Ewwxpl\n");
		for (a=0;a<xdf/2;a++)
			fprintf(fh,"%.14f %.14f %.14f %.14f\n",Kxf[a]*nu/utauavg,Euuxef[a]*0.5*(NZs-2.)/(ALPHA*nu*utauavg),Evvxef[a]*0.5*(NZs-2.)/(ALPHA*nu*utauavg),Ewwxef[a]*0.5*(NZs-2.)/(ALPHA*nu*utauavg));
		fclose(fh);
		
		sprintf(fn,"EY-zpl%d.dat",(int)(ZPTGT));
		fh=fopen(fn,"w");
		fprintf(fh,"# Kyplus Euuyplus Evvypl Ewwypl\n");
		for (a=0;a<ydf/2;a++)
			fprintf(fh,"%.14f %.14f %.14f %.14f\n",Kyf[a]*nu/utauavg,Euuyef[a]*0.5*(NZs-2.)/(BETA*nu*utauavg),Evvyef[a]*0.5*(NZs-2.)/(BETA*nu*utauavg),Ewwyef[a]*0.5*(NZs-2.)/(BETA*nu*utauavg));
		fclose(fh);
	}
}

int main()
{
	int i,j,k;
	int xdc,ydc,zdc;
	int xdf,ydf,zdf;
	int GR;
	
	GR=gratio;
	
	xdc = xdv;
	ydc = ydv;
	zdc = zdv;
	
	xdf = GR*xdc;
	ydf = GR*ydc;
	zdf = zds+1;
	
	e_spact(xdc, ydc, zdc, xdf, ydf, zdf, GR, TSTR, TEND);
  
  return 0;
}
