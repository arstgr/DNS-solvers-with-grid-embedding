#include "definitions.h"

void Vorticity(int xdc, int ydc, int zdc, int xdf, int ydf, int zdf, int GR, int st, int en)
{
	int i,j,k,a,zd;
	FILE *sv;
	FILE *tcp,*wstr;
	DP *mox,*moy,*moz,*fox,*foy,*foz,*moxt,*moyt,*mozt,*foxt,*foyt,*fozt;
	char fn[20];
	int num=(en-st)+1;
	DP *z,dplus,Ret;
	DP utauavg=0.,utvt,time,nu;
	int cnt=0,ct=0,rd;
	int up;
	DP ubavg=0.,ubt=0.,twall,val,tau;
	char stri[30];
	DP uttemp,uttemp1,ttemp;
	DP *du, *dudx;
	int NZ=2*zdf+(zdc-4);
	int NZs=2*zdf+GR*(zdc-4)+GR-1;
	
	zd=NZ;

	up=1+(NZ-1)/2;
	num/=DT;
	z=(DP *)calloc(NZ,sizeof(DP));

	fox = (DP *)calloc(NZ,sizeof(DP));
	foy = (DP *)calloc(NZ,sizeof(DP));
	foz = (DP *)calloc(NZ,sizeof(DP));

	foxt = (DP *)calloc(NZ,sizeof(DP));
	foyt = (DP *)calloc(NZ,sizeof(DP));
	fozt = (DP *)calloc(NZ,sizeof(DP));


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
	printf("vorticity: ut=%f nu=%f Ret=%f dplus=%f\n",utauavg,nu,Ret,dplus);
	
	for (k=1;k<zdf;k++)
          *(z+k) = k-0.5;
        for (k=2;k<(zdc-2);k++)
          *(z+zdf-2+k) = GR*(k-1)+(zdf-1-0.5);
        for (k=0;k<(zdf-1);k++)
          *(z+zdf+zdc-4+k) = k + GR*(zdc-3) + zdf-1.-0.5;
        *(z+2*zdf+zdc-5) = 2*zdf+GR*(zdc-4)+GR-3;

        for (k=0;k<NZ;k++)
          *(z+k) /= ((double)(NZs-2.));

        for (k=0;k<NZ;k++)
          *(z+k) *= 2.;

        for (k=0;k<NZ;k++)
          *(z+k) -= 1.;

	for (k=st;k<(en+1);k+=DT)
	{
		sprintf(fn,"vort-fluc-2d.%.3d.%.4d",0,k);
		sv=fopen(fn,"rb");
		fread(foxt,sizeof(DP),NZ,sv);
		fread(foyt,sizeof(DP),NZ,sv);
		fread(fozt,sizeof(DP),NZ,sv);
		fclose(sv);

		for (a=0;a<NZ;a++)
		{
			*(fox+a) += (*(foxt+a))*(*(foxt+a));
			*(foy+a) += (*(foyt+a))*(*(foyt+a));
			*(foz+a) += (*(fozt+a))*(*(fozt+a));
		}
	}

	for (k=0;k<NZ;k++)
	{
		*(fox+k) /= num;
		*(foy+k) /= num;
		*(foz+k) /= num;
	}

	sprintf(stri,TEXT);
//	Ret=1.;
	
	zd=NZ;

	sprintf(fn,"omega-rms-x-pl.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,sqrt(0.5*((*(fox+k))+(*(fox+zd-1-k))))*nu/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"omega-rms-x-h.dat");
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z/H\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1),sqrt(0.5*((*(fox+k))+(*(fox+zd-1-k))))*nu/(utauavg*utauavg));
	for (k=up-2;k>-1;k--)
                fprintf(tcp,"%f %.9f\n",2.-((*(z+k))+1),sqrt(0.5*((*(fox+k))+(*(fox+zd-1-k))))*nu/(utauavg*utauavg));
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"omega-rms-y-pl.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,sqrt(0.5*((*(foy+k))+(*(foy+zd-1-k))))*nu/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"omega-rms-y-h.dat");
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z/H\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1),sqrt(0.5*((*(foy+k))+(*(foy+zd-1-k))))*nu/(utauavg*utauavg));
	for (k=up-2;k>-1;k--)
                fprintf(tcp,"%f %.9f\n",2.-((*(z+k))+1),sqrt(0.5*((*(foy+k))+(*(foy+zd-1-k))))*nu/(utauavg*utauavg));
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"omega-rms-z-pl.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,sqrt(0.5*((*(foz+k))+(*(foz+zd-1-k))))*nu/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"omega-rms-z-h.dat");
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z/H\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1),sqrt(0.5*((*(foz+k))+(*(foz+zd-1-k))))*nu/(utauavg*utauavg));
	for (k=up-2;k>-1;k--)
                fprintf(tcp,"%f %.9f\n",2.-((*(z+k))+1),sqrt(0.5*((*(foz+k))+(*(foz+zd-1-k))))*nu/(utauavg*utauavg));
        fprintf(tcp,"\n");
        fclose(tcp);

	free(fox);
	free(foy);
	free(foz);

	free(foxt);
	free(foyt);
	free(fozt);
}
