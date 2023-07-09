#include "definitions.h"

void time_stat(int xdc, int ydc, int zdc, int xdf, int ydf, int zdf, int GR, int st, int en)
{
	int i,j,k,a;
	FILE *sv;
	FILE *tcp,*wstr;
	DP *uavg,*vavg,*wavg,*davg;
	DP *utin,*vtin,*wtin,*reys,*dtin;
	DP *uavgt,*vavgt,*wavgt,*davgt;
	DP *utint,*vtint,*wtint,*reyst,*dtint;
	DP *sku,*skv,*skw,*skut,*skvt,*skwt;
	DP *skuw,*skp,*skuwt,*skpt;
	DP *kuu,*kuv,*kuw,*kuut,*kuvt,*kuwt;
	DP *kuuw,*kup,*kuuwt,*kupt;
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
	int zd;
	
	zd=NZ;

	up=1+(zd-1)/2;
	num/=DT;
	z=(DP *)calloc(zd,sizeof(DP));

	uavg = (DP *)calloc(zd,sizeof(DP));
	vavg = (DP *)calloc(zd,sizeof(DP));
	wavg = (DP *)calloc(zd,sizeof(DP));
	davg = (DP *)calloc(zd,sizeof(DP));

	utin = (DP *)calloc(zd,sizeof(DP));
	vtin = (DP *)calloc(zd,sizeof(DP));
	wtin = (DP *)calloc(zd,sizeof(DP));
	reys = (DP *)calloc(zd,sizeof(DP));
	dtin = (DP *)calloc(zd,sizeof(DP));

	uavgt = (DP *)calloc(zd,sizeof(DP));
	vavgt = (DP *)calloc(zd,sizeof(DP));
	wavgt = (DP *)calloc(zd,sizeof(DP));
	davgt = (DP *)calloc(zd,sizeof(DP));

	utint = (DP *)calloc(zd,sizeof(DP));
	vtint = (DP *)calloc(zd,sizeof(DP));
	wtint = (DP *)calloc(zd,sizeof(DP));
	reyst = (DP *)calloc(zd,sizeof(DP));
	dtint = (DP *)calloc(zd,sizeof(DP));

	dudx = (DP *)calloc(up,sizeof(DP));
        du = (DP*)calloc(up,sizeof(DP));
        
        sku = (DP *)calloc(zd,sizeof(DP));
        skv = (DP *)calloc(zd,sizeof(DP));
        skw = (DP *)calloc(zd,sizeof(DP));
        skuw = (DP *)calloc(zd,sizeof(DP));
        skut = (DP *)calloc(zd,sizeof(DP));
        skvt = (DP *)calloc(zd,sizeof(DP));
        skwt = (DP *)calloc(zd,sizeof(DP));
        skuwt = (DP *)calloc(zd,sizeof(DP));
        
        skp = (DP *)calloc(zd,sizeof(DP));
        skpt = (DP *)calloc(zd,sizeof(DP));
        
        kuu = (DP *)calloc(zd,sizeof(DP));
        kuv = (DP *)calloc(zd,sizeof(DP));
        kuw = (DP *)calloc(zd,sizeof(DP));
        kuuw = (DP *)calloc(zd,sizeof(DP));
        kuut = (DP *)calloc(zd,sizeof(DP));
        kuvt = (DP *)calloc(zd,sizeof(DP));
        kuwt = (DP *)calloc(zd,sizeof(DP));
        kuuwt = (DP *)calloc(zd,sizeof(DP));
        
        kup = (DP *)calloc(zd,sizeof(DP));
        kupt = (DP *)calloc(zd,sizeof(DP));

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
	printf("ut=%f nu=%f\n",utauavg,nu);

	sv=fopen("flowrate.txt","r");
        ct=0;
        while (fscanf(sv,"%lf %lf\n",&ttemp,&ubt)!=EOF)
        {
                if ((ttemp>=TSTR)&&(ttemp<=TEND))
                {
                        ubavg+=ubt;
                        ct++;
                }
        }
        fclose(sv);
	ubavg=ubavg/((double)(ct*ydf*(NZs-2.)));
        printf("ub=%f\n",ubavg);
	
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
		sprintf(fn,"turb-fld-2d.%.3d.%.4d",0,k);
		sv=fopen(fn,"rb");
		fread(uavgt,sizeof(DP),zd,sv);
		fread(vavgt,sizeof(DP),zd,sv);
		fread(wavgt,sizeof(DP),zd,sv);
		fread(utint,sizeof(DP),zd,sv);
		fread(vtint,sizeof(DP),zd,sv);
		fread(wtint,sizeof(DP),zd,sv);
		fread(reyst,sizeof(DP),zd,sv);
		fread(davgt,sizeof(DP),zd,sv);
		fread(dtint,sizeof(DP),zd,sv);
		fread(skut,sizeof(DP),zd,sv);
		fread(skvt,sizeof(DP),zd,sv);
		fread(skwt,sizeof(DP),zd,sv);
		fread(skuwt,sizeof(DP),zd,sv);
		fread(skpt,sizeof(DP),zd,sv);
		fread(kuut,sizeof(DP),zd,sv);
		fread(kuvt,sizeof(DP),zd,sv);
		fread(kuwt,sizeof(DP),zd,sv);
		fread(kuuwt,sizeof(DP),zd,sv);
		fread(kupt,sizeof(DP),zd,sv);
		fclose(sv);

		for (a=0;a<zd;a++)
		{
			*(uavg+a) += *(uavgt+a);
			*(vavg+a) += *(vavgt+a);
			*(wavg+a) += *(wavgt+a);
			*(davg+a) += *(davgt+a);

			*(utin+a) += (*(utint+a))*(*(utint+a));
			*(vtin+a) += (*(vtint+a))*(*(vtint+a));
			*(wtin+a) += (*(wtint+a))*(*(wtint+a));
			*(dtin+a) += (*(dtint+a))*(*(dtint+a));
			*(reys+a) += (*(reyst+a));
			*(sku+a) += (*(skut+a));
			*(skv+a) += (*(skvt+a));
			*(skw+a) += (*(skwt+a));
			*(skuw+a) += (*(skuwt+a));
			*(skp+a) += (*(skpt+a));
			*(kuu+a) += (*(kuut+a));
			*(kuv+a) += (*(kuvt+a));
			*(kuw+a) += (*(kuwt+a));
			*(kuuw+a) += (*(kuuwt+a));
			*(kup+a) += (*(kupt+a));
		}
	}

	for (k=0;k<zd;k++)
	{
		*(uavg+k) /= num;
		*(vavg+k) /= num;
		*(wavg+k) /= num;
		*(davg+k) /= num;

		*(utin+k) /= num;
		*(vtin+k) /= num;
		*(wtin+k) /= num;
		*(reys+k) /= num;
		*(dtin+k) /= num;
		
		*(sku+k) /= num;
		*(skv+k) /= num;
		*(skw+k) /= num;
		*(skuw+k) /= num;
		*(skp+k) /= num;
		*(kuu+k) /= num;
		*(kuv+k) /= num;
		*(kuw+k) /= num;
		*(kuuw+k) /= num;
		*(kup+k) /= num;
	}

	for (k=0;k<up;k++)
                *(du+k) = 0.5*((*(uavg+k))+(*(uavg+zd-1-k)));

        *(dudx+0)=(-2.*(*(du+1))+3.*(*(du+2))-(*(du+3)));//((*(z+2)+1.)-(*(z+1)+1.));
        printf("ut=%f tw=%f\n",sqrt(nu*(*(dudx+0))),nu*(*(dudx+0)));
        *(dudx+1)=(-3.*(*(du+1))+4.*(*(du+2))-(*(du+3)))/2.;///(2.*((*(z+2)+1.)-(*(z+1)+1.)));
//        for (k=2;k<(up-1);k++)
//                *(dudx+k)=(*(du+k+1)-*(du+k-1))/2.;//(*(z+k+1)-*(z+k-1));
        for (k=2;k<(up-2);k++)
                *(dudx+k)=(-3.*(*(du+k))+4.*(*(du+k+1))-(*(du+k+2)))/2.;
        *(dudx+up-2)=(*(du+up-2)-*(du+up-3));
        *(dudx+up-1)=(*(du+up-1)-*(du+up-2));//(*(z+up-1)-*(z+up-2));


	sprintf(stri,TEXT);
//	Ret=1.;

	sprintf(fn,"Uavg.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"U<sup>+</sup>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(uavg+k))+(*(uavg+zd-1-k)))/utauavg);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Uavg-Ubulk.dat");
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z/H\", \"U/U<sub>bulk</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
        for (k=0;k<zd;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(uavg+k))/ubavg);
        fprintf(tcp,"\n");
        fclose(tcp);


	sprintf(fn,"Us-b.dat");
        tcp=fopen(fn,"w");
//        fprintf(tcp,"TITLE= \"DNS Results\"\n");
//        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"U<sup>+</sup>\"\n");
//        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
        val=((0.5*((*(uavg+1))+(*(uavg+zd-2))))/utauavg)-(*(z+1)+1.)*Ret;
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(uavg+k))+(*(uavg+zd-1-k)))/utauavg-val);
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"Vavg.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"V<sup>+</sup>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(vavg+k))+(*(vavg+zd-1-k)))/utauavg);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Wavg.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"W<sup>+</sup>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(wavg+k))+(*(wavg+zd-1-k)))/utauavg);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Davg.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<greek>r</greek><sup>+</sup>\"\n");
	fprintf(tcp,"ZONE T=\"D AVG\", I=%d, F=POINT\n",up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(davg+k))+(*(davg+zd-1-k)))/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"utin.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"u<sup>+</sup><sub>rms</sub>,v<sup>+</sup><sub>rms</sub>,w<sup>+</sup><sub>rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,sqrt(0.5*((*(utin+k))+(*(utin+zd-1-k))))/utauavg);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"utin-h.dat");
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z/H\", \"u<sup>+</sup><sub>rms</sub>,v<sup>+</sup><sub>rms</sub>,w<sup>+</sup><sub>rms</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1),sqrt(0.5*((*(utin+k))+(*(utin+zd-1-k))))/utauavg);
	for (k=up-2;k>-1;k--)
                fprintf(tcp,"%f %.9f\n",2.-((*(z+k))+1),sqrt(0.5*((*(utin+k))+(*(utin+zd-1-k))))/utauavg);
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"vtin.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"u<sup>+</sup><sub>rms</sub>,v<sup>+</sup><sub>rms</sub>,w<sup>+</sup><sub>rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,sqrt(0.5*((*(vtin+k))+(*(vtin+zd-1-k))))/utauavg);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"vtin-h.dat");
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z/H\", \"u<sup>+</sup><sub>rms</sub>,v<sup>+</sup><sub>rms</sub>,w<sup>+</sup><sub>rms</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1),sqrt(0.5*((*(vtin+k))+(*(vtin+zd-1-k))))/utauavg);
	for (k=up-2;k>-1;k--)
                fprintf(tcp,"%f %.9f\n",2.-((*(z+k))+1),sqrt(0.5*((*(vtin+k))+(*(vtin+zd-1-k))))/utauavg);
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"wtin.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"u<sup>+</sup><sub>rms</sub>,v<sup>+</sup><sub>rms</sub>,w<sup>+</sup><sub>rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,sqrt(0.5*((*(wtin+k))+(*(wtin+zd-1-k))))/utauavg);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"wtin-h.dat");
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z/H\", \"u<sup>+</sup><sub>rms</sub>,v<sup>+</sup><sub>rms</sub>,w<sup>+</sup><sub>rms</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1),sqrt(0.5*((*(wtin+k))+(*(wtin+zd-1-k))))/utauavg);
	for (k=up-2;k>-1;k--)
                fprintf(tcp,"%f %.9f\n",2.-((*(z+k))+1),sqrt(0.5*((*(wtin+k))+(*(wtin+zd-1-k))))/utauavg);
        fprintf(tcp,"\n");
        fclose(tcp);
        
        sprintf(fn,"KE-sum.dat");
        tcp=fopen(fn,"w");
//        fprintf(tcp,"TITLE= \"DNS Results\"\n");
//        fprintf(tcp,"VARIABLES = \"Z/H\", \"K<sup>+</sup>\"\n");
//        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
        for (k=0;k<zd;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1),((*(utin+k))+(*(vtin+k))+(*(wtin+k)))/(utauavg*utauavg));
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"rests.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"<uw><sup>+</sup>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),0.5*(-(*(reys+k))+(*(reys+zd-1-k)))/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"rests-pl.dat");
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<uw><sup>+</sup>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(-(*(reys+k))+(*(reys+zd-1-k)))/(utauavg*utauavg));
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"rests-h.dat");
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z/H\", \"<uw><sup>+</sup>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1),-0.5*(-(*(reys+k))+(*(reys+zd-1-k)))/(utauavg*utauavg));
	for (k=up-2;k>-1;k--)
                fprintf(tcp,"%f %.9f\n",2.-((*(z+k))+1),0.5*(-(*(reys+k))+(*(reys+zd-1-k)))/(utauavg*utauavg));
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"dtin.dat");
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"p<sup>'+</sup><sub>rms</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,sqrt(0.5*((*(dtin+k))+(*(dtin+zd-1-k))))/(utauavg*utauavg));
        fprintf(tcp,"\n");
        fclose(tcp);

        sprintf(fn,"dtin-long.dat");
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"p<sup>'+</sup><sub>rms</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
        for (k=0;k<zd;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1)-1.,sqrt((*(dtin+k)))/(utauavg*utauavg));
        fprintf(tcp,"\n");
        fclose(tcp);	

	sprintf(fn,"dudx.dat");
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<greek>n</greek>dU/dx/u<sub><greek>t</greek></sub><sup>2</sup>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,nu*(*(dudx+k))/(utauavg*utauavg));
        fprintf(tcp,"\n");
        fclose(tcp);

        sprintf(fn,"dudx+tw-log.dat");
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<u'w'><sup>+</sup>+<greek>n</greek>dU/dx/u<sub><greek>t</greek></sub><sup>2</sup>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(-(*(reys+k))+(*(reys+zd-1-k)))/(utauavg*utauavg)+nu*(*(dudx+k))/(utauavg*utauavg));
        fprintf(tcp,"\n");
        fclose(tcp);

        sprintf(fn,"dudx+tw-nolog.dat");
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<u'w'><sup>+</sup>+<greek>n</greek>dU/dx/u<sub><greek>t</greek></sub><sup>2</sup>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1),0.5*(-(*(reys+k))+(*(reys+zd-1-k)))/(utauavg*utauavg)+nu*(*(dudx+k))/(utauavg*utauavg));
        fprintf(tcp,"\n");
        fclose(tcp);
        
        sprintf(fn,"skew-u.dat");
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"S<sub>u</sub>, S<sub>v</sub>, S<sub>w</sub>, S<sub>p</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(sku+k))+(*(sku+zd-1-k))));
        fprintf(tcp,"\n");
        fclose(tcp);
        
        sprintf(fn,"skew-u-long.dat");
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"S<sub>u</sub>, S<sub>v</sub>, S<sub>w</sub>, S<sub>p</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
        for (k=0;k<zd;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(sku+k)));
        fprintf(tcp,"\n");
        fclose(tcp);
        
        sprintf(fn,"skew-v.dat");
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"S<sub>u</sub>, S<sub>v</sub>, S<sub>w</sub>, S<sub>p</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(skv+k))+(*(skv+zd-1-k))));
        fprintf(tcp,"\n");
        fclose(tcp);
        
        sprintf(fn,"skew-v-long.dat");
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"S<sub>u</sub>, S<sub>v</sub>, S<sub>w</sub>, S<sub>p</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
        for (k=0;k<zd;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(skv+k)));
        fprintf(tcp,"\n");
        fclose(tcp);
        
        sprintf(fn,"skew-w.dat");
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"S<sub>u</sub>, S<sub>v</sub>, S<sub>w</sub>, S<sub>p</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(skw+k))-(*(skw+zd-1-k))));
        fprintf(tcp,"\n");
        fclose(tcp);
        
        sprintf(fn,"skew-w-long.dat");
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"S<sub>u</sub>, S<sub>v</sub>, S<sub>w</sub>, S<sub>p</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
        for (k=0;k<zd;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(skw+k)));
        fprintf(tcp,"\n");
        fclose(tcp);
        
        sprintf(fn,"skew-p.dat");
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"S<sub>u</sub>, S<sub>v</sub>, S<sub>w</sub>, S<sub>p</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(skp+k))+(*(skp+zd-1-k))));
        fprintf(tcp,"\n");
        fclose(tcp);
        
        sprintf(fn,"skew-p-long.dat");
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"S<sub>u</sub>, S<sub>v</sub>, S<sub>w</sub>, S<sub>p</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
        for (k=0;k<zd;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(skp+k)));
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"skew-all.dat");
        tcp=fopen(fn,"w");
        for (k=0;k<up;k++)
                fprintf(tcp,"%d %f %.9f %.9f %.9f %.9f\n",k,((*(z+k))+1)*Ret,0.5*((*(sku+k))+(*(sku+zd-1-k))),0.5*((*(skw+k))-(*(skw+zd-1-k))),0.5*((*(skv+k))+(*(skv+zd-1-k))),0.5*((*(skp+k))+(*(skp+zd-1-k))));
        fprintf(tcp,"\n");
        fclose(tcp);        

        sprintf(fn,"skew-uw-long.dat");
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"S<sub>uw</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
        for (k=0;k<zd;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(skuw+k)));
        fprintf(tcp,"\n");
        fclose(tcp);
        
        
        sprintf(fn,"kurt-u.dat");
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"F<sub>u</sub>, F<sub>v</sub>, F<sub>w</sub>, F<sub>p</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(kuu+k))+(*(kuu+zd-1-k))));
        fprintf(tcp,"\n");
        fclose(tcp);
        
        sprintf(fn,"kurt-u-long.dat");
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"F<sub>u</sub>, F<sub>v</sub>, F<sub>w</sub>, F<sub>p</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
        for (k=0;k<zd;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(kuu+k)));
        fprintf(tcp,"\n");
        fclose(tcp);
        
        sprintf(fn,"kurt-v.dat");
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"F<sub>u</sub>, F<sub>v</sub>, F<sub>w</sub>, F<sub>p</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(kuv+k))+(*(kuv+zd-1-k))));
        fprintf(tcp,"\n");
        fclose(tcp);
        
        sprintf(fn,"kurt-v-long.dat");
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"F<sub>u</sub>, F<sub>v</sub>, F<sub>w</sub>, F<sub>p</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
        for (k=0;k<zd;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(kuv+k)));
        fprintf(tcp,"\n");
        fclose(tcp);
        
        sprintf(fn,"kurt-w.dat");
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"F<sub>u</sub>, F<sub>v</sub>, F<sub>w</sub>, F<sub>p</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(kuw+k))+(*(kuw+zd-1-k))));
        fprintf(tcp,"\n");
        fclose(tcp);
        
        sprintf(fn,"kurt-w-long.dat");
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"F<sub>u</sub>, F<sub>v</sub>, F<sub>w</sub>, F<sub>p</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
        for (k=0;k<zd;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(kuw+k)));
        fprintf(tcp,"\n");
        fclose(tcp);
        
        sprintf(fn,"kurt-p.dat");
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"F<sub>u</sub>, F<sub>v</sub>, F<sub>w</sub>, F<sub>p</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(kup+k))+(*(kup+zd-1-k))));
        fprintf(tcp,"\n");
        fclose(tcp);
        
        sprintf(fn,"kurt-p-long.dat");
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"F<sub>u</sub>, F<sub>v</sub>, F<sub>w</sub>, F<sub>p</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
        for (k=0;k<zd;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(kup+k)));
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"kurt-all.dat");
        tcp=fopen(fn,"w");
        for (k=0;k<up;k++)
                fprintf(tcp,"%d %f %.9f %.9f %.9f %.9f\n",k,((*(z+k))+1)*Ret,0.5*((*(kuu+k))+(*(kuu+zd-1-k))),0.5*((*(kuw+k))+(*(kuw+zd-1-k))),0.5*((*(kuv+k))+(*(kuv+zd-1-k))),0.5*((*(kup+k))+(*(kup+zd-1-k))));
        fprintf(tcp,"\n");
        fclose(tcp);
        
        sprintf(fn,"kurt-uw-long.dat");
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"F<sub>uw</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
        for (k=0;k<zd;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(kuuw+k)));
        fprintf(tcp,"\n");
        fclose(tcp);

	
	free(uavg);
	free(vavg);
	free(wavg);
	free(davg);
	free(utin);
	free(vtin);
	free(wtin);
	free(reys);
	free(dtin);

	free(uavgt);
	free(vavgt);
	free(wavgt);
	free(davgt);
	free(utint);
	free(vtint);
	free(wtint);
	free(reyst);
	free(dtint);
	free(du);
	free(dudx);

	free(z);
	
	free(sku);
	free(skv);
	free(skw);
	free(skuw);
	free(skut);
	free(skvt);
	free(skwt);
	free(skuwt);
	
	free(skp);
	free(skpt);
	
	free(kuu);
	free(kuv);
	free(kuw);
	free(kuuw);
	free(kuut);
	free(kuvt);
	free(kuwt);
	free(kuuwt);
	
	free(kup);
	free(kupt);
}
