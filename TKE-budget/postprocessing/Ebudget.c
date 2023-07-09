#include "definitions.h"

void Ebudget(int xdc, int ydc, int zdc, int xdf, int df, int zdf, int GR,int st, int en)
{
	int i,j,k,a;
	FILE *sv,*sw;
	FILE *tcp,*wstr;
	char fn[20];
	int num=(en-st)+1;
	DP *z,dplus,Ret;
	DP utauavg=0.,utvt,time,nu;
	int cnt=0,ct=0,rd;
	int up;
	DP ubavg=0.,ubt=0.,twall,val,tau;
	char stri[30];
	DP uttemp,uttemp1,ttemp;
	DP *Pii,*Puu,*Pvv,*Pww,*Puw,*Paiuu,*Paivv,*Paiww,*Paiuw;
	DP *Piit,*Puut,*Pvvt,*Pwwt,*Puwt,*Paiuut,*Paivvt,*Paiwwt,*Paiuwt;
	DP *Tiir,*Tuur,*Tvvr,*Twwr,*Tuwr;
	DP *Tiirt,*Tuurt,*Tvvrt,*Twwrt,*Tuwrt;
	DP *Tiipai,*Tuupai,*Tvvpai,*Twwpai,*Tuwpai;
	DP *Tiipait,*Tuupait,*Tvvpait,*Twwpait,*Tuwpait;
	DP *Tiis,*Tuus,*Tvvs,*Twws,*Tuws;
	DP *Tiist,*Tuust,*Tvvst,*Twwst,*Tuwst;
	DP *epsilon,*Epuu,*Epvv,*Epww,*Epuw;
	DP *epsilont,*Epuut,*Epvvt,*Epwwt,*Epuwt;
	DP *KEt,*KE;
	double *Tiis2t,*Tiis2,*Tuus2t,*Tuus2,*Tvvs2t,*Tvvs2,*Twws2t,*Twws2,*Tuws2,*Tuws2t;
	double *Epii2,*Epii2t,*Epuu2,*Epuu2t,*Epvv2,*Epvv2t,*Epww2,*Epww2t,*Epuw2,*Epuw2t;
	int NZ=2*zdf+(zdc-4);
        int NZs=2*zdf+GR*(zdc-4)+GR-1;
	int zd;
	
	zd=NZ;
	
	Pii = (DP *)calloc(zd,sizeof(DP));
	Puu = (DP *)calloc(zd,sizeof(DP));
	Pvv = (DP *)calloc(zd,sizeof(DP));
	Pww = (DP *)calloc(zd,sizeof(DP));
	Puw = (DP *)calloc(zd,sizeof(DP));
	
	Paiuu = (DP *)calloc(zd,sizeof(DP));
	Paivv = (DP *)calloc(zd,sizeof(DP));
	Paiww = (DP *)calloc(zd,sizeof(DP));
	Paiuw = (DP *)calloc(zd,sizeof(DP));
	
	Piit = (DP *)calloc(zd,sizeof(DP));
	Puut = (DP *)calloc(zd,sizeof(DP));
	Pvvt = (DP *)calloc(zd,sizeof(DP));
	Pwwt = (DP *)calloc(zd,sizeof(DP));
	Puwt = (DP *)calloc(zd,sizeof(DP));
	
	Paiuut = (DP *)calloc(zd,sizeof(DP));
	Paivvt = (DP *)calloc(zd,sizeof(DP));
	Paiwwt = (DP *)calloc(zd,sizeof(DP));
	Paiuwt = (DP *)calloc(zd,sizeof(DP));
	
	Tiir = (DP *)calloc(zd,sizeof(DP));
	Tuur = (DP *)calloc(zd,sizeof(DP));
	Tvvr = (DP *)calloc(zd,sizeof(DP));
	Twwr = (DP *)calloc(zd,sizeof(DP));
	Tuwr = (DP *)calloc(zd,sizeof(DP));
	
	Tiirt = (DP *)calloc(zd,sizeof(DP));
	Tuurt = (DP *)calloc(zd,sizeof(DP));
	Tvvrt = (DP *)calloc(zd,sizeof(DP));
	Twwrt = (DP *)calloc(zd,sizeof(DP));
	Tuwrt = (DP *)calloc(zd,sizeof(DP));
	
	Tiipai = (DP *)calloc(zd,sizeof(DP));
	Tuupai = (DP *)calloc(zd,sizeof(DP));
	Tvvpai = (DP *)calloc(zd,sizeof(DP));
	Twwpai = (DP *)calloc(zd,sizeof(DP));
	Tuwpai = (DP *)calloc(zd,sizeof(DP));
	
	Tiipait = (DP *)calloc(zd,sizeof(DP));
	Tuupait = (DP *)calloc(zd,sizeof(DP));
	Tvvpait = (DP *)calloc(zd,sizeof(DP));
	Twwpait = (DP *)calloc(zd,sizeof(DP));
	Tuwpait = (DP *)calloc(zd,sizeof(DP));
	
	Tiis = (DP *)calloc(zd,sizeof(DP));
	Tuus = (DP *)calloc(zd,sizeof(DP));
	Tvvs = (DP *)calloc(zd,sizeof(DP));
	Twws = (DP *)calloc(zd,sizeof(DP));
	Tuws = (DP *)calloc(zd,sizeof(DP));
	
	Tiist = (DP *)calloc(zd,sizeof(DP));
	Tuust = (DP *)calloc(zd,sizeof(DP));
	Tvvst = (DP *)calloc(zd,sizeof(DP));
	Twwst = (DP *)calloc(zd,sizeof(DP));
	Tuwst = (DP *)calloc(zd,sizeof(DP));
	
	epsilon = (DP *)calloc(zd,sizeof(DP));
	Epuu = (DP *)calloc(zd,sizeof(DP));
	Epvv = (DP *)calloc(zd,sizeof(DP));
	Epww = (DP *)calloc(zd,sizeof(DP));
	Epuw = (DP *)calloc(zd,sizeof(DP));
	
	epsilont = (DP *)calloc(zd,sizeof(DP));
	Epuut = (DP *)calloc(zd,sizeof(DP));
	Epvvt = (DP *)calloc(zd,sizeof(DP));
	Epwwt = (DP *)calloc(zd,sizeof(DP));
	Epuwt = (DP *)calloc(zd,sizeof(DP));
	
	Paiuu = (DP *)calloc(zd,sizeof(DP));
	Paivv = (DP *)calloc(zd,sizeof(DP));
	Paiww = (DP *)calloc(zd,sizeof(DP));
	Paiuw = (DP *)calloc(zd,sizeof(DP));
	
	Paiuut = (DP *)calloc(zd,sizeof(DP));
	Paivvt = (DP *)calloc(zd,sizeof(DP));
	Paiwwt = (DP *)calloc(zd,sizeof(DP));
	Paiuwt = (DP *)calloc(zd,sizeof(DP));
	
	KEt = (DP *)calloc(zd,sizeof(DP));
	KE = (DP *)calloc(zd,sizeof(DP));

	Epii2 = (DP *)calloc(zd,sizeof(DP));
	Epii2t = (DP *)calloc(zd,sizeof(DP));
	Epuu2 = (DP *)calloc(zd,sizeof(DP));
	Epuu2t = (DP *)calloc(zd,sizeof(DP));
	Epvv2 = (DP *)calloc(zd,sizeof(DP));
	Epvv2t = (DP *)calloc(zd,sizeof(DP));
	Epww2 = (DP *)calloc(zd,sizeof(DP));
	Epww2t = (DP *)calloc(zd,sizeof(DP));
	Epuw2 = (DP *)calloc(zd,sizeof(DP));
	Epuw2t = (DP *)calloc(zd,sizeof(DP));

	Tiis2 = (DP *)calloc(zd,sizeof(DP));
	Tiis2t = (DP *)calloc(zd,sizeof(DP));
	Tuus2 = (DP *)calloc(zd,sizeof(DP));
	Tuus2t = (DP *)calloc(zd,sizeof(DP));
	Tvvs2 = (DP *)calloc(zd,sizeof(DP));
	Tvvs2t = (DP *)calloc(zd,sizeof(DP));
	Twws2 = (DP *)calloc(zd,sizeof(DP));
	Twws2t = (DP *)calloc(zd,sizeof(DP));
	Tuws2 = (DP *)calloc(zd,sizeof(DP));
	Tuws2t = (DP *)calloc(zd,sizeof(DP));

	up=1+(zd-1)/2;
	num/=DT;
	z=(DP *)calloc(zd,sizeof(DP));


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
	printf("In Budget ut=%f nu=%f\n",utauavg,nu);

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
		sprintf(fn,"KE-Budget-2d.%.3d.%.4d",0,k);
		sv=fopen(fn,"rb");
		
		fread(Piit,sizeof(DP),zd,sv);
		fread(Puut,sizeof(DP),zd,sv);
		fread(Pvvt,sizeof(DP),zd,sv);
		fread(Pwwt,sizeof(DP),zd,sv);
		fread(Puwt,sizeof(DP),zd,sv);
		
		fread(Tiirt,sizeof(DP),zd,sv);
		fread(Tuurt,sizeof(DP),zd,sv);
		fread(Tvvrt,sizeof(DP),zd,sv);
		fread(Twwrt,sizeof(DP),zd,sv);
		fread(Tuwrt,sizeof(DP),zd,sv);
		
		fread(Tiipait,sizeof(DP),zd,sv);
		fread(Tuupait,sizeof(DP),zd,sv);
		fread(Tvvpait,sizeof(DP),zd,sv);
		fread(Twwpait,sizeof(DP),zd,sv);
		fread(Tuwpait,sizeof(DP),zd,sv);
		
		fread(Paiuut,sizeof(DP),zd,sv);
		fread(Paivvt,sizeof(DP),zd,sv);
		fread(Paiwwt,sizeof(DP),zd,sv);
		fread(Paiuwt,sizeof(DP),zd,sv);
		
		fread(Tiist,sizeof(DP),zd,sv);
		fread(Tuust,sizeof(DP),zd,sv);
		fread(Tvvst,sizeof(DP),zd,sv);
		fread(Twwst,sizeof(DP),zd,sv);
		fread(Tuwst,sizeof(DP),zd,sv);
		
		fread(epsilont,sizeof(DP),zd,sv);
		fread(Epuut,sizeof(DP),zd,sv);
		fread(Epvvt,sizeof(DP),zd,sv);
		fread(Epwwt,sizeof(DP),zd,sv);
		fread(Epuwt,sizeof(DP),zd,sv);
		
		fread(KEt,sizeof(DP),zd,sv);		
		fclose(sv);

		for (a=0;a<zd;a++)
		{
			*(Pii+a) += (*(Piit+a));
			*(Puu+a) += (*(Puut+a));
			*(Pvv+a) += (*(Pvvt+a));
			*(Pww+a) += (*(Pwwt+a));
			*(Puw+a) += (*(Puwt+a));
			
			*(Tiir+a) += (*(Tiirt+a));
			*(Tuur+a) += (*(Tuurt+a));
			*(Tvvr+a) += (*(Tvvrt+a));
			*(Twwr+a) += (*(Twwrt+a));
			*(Tuwr+a) += (*(Tuwrt+a));
			
			*(Tiipai+a) += (*(Tiipait+a));
			*(Tuupai+a) += (*(Tuupait+a));
			*(Tvvpai+a) += (*(Tvvpait+a));
			*(Twwpai+a) += (*(Twwpait+a));
			*(Tuwpai+a) += (*(Tuwpait+a));
			
			*(Tiis+a) += (*(Tiist+a));
			*(Tuus+a) += (*(Tuust+a));
			*(Tvvs+a) += (*(Tvvst+a));
			*(Twws+a) += (*(Twwst+a));
			*(Tuws+a) += (*(Tuwst+a));
			
			*(epsilon+a) += (*(epsilont+a));
			*(Epuu+a) += (*(Epuut+a));
			*(Epvv+a) += (*(Epvvt+a));
			*(Epww+a) += (*(Epwwt+a));
			*(Epuw+a) += (*(Epuwt+a));

			*(Paiuu+a) += (*(Paiuut+a));
			*(Paivv+a) += (*(Paivvt+a));
			*(Paiww+a) += (*(Paiwwt+a));
			*(Paiuw+a) += (*(Paiuwt+a));
			
			*(KE+a) += (*(KEt+a));
		}
	}

	for (a=0;a<zd;a++)
	{
		*(Pii+a) /= num;
		*(Puu+a) /= num;
		*(Pvv+a) /= num;
		*(Pww+a) /= num;
		*(Puw+a) /= num;
			
		*(Tiir+a) /= num;
		*(Tuur+a) /= num;
		*(Tvvr+a) /= num;
		*(Twwr+a) /= num;
		*(Tuwr+a) /= num;
			
		*(Tiipai+a) /= num;
		*(Tuupai+a) /= num;
		*(Tvvpai+a) /= num;
		*(Twwpai+a) /= num;
		*(Tuwpai+a) /= num;
			
		*(Tiis+a) /= num;
		*(Tuus+a) /= num;
		*(Tvvs+a) /= num;
		*(Twws+a) /= num;
		*(Tuws+a) /= num;
			
		*(epsilon+a) /= num;
		*(Epuu+a) /= num;
		*(Epvv+a) /= num;
		*(Epww+a) /= num;
		*(Epuw+a) /= num;

		*(Paiuu+a) /= num;
		*(Paivv+a) /= num;
		*(Paiww+a) /= num;
		*(Paiuw+a) /= num;
		
		*(KE+a) /= num;
	}
	
	sprintf(stri,TEXT);
//	Ret=1.;

//	strcat(stri," P<sub>ii</sub>");
	sprintf(fn,"Pii.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Pii+k))+(*(Pii+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," P<sub>uu</sub>");
	sprintf(fn,"Puu.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Puu+k))+(*(Puu+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," P<sub>vv</sub>");
	sprintf(fn,"Pvv.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Pvv+k))+(*(Pvv+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," P<sub>ww</sub>");
	sprintf(fn,"Pww.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Pww+k))+(*(Pww+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," P<sub>uw</sub>");
	sprintf(fn,"Puw.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Puw+k))-(*(Puw+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," T<sub>ii</sub><sup>(r)</sup>");
	sprintf(fn,"Tiir.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tiir+k))+(*(Tiir+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," T<sub>uu</sub><sup>(r)</sup>");
	sprintf(fn,"Tuur.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tuur+k))+(*(Tuur+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," T<sub>vv</sub><sup>(r)</sup>");
	sprintf(fn,"Tvvr.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tvvr+k))+(*(Tvvr+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," T<sub>ww</sub><sup>(r)</sup>");
	sprintf(fn,"Twwr.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Twwr+k))+(*(Twwr+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," T<sub>uw</sub><sup>(r)</sup>");
	sprintf(fn,"Tuwr.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tuwr+k))-(*(Tuwr+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," T<sub>ii</sub><sup>(<greek>p</greek>)</sup>");
	sprintf(fn,"Tiipai.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tiipai+k))+(*(Tiipai+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," T<sub>uu</sub><sup>(<greek>p</greek>)</sup>");
	sprintf(fn,"Tuupai.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tuupai+k))+(*(Tuupai+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," T<sub>vv</sub><sup>(<greek>p</greek>)</sup>");
	sprintf(fn,"Tvvpai.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tvvpai+k))+(*(Tvvpai+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," T<sub>ww</sub><sup>(<greek>p</greek>)</sup>");
	sprintf(fn,"Twwpai.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Twwpai+k))+(*(Twwpai+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," T<sub>uw</sub><sup>(<greek>p</greek>)</sup>");
	sprintf(fn,"Tuwpai.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tuwpai+k))-(*(Tuwpai+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," T<sub>ii</sub><sup>(s)</sup>");
	sprintf(fn,"Tiis.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tiis+k))+(*(Tiis+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," T<sub>uu</sub><sup>(s)</sup>");
	sprintf(fn,"Tuus.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tuus+k))+(*(Tuus+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," T<sub>vv</sub><sup>(s)</sup>");
	sprintf(fn,"Tvvs.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tvvs+k))+(*(Tvvs+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," T<sub>ww</sub><sup>(s)</sup>");
	sprintf(fn,"Twws.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Twws+k))+(*(Twws+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," T<sub>uw</sub><sup>(s)</sup>");
	sprintf(fn,"Tuws.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tuws+k))-(*(Tuws+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," <greek>e</greek><sub>ii</sub>");
	sprintf(fn,"Epii.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,-0.5*((*(epsilon+k))+(*(epsilon+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," <greek>e</greek><sub>uu</sub>");
	sprintf(fn,"Epuu.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,-0.5*((*(Epuu+k))+(*(Epuu+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," <greek>e</greek><sub>vv</sub>");
	sprintf(fn,"Epvv.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,-0.5*((*(Epvv+k))+(*(Epvv+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," <greek>e</greek><sub>ww</sub>");
	sprintf(fn,"Epww.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,-0.5*((*(Epww+k))+(*(Epww+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," <greek>e</greek><sub>uw</sub>");
	sprintf(fn,"Epuw.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,-0.5*((*(Epuw+k))-(*(Epuw+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," <greek>P</greek><sub>uu</sub>");
	sprintf(fn,"Paiuu.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Paiuu+k))+(*(Paiuu+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," <greek>P</greek><sub>vv</sub>");
	sprintf(fn,"Paivv.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Paivv+k))+(*(Paivv+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," <greek>P</greek><sub>ww</sub>");
	sprintf(fn,"Paiww.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Paiww+k))+(*(Paiww+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," <greek>P</greek><sub>uw</sub>");
	sprintf(fn,"Paiuw.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Paiuw+k))-(*(Paiuw+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"KE.dat");
	tcp=fopen(fn,"w");
//	fprintf(tcp,"TITLE= \"DNS Results\"\n");
//	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"K<sup>+</sup>\"\n");
//	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(KE+k))/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"tuu.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tuus+k))+(*(Tuus+zd-1-k))+(*(Tuur+k))+(*(Tuur+zd-1-k))+(*(Tuupai+k))+(*(Tuupai+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"tvv.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tvvs+k))+(*(Tvvs+zd-1-k))+(*(Tvvr+k))+(*(Tvvr+zd-1-k))+(*(Tvvpai+k))+(*(Tvvpai+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"tww.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Twws+k))+(*(Twws+zd-1-k))+(*(Twwr+k))+(*(Twwr+zd-1-k))+(*(Twwpai+k))+(*(Twwpai+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"tii.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tiis+k))+(*(Tiis+zd-1-k))+(*(Tiir+k))+(*(Tiir+zd-1-k))+(*(Tiipai+k))+(*(Tiipai+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
	free(Pii);
	free(Puu);
	free(Pvv);
	free(Pww);
	free(Puw);
	
	free(Paiuu);
	free(Paivv);
	free(Paiww);
	free(Paiuw);
	
	free(Piit);
	free(Puut);
	free(Pvvt);
	free(Pwwt);
	free(Puwt);
	
	free(Paiuut);
	free(Paivvt);
	free(Paiwwt);
	free(Paiuwt);
	
	free(Tiir);
	free(Tuur);
	free(Tvvr);
	free(Twwr);
	free(Tuwr);
	
	free(Tiirt);
	free(Tuurt);
	free(Tvvrt);
	free(Twwrt);
	free(Tuwrt);
	
	free(Tiipai);
	free(Tuupai);
	free(Tvvpai);
	free(Twwpai);
	free(Tuwpai);
	
	free(Tiipait);
	free(Tuupait);
	free(Tvvpait);
	free(Twwpait);
	free(Tuwpait);
	
	free(Tiis);
	free(Tuus);
	free(Tvvs);
	free(Twws);
	free(Tuws);	
	
	free(Tiist);
	free(Tuust);
	free(Tvvst);
	free(Twwst);
	free(Tuwst);
	
	free(epsilon);
	free(Epuu);
	free(Epvv);
	free(Epww);
	free(Epuw);
	
	free(epsilont);
	free(Epuut);
	free(Epvvt);
	free(Epwwt);
	free(Epuwt);

	free(Epii2);
	free(Epii2t);
	free(Epuu2);
	free(Epuu2t);
	free(Epvv2);
	free(Epvv2t);
	free(Epww2);
	free(Epww2t);

	free(Tiis2);
	free(Tiis2t);
	free(Tuus2);
	free(Tuus2t);
	free(Tvvs2);
	free(Tvvs2t);
	free(Twws2);
	free(Twws2t);
	
	free(z);
}