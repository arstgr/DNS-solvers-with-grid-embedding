#include "definitions.h"

int main()
{
	int xdc,ydc,zdc,TM;
	int xdf,ydf,zdf;
	int GR;
	int i,j,k,q,a,b,c,nut;
	int tagl=1,tagr=2;
	int cntr;
	int ss,up;
	FILE *epp;
	const DP C=1.;
	DP dx=1.;
	DP Fx=0.,Fy=0.,Fz=0.;
	DP Gx,Gy,Gz;
	DP tau=0.;
	DP rhozero=1.,rhozeroinv;
	DP t1,t2,t3,t4,t5;
	char ch;
	FILE *mysave,*tm;
	FILE *sv;
	char fn[20],fn2[20];
	DP utime,umax,dt;
	DP utavg;
	
	GR=gratio;

	zdc=zdv;
	ydc=ydv;
	xdc=xdv;
	
	xdf = GR*xdc;
	ydf = GR*ydc;
	zdf = zds+1;	

	time_stat(xdc,ydc,zdc, xdf, ydf, zdf, GR,TSTR,TEND);
//	CFriction(xdc,ydc,zdc, xdf, ydf, zdf, GR,TSTR,TEND);
	Vorticity(xdc,ydc,zdc,xdf,ydf,zdf,GR,TSTR,TEND);
	Ebudget(xdc,ydc,zdc,xdf,ydf,zdf,GR,TSTR,TEND);

	return(0);
}
