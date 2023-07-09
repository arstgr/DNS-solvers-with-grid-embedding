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
#include "fun_defs.h"

int main(int argc, char *argv[])
{
	
	int xdc,ydc,zdc,TM;
	int xdf,ydf,zdf;
	int GR;
	int i,j,k,q,a,b,c,nut;
	int cntr;
	long ss,up;
	POINTER V;
	PDATA pd;
	DP Fx=0.,Fy=0.,Fz=0.;
	DP Gx,Gy,Gz;
	FILE *mysave,*tm;
	FILE *sv;
	char fn[20],fn2[20];
	char *ch;
	MPI_Status status;
	int tst,ttt;
	DP U,Q,ux,dens,mdt;
	DP wup,wdo,dummy,davg,mu;
	DP pgrad,dl,dr;
	DP tauc,tauf;
	DP Reb, dplus,Lplus,ubulk;
	int slp,extx,exty;
	int rank,nproc;

	const double ut=0.0400257;
	const double Ret=216.13878;
	Reb=Rebs;
	ubulk=ub;
	
	GR=gratio;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	pd.myrank=rank;
	pd.numproc=nproc;
	
	if (pd.myrank != 0)
		pd.left = pd.myrank - 1;
	else
		pd.left = pd.numproc - 1;
	
	if (pd.myrank != (pd.numproc-1))
		pd.right = pd.myrank + 1;
	else
		pd.right = 0;

	zdc=zdv;
	ydc=ydv;
	xdc=1+(xdv/pd.numproc);
	
	xdf = GR*(xdc-1)+1;
	ydf = GR*ydc;
	zdf = zds+1;	
	
	V.velc=NULL;
	V.velfl=NULL;
	V.velfu=NULL;

	V.velc=(DP *)calloc((xdc+3)*ydc*zdc*3,sizeof(DP));
	V.denc=(DP *)calloc((xdc+3)*ydc*zdc,sizeof(DP));
	if (V.velc==NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.velc: %s\n",strerror(errno));
	}
	
	V.velfu=(DP *)calloc((xdf+3)*ydf*zdf*3,sizeof(DP));
	V.denfu=(DP *)calloc((xdf+3)*ydf*zdf,sizeof(DP));
	if (V.velfu==NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.velfu: %s\n",strerror(errno));
	}
	
	V.velfl=(DP *)calloc((xdf+3)*ydf*zdf*3,sizeof(DP));
	V.denfl=(DP *)calloc((xdf+3)*ydf*zdf,sizeof(DP));
	if (V.velfl==NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.velfl: %s\n",strerror(errno));
	}

	for (tst=TSTR;tst<=TEND;tst+=DT)
	{
		V=reading(xdc, ydc, zdc, xdf,ydf,zdf, V, pd, tst, ubulk);
//		V=reading2(xdc, ydc, zdc, xdf,ydf,zdf, V, pd, tst);
		V=budget2D(xdc,ydc,zdc, xdf,ydf,zdf,GR,V,pd,tst);
		if (!(tst%200))
		{
			uv_fluc(xdf, ydf, zdf, V, pd, tst, 4);
		}
	}


	free(V.velc);
	free(V.denc);
	
	free(V.velfu);
	free(V.denfu);
	
	free(V.velfl);
	free(V.denfl);

	MPI_Finalize();

	return(0);
}


