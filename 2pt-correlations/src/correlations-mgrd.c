/********************************************************************************
 * Amirreza Rastegari							                             	*
 * arstgri@gmail.com						                             		*
 * This code is used to calculate the 2 point velocity correlations in parallel *
 * Formulation:								                                 	*
 * Mathiu & Scott, An introduction to turbulent flows, Cambridge university	    *
 * press, 2000, pp. 62								                            *
 ********************************************************************************/
#include "definitions.h"
#include "fun_defs.h"

int main(int argc, char *argv[])
{
	int i,j,k,q,a,b,c,nut,TM;
	int xdf,ydf,zdf,xdc,ydc,zdc,GR;
	int cntr;
	long ss,up;
	POINTER V;
	PDATA pd;
	FILE *mysave,*tm;
	FILE *sv;
	char fn[60],fn2[60];
	char *ch;
	MPI_Status status;
	int tst,ttt;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &(pd.numproc));
	MPI_Comm_rank(MPI_COMM_WORLD, &(pd.myrank));
	
	if (pd.myrank != 0)
		pd.left = pd.myrank - 1;
	else
		pd.left = pd.numproc - 1;
	
	if (pd.myrank != (pd.numproc-1))
		pd.right = pd.myrank + 1;
	else
		pd.right = 0;
	
	GR=gratio;
	
	zdc=zdv;
	ydc=ydv;
	xdc=(xdv/pd.numproc);
	
	xdf = GR*xdc;
	ydf = GR*ydc;
	zdf = zds+1;	

	V.velc=(DP *)calloc(xdc*ydc*zdc*3,sizeof(DP));
	V.velfl=(DP *)calloc(xdf*ydf*zdf*3,sizeof(DP));
	V.velfu=(DP *)calloc(xdf*ydf*zdf*3,sizeof(DP));
	
	V.uavgc=(DP *)calloc(zdc,sizeof(DP));
	V.vavgc=(DP *)calloc(zdc,sizeof(DP));
	V.wavgc=(DP *)calloc(zdc,sizeof(DP));
	
	V.uavgfl=(DP *)calloc(zdf,sizeof(DP));
	V.vavgfl=(DP *)calloc(zdf,sizeof(DP));
	V.wavgfl=(DP *)calloc(zdf,sizeof(DP));
	
	V.uavgfu=(DP *)calloc(zdf,sizeof(DP));
	V.vavgfu=(DP *)calloc(zdf,sizeof(DP));
	V.wavgfu=(DP *)calloc(zdf,sizeof(DP));

	V.velrecvc=(DP *)calloc(xdc*ydc*zdc*3,sizeof(DP));
	V.velrecvfl=(DP *)calloc(xdf*ydf*zdf*3,sizeof(DP));
	V.velrecvfu=(DP *)calloc(xdf*ydf*zdf*3,sizeof(DP));
	
	V.ycorlc = (DP *)calloc(ydc*zdc*3,sizeof(double));
	V.ycorlfl = (DP *)calloc(ydf*zdf*3,sizeof(double));
	V.ycorlfu = (DP *)calloc(ydf*zdf*3,sizeof(double));
	
	V.xcorlc = (DP *)calloc(xdv*zdc*3,sizeof(double));
	V.xcorlfl = (DP *)calloc(xdv*GR*zdf*3,sizeof(double));
	V.xcorlfu = (DP *)calloc(xdv*GR*zdf*3,sizeof(double));

	V=mean(xdc, ydc, zdc, xdf, ydf, zdf, V, pd);
	MPI_Barrier(MPI_COMM_WORLD);
        printf("mean calculated\n");
        MPI_Barrier(MPI_COMM_WORLD);
	
	for (tst=TSTR;tst<=TEND;tst+=DT)
	{
		V=reading(xdc, ydc, zdc, xdf, ydf, zdf, V, pd, tst);
	//	MPI_Barrier(MPI_COMM_WORLD);
	//	printf("read the instantaneous field\n");
	//	MPI_Barrier(MPI_COMM_WORLD);	
		V=mean_subtract(xdc, ydc, zdc, xdf, ydf, zdf, V, pd);
	//	MPI_Barrier(MPI_COMM_WORLD);
	//	printf("subtracted the mean\n");
	//	MPI_Barrier(MPI_COMM_WORLD);
		V=ycorrel(xdc, ydc, zdc, xdf, ydf, zdf, V, pd);
	//	MPI_Barrier(MPI_COMM_WORLD);
	//	printf("ycorrl calculated\n");
        //      MPI_Barrier(MPI_COMM_WORLD);
		V=xcorrel(xdc, ydc, zdc, xdf, ydf, zdf, V, pd);
	//	MPI_Barrier(MPI_COMM_WORLD);
        //      printf("xcorrl calculated\n");
        //      MPI_Barrier(MPI_COMM_WORLD);
		
		if (pd.myrank==0)
		{
			sprintf(fn,"xcorl-c.%.4d",tst);
			sv=fopen(fn,"wb");
			fwrite(V.xcorlc,sizeof(double),(xdv*zdc*3),sv);
			fclose(sv);
			
			sprintf(fn,"xcorl-fl.%.4d",tst);
			sv=fopen(fn,"wb");
			fwrite(V.xcorlfl,sizeof(double),(xdv*GR*zdf*3),sv);
			fclose(sv);
			
			sprintf(fn,"xcorl-fu.%.4d",tst);
			sv=fopen(fn,"wb");
			fwrite(V.xcorlfu,sizeof(double),(xdv*GR*zdf*3),sv);
			fclose(sv);

			sprintf(fn,"ycorl-c.%.4d",tst);
			sv=fopen(fn,"wb");
			fwrite(V.ycorlc,sizeof(double),(ydc*zdc*3),sv);
			fclose(sv);
			
			sprintf(fn,"ycorl-fl.%.4d",tst);
			sv=fopen(fn,"wb");
			fwrite(V.ycorlfl,sizeof(double),(ydf*zdf*3),sv);
			fclose(sv);
			
			sprintf(fn,"ycorl-fu.%.4d",tst);
			sv=fopen(fn,"wb");
			fwrite(V.ycorlfu,sizeof(double),(ydf*zdf*3),sv);
			fclose(sv);
		}
		MPI_Barrier(MPI_COMM_WORLD); 
	}


	free(V.velc);
	free(V.velfl);
	free(V.velfu);
	
	free(V.uavgc);
	free(V.vavgc);
	free(V.wavgc);
	
	free(V.uavgfl);
	free(V.vavgfl);
	free(V.wavgfl);
	
	free(V.uavgfu);
	free(V.vavgfu);
	free(V.wavgfu);
	
	free(V.velrecvc);
	free(V.velrecvfl);
	free(V.velrecvfu);
	
	free(V.xcorlc);
	free(V.xcorlfl);
	free(V.xcorlfu);
	
	free(V.ycorlc);
	free(V.ycorlfl);
	free(V.ycorlfu);

	MPI_Finalize();

	return(0);
}


