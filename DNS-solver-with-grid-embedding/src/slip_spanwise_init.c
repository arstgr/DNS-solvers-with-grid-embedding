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

POINTER slip_spanwise_init(int xdf, int ydf, int zdf, POINTER V, PDATA pd)
{
	int i,j,jf,k,a,extx,exty;
	int gap, width;
	int slip=-2,solid=-1;;
	int coord[2], *sgrid, strt=1;
	int nextpcs;

	MPI_Cart_coords(cart_grid, pd.myrank, 2, coord);
	extx = XDIM*(xdf-1);
	sgrid = (int *)calloc(extx, sizeof(int));
	
/******************************************************************************/
	gap=gxplussize;
	width=wxplussize;
	strt=width/2;
/******************************************************************************/
	nextpcs = coord[0]+1;
	if (nextpcs==XDIM)
	  nextpcs = 0;

	for (i=0;i<extx;i++)
	  sgrid[i]=solid;
	
	for (i=width/2;i<extx;i+=(gap+width))
	  for (a=0;a<gap;a++)
	     sgrid[i+a] = slip;
	
	extx=xdf-1;
	exty=ydf;
	k=0;
	for (j=0;j<exty;j++)
	{
		for (i=0;i<extx;i++)
		{
			jf = coord[0]*(xdf-1) + i;
			Fs(V.sls,i,j,k,xdf,ydf,zdf)=sgrid[jf];
		}
	}
	i=xdf-1;
	for (j=0;j<exty;j++)
	{
	  jf = nextpcs*(xdf-1);
	  Fs(V.sls,i,j,k,xdf,ydf,zdf)=sgrid[jf];
	}
	  
	k=zdf-1;
	for (j=0;j<exty;j++)
	{
		for (i=0;i<extx;i++)
		{
			jf = coord[0]*(xdf-1) + i;
			Fs(V.sus,i,j,k,xdf,ydf,zdf)=sgrid[jf];
		}
	}
	i=xdf-1;
	for (j=0;j<exty;j++)
	{
	  jf = nextpcs*(xdf-1);
	  Fs(V.sus,i,j,k,xdf,ydf,zdf)=sgrid[jf];
	}

	free(sgrid);
//	fprintf(stderr,"in function slip_ridges_init, c error handler is %s\n",strerror(errno));
	return V;
}