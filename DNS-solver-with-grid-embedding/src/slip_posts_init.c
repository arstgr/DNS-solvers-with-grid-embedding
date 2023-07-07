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

POINTER slip_posts_init(int xdf, int ydf, int zdf, POINTER V, PDATA pd)
{
        int i,j,k,a,b,extx,exty,extn,t1,t2;
        int gx,gy, wx,wy;
        int slip=-2;
        MPI_Status status;
        int solid=-1,periodic=1,fluid=0;
	int *sxgrid,*sygrid;
	int extxn,extyn;
	int coord[2],jf;
	int nextxpcs,nextypcs;
	
	MPI_Cart_coords(cart_grid, pd.myrank, 2, coord);
	
	extxn = XDIM*(xdf-1);
	extyn = YDIM*(ydf-1);
	
	sxgrid = (int *)calloc(extxn,sizeof(int));
	sygrid = (int *)calloc(extyn,sizeof(int));
	
	nextxpcs = coord[0]+1;
	if (nextxpcs==XDIM)
	  nextxpcs=0;
	
	nextypcs = coord[1]+1;
	if (nextypcs==YDIM)
	  nextypcs = 0;

        gx=gxplussize;
        wx=wxplussize;

        gy=gyplussize;
        wy=wyplussize;
	
	for (i=0;i<extxn;i++)
	  sxgrid[i] = solid;
	for (j=0;j<extyn;j++)
	  sygrid[j] = solid;
	
	for (i=wx/2;i<extxn;i+=(gx+wx))
	  for (a=0;a<gx;a++)
	     sxgrid[i+a] = slip;
	  
	for (j=wy/2;j<extyn;j+=(gy+wy))
	  for (a=0;a<gy;a++)
	     sygrid[j+a] = slip;
	
	extx = xdf-1;
	exty = ydf-1;
	
	k=0;
	for (i=0;i<(xdf-1);i++)
	{
		for (j=0;j<ydf;j++)
		{
			jf = coord[0]*(xdf-1) + i;
			Fs(V.sls,i,j,k,xdf,ydf,zdf)=sxgrid[jf];
		}
	}
	i=xdf-1;
	jf = nextxpcs*(xdf-1);
	for (j=0;j<ydf;j++)
	  Fs(V.sls,i,j,k,xdf,ydf,zdf)=sxgrid[jf];
	
	for (i=0;i<xdf;i++)
	{
		for (j=0;j<(ydf-1);j++)
		{
			jf = coord[1]*(ydf-1) + j;
			if (sygrid[jf] == slip)
			  Fs(V.sls,i,j,k,xdf,ydf,zdf) = sygrid[jf];
		}
	}
	j=ydf-1;
	jf = nextypcs*(ydf-1);
	if (sygrid[jf]==slip)
	{
	  for (i=0;i<xdf;i++)
	    Fs(V.sls,i,j,k,xdf,ydf,zdf) = sygrid[jf];
	}
	
	k=zdf-1;
	for (i=0;i<(xdf-1);i++)
	{
		for (j=0;j<ydf;j++)
		{
			jf = coord[0]*(xdf-1) + i;
			Fs(V.sus,i,j,k,xdf,ydf,zdf)=sxgrid[jf];
		}
	}
	i=xdf-1;
	jf = nextxpcs*(xdf-1);
	for (j=0;j<ydf;j++)
	  Fs(V.sus,i,j,k,xdf,ydf,zdf)=sxgrid[jf];
	
	for (i=0;i<xdf;i++)
	{
		for (j=0;j<(ydf-1);j++)
		{
			jf = coord[1]*(ydf-1) + j;
			if (sygrid[jf] == slip)
			  Fs(V.sus,i,j,k,xdf,ydf,zdf) = sygrid[jf];
		}
	}
	j=ydf-1;
	jf = nextypcs*(ydf-1);
	if (sygrid[jf]==slip)
	{
	  for (i=0;i<xdf;i++)
	    Fs(V.sus,i,j,k,xdf,ydf,zdf) = sygrid[jf];
	}
	    
	free(sxgrid);
	free(sygrid);
//	fprintf(stderr,"in function slip_posts_init, c error handler is %s\n",strerror(errno));

        return V;
}