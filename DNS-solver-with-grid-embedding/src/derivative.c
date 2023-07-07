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

/* Calculation of space derivatives for the bicubic spatial interpolation */
/* All derivatives are 4th order accurate */
void derivative(int xd, int yd, double *addr)
{
	int i,j,a;
	int jp,jpp,jn,jnn;
	int ip,ipp,in,inn;
	int extx = xd;
	
/*      
*addr is filled inside the xcomputations but its derivatives are calculated here
	0: values of f_i
	1: df_i/dx
	2: df_i/dy
	3: d2f_i/dxdy  */
  
	for (i=0;i<2;i++)
	{
		for (j=0;j<yd;j++)
		{
			jp=j-1;
			jpp=j-2;
			jn=j+1;
			jnn=j+2;
			for (a=0;a<23;a++)
			{
				Fbi(addr,(i-2),j,2,a,xd,yd,4,23) = (2./3.)*(Fbi(addr,(i-2),jn,0,a,xd,yd,4,23) - Fbi(addr,(i-2),jp,0,a,xd,yd,4,23)) - (1./12.)*(Fbi(addr,(i-2),jnn,0,a,xd,yd,4,23) - Fbi(addr,(i-2),jpp,0,a,xd,yd,4,23));
				Fbi(addr,(i+extx+1),j,2,a,xd,yd,4,23) = (2./3.)*(Fbi(addr,(i+extx+1),jn,0,a,xd,yd,4,23) - Fbi(addr,(i+extx+1),jp,0,a,xd,yd,4,23)) - (1./12.)*(Fbi(addr,(i+extx+1),jnn,0,a,xd,yd,4,23) - Fbi(addr,(i+extx+1),jpp,0,a,xd,yd,4,23));
			}      
		}
	}
  
	for (i=0;i<(extx+1);i++)
	{
		for (j=0;j<yd;j++)
		{
			jp=j-1;
			jpp=j-2;
			jn=j+1;
			jnn=j+2;
				
			ip=i-1;
			ipp=i-2;
			in=i+1;
			inn=i+2;

			for (a=0;a<23;a++)
			{
				Fbi(addr,i,j,1,a,xd,yd,4,23) = (2./3.)*(Fbi(addr,in,j,0,a,xd,yd,4,23) - Fbi(addr,ip,j,0,a,xd,yd,4,23)) - (1./12.)*(Fbi(addr,inn,j,0,a,xd,yd,4,23) - Fbi(addr,ipp,j,0,a,xd,yd,4,23));
				Fbi(addr,i,j,2,a,xd,yd,4,23) = (2./3.)*(Fbi(addr,i,jn,0,a,xd,yd,4,23) - Fbi(addr,i,jp,0,a,xd,yd,4,23)) - (1./12.)*(Fbi(addr,i,jnn,0,a,xd,yd,4,23) - Fbi(addr,i,jpp,0,a,xd,yd,4,23));
			}
		}
	}
     
	for (i=0;i<(extx+1);i++)
	{
		for (j=0;j<yd;j++)
		{
			jp=j-1;
			jpp=j-2;
			jn=j+1;
			jnn=j+2;

			ip=i-1;
			ipp=i-2;
			in=i+1;
			inn=i+2;
		    
			for (a=0;a<23;a++)
				Fbi(addr,i,j,3,a,xd,yd,4,23) = (2./3.)*(Fbi(addr,in,j,2,a,xd,yd,4,23) - Fbi(addr,ip,j,2,a,xd,yd,4,23)) - (1./12.)*(Fbi(addr,inn,j,2,a,xd,yd,4,23) - Fbi(addr,ipp,j,2,a,xd,yd,4,23));
		}
	}
// 	fprintf(stderr,"in function derivatie, c error handler is %s\n",strerror(errno));

}
