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

/* Cubic Hermit Interpolation in time */
POINTER tinterpolate_hrmt(int xdc, int ydc, int zdc, POINTER V, int iter, int GR)
{
	int i,j,a;
	double c1,c2,c3,c4;
	double tmp = (((double)iter)/((double)GR));
	
	for (i=-2;i<(xdc+3);i++)
	{
		for (j=-2;j<(ydc+3);j++)
		{
			for (a=0;a<23;a++)
			{
				c1 = Fbi(V.linp,i,j,0,a,xdc,ydc,4,23);
				c2 = (Fbi(V.lin,i,j,0,a,xdc,ydc,4,23)-Fbi(V.linpp,i,j,0,a,xdc,ydc,4,23))*0.5;
				
				c3 = Fbi(V.lin,i,j,0,a,xdc,ydc,4,23);
				c4 = (3.*Fbi(V.lin,i,j,0,a,xdc,ydc,4,23)-4.*Fbi(V.linp,i,j,0,a,xdc,ydc,4,23) + Fbi(V.linpp,i,j,0,a,xdc,ydc,4,23))*0.5;
				
				Fbi(V.lint,i,j,0,a,xdc,ydc,4,23) = ((2.*tmp*tmp*tmp - 3.* tmp*tmp +1.)*c1 + (tmp*tmp*tmp - 2.* tmp*tmp + tmp)*c2 + (-2.*tmp*tmp*tmp + 3.*tmp*tmp)*c3 + (tmp*tmp*tmp - tmp*tmp)*c4);
			  }
		}
	}
	for (i=-2;i<(xdc+3);i++)
	{
		for (j=-2;j<(ydc+3);j++)
		{
			for (a=0;a<23;a++)
			{
				c1 = Fbi(V.uinp,i,j,0,a,xdc,ydc,4,23);
				c2 = (Fbi(V.uin,i,j,0,a,xdc,ydc,4,23)-Fbi(V.uinpp,i,j,0,a,xdc,ydc,4,23))*0.5;
				
				c3 = Fbi(V.uin,i,j,0,a,xdc,ydc,4,23);
				c4 = (3.*Fbi(V.uin,i,j,0,a,xdc,ydc,4,23)-4.*Fbi(V.uinp,i,j,0,a,xdc,ydc,4,23) + Fbi(V.uinpp,i,j,0,a,xdc,ydc,4,23))*0.5;
				
				Fbi(V.uint,i,j,0,a,xdc,ydc,4,23) = (2.*tmp*tmp*tmp - 3.* tmp*tmp +1.)*c1 + (tmp*tmp*tmp - 2.* tmp*tmp + tmp)*c2 + (-2.*tmp*tmp*tmp + 3.*tmp*tmp)*c3 + (tmp*tmp*tmp - tmp*tmp)*c4;
			}
		}
	}
	return V;	
}
