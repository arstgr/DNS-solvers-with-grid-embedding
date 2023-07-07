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

/* Linear Interpolation in time */
POINTER tinterpolate_ln(int xdc, int ydc, int zdc, POINTER V, int iter, int GR)
{
	int i,j,a; 
	const double time = ((double)(((double)iter)/((double)GR)));
	
	for (i=-2;i<(xdc+3);i++)
	{
		for (j=-2;j<(ydc+3);j++)
		{
			for (a=0;a<23;a++)
			{
				Fbi(V.lint,i,j,0,a,xdc,ydc,4,23) = (Fbi(V.linp,i,j,0,a,xdc,ydc,4,23)+time*(Fbi(V.lin,i,j,0,a,xdc,ydc,4,23)-Fbi(V.linp,i,j,0,a,xdc,ydc,4,23)));
				Fbi(V.uint,i,j,0,a,xdc,ydc,4,23) = (Fbi(V.uinp,i,j,0,a,xdc,ydc,4,23)+time*(Fbi(V.uin,i,j,0,a,xdc,ydc,4,23)-Fbi(V.uinp,i,j,0,a,xdc,ydc,4,23)));
			}
		}
	}
	return V;	
}