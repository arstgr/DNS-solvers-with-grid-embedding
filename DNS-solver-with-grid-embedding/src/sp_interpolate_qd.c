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

/* Performs spatial quadratic interpolation using lagrange polynomials 
 * ind[3] index of 3 layers of grid in x these are values of indices
 * jnd[3] and knd[3] are also indices in y and z
 * xd, yd, zd are the extent of the grid, either coarse or fine
 * -1 <= dx, dy, dz <= 1 is the position of the point with respect to grid points
 * adr address of the grid
 * ans[19] the interpolated value */ 
void sp_interpolate_qd(int *ind, int *jnd, int *knd, int xd, int yd, int zd, double dx, double dy, double dz, double *adr, double *ans)
{
	int a, i, j, k;
	double cz[3][3][19]={0.}, cy[3][19]={0.};
	
	for (i=0;i<3;i++)
	{
		for (j=0;j<3;j++)
		{
			for (a=0;a<19;a++)
			{
				cz[i][j][a] = 0.5*dz*(dz-1.)*Fb(adr,(ind[i]+2),(jnd[j]+2),knd[0],a, xd, yd,zd,19) - (dz-1.)*(dz+1.)*Fb(adr,(ind[i]+2),(jnd[j]+2),knd[1],a,xd,yd,zd,19) + 0.5*(dz+1.)*dz*Fb(adr,(ind[i]+2),(jnd[j]+2),knd[2],a,xd,yd,zd,19);
			}
		}
	}
	for (i=0;i<3;i++)
	{
		for (a=0;a<19;a++)
		{
			cy[i][a] = 0.5*dy*(dy-1.)*cz[i][0][a] - (dy-1.)*(dy+1.)*cz[i][1][a] + 0.5*(dy+1.)*dy*cz[i][2][a];
		}
	}
	for (a=0;a<19;a++)
		ans[a] = 0.5*dx*(dx-1.)*cy[0][a] - (dx-1.)*(dx+1.)*cy[1][a] + 0.5*(dx+1.)*dx*cy[2][a];
}