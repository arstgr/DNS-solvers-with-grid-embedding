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

/* bicubic interpolation within a grid square. Input quantities are p1[0.3], .., p4[0..3], x1l, x1u the lower and upper 
 * coordinates of the grid square in 1 direction, x2l, x2u likewise for the 2 direction and x1 and x2 the coordinates of the 
 * desired point for the interpolation. The interpolated function value is returned as ans. This function calls bcucof */
//__inline double bcuint(int xdc, int ydc, int GR, int i, int j, double *addr, int a)
double bcuint(int xdc, int ydc, int GR, int i, int j, double *addr, int a)
{
  int ii;
  double t,u;
  double c[4][4]={0.};
  double ans=0.;
  int ic, jc, tmp;
  
  ic = i/GR; jc =j/GR;
  
  bcucof(xdc, ydc, ic, jc, a, GR, addr, &c[0][0]);
  
  tmp = i%GR; t = ((double)tmp)/((double)GR);
  tmp = j%GR; u = ((double)tmp)/((double)GR);
  
  for (ii=3;ii>-1;ii--)
  {
    ans=t*ans+((c[ii][3]*u+c[ii][2])*u+c[ii][1])*u+c[ii][0];
  }
//  fprintf(stderr,"in function bcuint, c error handler is %s\n",strerror(errno));
  return ans;
}
