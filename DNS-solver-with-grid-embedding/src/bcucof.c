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

/* Given arrays of p1[0..3], .., p4[0..3] containing function, its x and y derivative and the x-y derivatie,
 * at 4 points of a rectangular grid cell (numbered counter-clockwise from the lower left) and given d1 and
 * d2 the length of the grid in 1 and 2 directions, this function returns c[0..3][0..3] that is used by
 * function bcuint for bicubic interpolation.  */ 
//__inline void bcucof(int xdc, int ydc, int ic, int jc, int a, int GR, double *addr, double *c)
void bcucof(int xdc, int ydc, int ic, int jc, int a, int GR, double *addr, double *c)
{
  const int wt_d[16][16] = {
     {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0},
     {-3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0},
     {2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0},
     {0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0},
     {0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1},
     {0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1},
     {-3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0},
     {9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2},
     {-6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2},
     {2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0},
     {-6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1},
     {4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1}
   };
   
   int i,j,k,l;
   int in,jn;
   double xx,x[16]={0.};
   double dx;
   
   in = ic+1;
   jn = jc+1;
   
   dx=1.;
   
   x[0] = Fbi(addr,ic,jc,0,a,xdc,ydc,4,23); x[4] = dx*Fbi(addr,ic,jc,1,a,xdc,ydc,4,23); x[8] = dx*Fbi(addr,ic,jc,2,a,xdc,ydc,4,23); x[12] = dx*dx*Fbi(addr,ic,jc,3,a,xdc,ydc,4,23);
   x[1] = Fbi(addr,in,jc,0,a,xdc,ydc,4,23); x[5] = dx*Fbi(addr,in,jc,1,a,xdc,ydc,4,23); x[9] = dx*Fbi(addr,in,jc,2,a,xdc,ydc,4,23); x[13] = dx*dx*Fbi(addr,in,jc,3,a,xdc,ydc,4,23);
   x[2] = Fbi(addr,in,jn,0,a,xdc,ydc,4,23); x[6] = dx*Fbi(addr,in,jn,1,a,xdc,ydc,4,23); x[10]= dx*Fbi(addr,in,jn,2,a,xdc,ydc,4,23); x[14] = dx*dx*Fbi(addr,in,jn,3,a,xdc,ydc,4,23);
   x[3] = Fbi(addr,ic,jn,0,a,xdc,ydc,4,23); x[7] = dx*Fbi(addr,ic,jn,1,a,xdc,ydc,4,23); x[11]= dx*Fbi(addr,ic,jn,2,a,xdc,ydc,4,23); x[15] = dx*dx*Fbi(addr,ic,jn,3,a,xdc,ydc,4,23);
   
   for (i = 0; i <= 15; i++) {
     xx = 0.0;
     for (k = 0; k <= 15; k++) xx += wt_d[i][k]*x[k];
     c[i] = xx;
   }
}   