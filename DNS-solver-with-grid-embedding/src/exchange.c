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

/* Exchanging the pre- and post- collision distributions */
POINTER exchange(POINTER V)
{
	double *adr,*taddrs;
  
	taddrs=V.f;
	V.f=V.ftemp;
	V.ftemp=taddrs;

  return V;
}