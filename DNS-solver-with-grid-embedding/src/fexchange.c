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

POINTER fexchange(POINTER V)
{
	double *adr,*taddrs;
	
	taddrs=V.suf;
	V.suf=V.suftemp;
	V.suftemp=taddrs;
	
	taddrs=V.slf;
	V.slf=V.slftemp;
	V.slftemp=taddrs;

  return V;
}