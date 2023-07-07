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

POINTER initializer(int xdf, int ydf, int zdf, double tauf, int xdc, int ydc, int zdc, double tauc, int xd, int yd, int zd, double tau, int GR, POINTER V, PDATA pd, int cntr)
{
	int i,j,k,a;
	int ind[3],jnd[3],knd[3];
	double *adr;
	double dx,dy,dz;
	int iref, jref, kref;
	double dxf,dxc,dxr;
	int NZr=zd,NZs=2*zdf+GR*(zdc-4)+GR-1;
	double fintpl[19];
	char fn[40];
	int extx, exty,coord[2];
	double ux,uy,uz,dens,feq,fneq,ueq[3],v,u;
	MPI_Datatype distf,block,sarray,larray;
	MPI_Datatype cblock,sblock;
	MPI_Datatype send_to_dl,send_to_dr,send_to_dw,send_to_le,send_to_ri,send_to_ul,send_to_up,send_to_ur;
	MPI_Datatype recv_fr_dl,recv_fr_dr,recv_fr_dw,recv_fr_le,recv_fr_ri,recv_fr_ul,recv_fr_up,recv_fr_ur;
	
	MPI_Datatype FN_send_to_dl,FN_send_to_dw,FN_send_to_le;
	MPI_Datatype FN_recv_fr_ri,FN_recv_fr_up,FN_recv_fr_ur;
	
	MPI_Datatype sFN_send_to_dl,sFN_send_to_dw,sFN_send_to_le;
	MPI_Datatype sFN_recv_fr_ri,sFN_recv_fr_up,sFN_recv_fr_ur;
	
	int gsizes[2],lsizes[2],psizes[2],array_start[2],memsizes[2];
	const int NPCS=pd.numproc;
	
	MPI_File fh;
	MPI_Status tstatus;
	MPI_Offset offset;
	
	dxf = 2./((double)(NZs-2.));
	dxc = GR*dxf;
	dxr = 2./((double)(zd-2.));
	
	extx = (xd-1)*NPCS/XDIM;
	exty = yd/YDIM;
	
	printf("dxf=%f dxc=%f dxr=%f\n",dxf,dxc,dxr);
	
	adr = (double *)calloc((extx+6)*(exty+6)*zd*19,sizeof(double));
	
	MPI_Type_contiguous(19,MPI_DOUBLE,&distf);
	MPI_Type_commit(&distf);
	MPI_Type_contiguous(zd,distf,&block);
	MPI_Type_commit(&block);
	
	MPI_Type_contiguous(zdc,distf,&cblock);
	MPI_Type_commit(&cblock);
	
	MPI_Type_contiguous(zdf,distf,&sblock);
	MPI_Type_commit(&sblock);
	
	MPI_Cart_coords(cart_grid, pd.myrank, 2, coord);
	
	gsizes[0] = (xd-1)*NPCS;
	gsizes[1] = yd;
/****************************************************************************************************************/
	
	lsizes[0] = extx;
	lsizes[1] = exty;
	
	array_start[0] = coord[0] * lsizes[0];
	array_start[1] = coord[1] * lsizes[1];

	MPI_Type_create_subarray(2, gsizes, lsizes, array_start, MPI_ORDER_C, block, &sarray);
	MPI_Type_commit(&sarray);
	
	lsizes[0] = extx;
	lsizes[1] = exty;
	
	memsizes[0] = extx +6;
	memsizes[1] = exty +6;
	array_start[0] = 2;
	array_start[1] = 2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &larray);
	MPI_Type_commit(&larray);
/*****************************************************************************************************************/	
	lsizes[0] = extx;
	lsizes[1] = 2;
	array_start[0] = 2;
	array_start[1] = exty +2 -2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &send_to_up);
	MPI_Type_commit(&send_to_up);
	
	lsizes[0] = extx;
	lsizes[1] = 2;
	array_start[0] = 2;
	array_start[1] = 0;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &recv_fr_dw);
	MPI_Type_commit(&recv_fr_dw);
	
	lsizes[0] = extx;
	lsizes[1] = 4;
	array_start[0] = 2;
	array_start[1] = 2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &send_to_dw);
	MPI_Type_commit(&send_to_dw);
	
	lsizes[0] = extx;
	lsizes[1] = 4;
	array_start[0] = 2;
	array_start[1] = exty+2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &recv_fr_up);
	MPI_Type_commit(&recv_fr_up);
	
/********************************************************************************************************************************/
	lsizes[0] = 2;
	lsizes[1] = exty;
	array_start[0] = extx-2 +2;
	array_start[1] = 2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &send_to_ri);
	MPI_Type_commit(&send_to_ri);
	
	lsizes[0] = 2;
	lsizes[1] = exty;
	array_start[0] = 0;
	array_start[1] = 2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &recv_fr_le);
	MPI_Type_commit(&recv_fr_le);
	
	lsizes[0] = 4;
	lsizes[1] = exty;
	array_start[0] = 2;
	array_start[1] = 2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &send_to_le);
	MPI_Type_commit(&send_to_le);
	
	lsizes[0] = 4;
	lsizes[1] = exty;
	array_start[0] = extx + 2;
	array_start[1] = 2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &recv_fr_ri);
	MPI_Type_commit(&recv_fr_ri);
	
/*********************************************************************************************************************************/
	lsizes[0] = 2;
	lsizes[1] = 2;
	array_start[0] = extx-2 +2;
	array_start[1] = exty-2 +2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &send_to_ur);
	MPI_Type_commit(&send_to_ur);
	
	lsizes[0] = 4;
	lsizes[1] = 4;
	array_start[0] = extx+2;
	array_start[1] = exty+2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &recv_fr_ur);
	MPI_Type_commit(&recv_fr_ur);
	
	lsizes[0] = 4;
	lsizes[1] = 2;
	array_start[0] = 2;
	array_start[1] = exty -2 +2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &send_to_ul);
	MPI_Type_commit(&send_to_ul);
	
	lsizes[0] = 2;
	lsizes[1] = 4;
	array_start[0] = 0;
	array_start[1] = exty+2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &recv_fr_ul);
	MPI_Type_commit(&recv_fr_ul);
	
	lsizes[0] = 2;
	lsizes[1] = 4;
	array_start[0] = extx-2 +2;
	array_start[1] = 2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &send_to_dr);
	MPI_Type_commit(&send_to_dr);
	
	lsizes[0] = 4;
	lsizes[1] = 2;
	array_start[0] = extx +2;
	array_start[1] = 0;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &recv_fr_dr);
	MPI_Type_commit(&recv_fr_dr);
	
	lsizes[0] = 4;
	lsizes[1] = 4;
	array_start[0] = 2;
	array_start[1] = 2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &send_to_dl);
	MPI_Type_commit(&send_to_dl);
	
	lsizes[0] = 2;
	lsizes[1] = 2;
	array_start[0] = 0;
	array_start[1] = 0;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &recv_fr_dl);
	MPI_Type_commit(&recv_fr_dl);
	
/*****************************************************************************************************************************************/

	memsizes[0] = xdc;
	memsizes[1] = ydc;
	
	lsizes[0] = xdc-1;
	lsizes[1] = 1;
	array_start[0] = 0;
	array_start[1] = 0;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, cblock, &FN_send_to_dw);
	MPI_Type_commit(&FN_send_to_dw);
	
	lsizes[0] = xdc-1;
	lsizes[1] = 1;
	array_start[0] = 0;
	array_start[1] = ydc-1;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, cblock, &FN_recv_fr_up);
	MPI_Type_commit(&FN_recv_fr_up);
	
	lsizes[0] = 1;
	lsizes[1] = ydc-1;
	array_start[0] = 0;
	array_start[1] = 0;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, cblock, &FN_send_to_le);
	MPI_Type_commit(&FN_send_to_le);
	
	lsizes[0] = 1;
	lsizes[1] = ydc-1;
	array_start[0] = xdc-1;
	array_start[1] = 0;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, cblock, &FN_recv_fr_ri);
	MPI_Type_commit(&FN_recv_fr_ri);
	
	lsizes[0] = 1;
	lsizes[1] = 1;
	array_start[0] = 0;
	array_start[1] = 0;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, cblock, &FN_send_to_dl);
	MPI_Type_commit(&FN_send_to_dl);
	
	lsizes[0] = 1;
	lsizes[1] = 1;
	array_start[0] = xdc-1;
	array_start[1] = ydc-1;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, cblock, &FN_recv_fr_ur);
	MPI_Type_commit(&FN_recv_fr_ur);
	
/*****************************************************************************************************************************************/

	memsizes[0] = xdf;
	memsizes[1] = ydf;
	
	lsizes[0] = xdf-1;
	lsizes[1] = 1;
	array_start[0] = 0;
	array_start[1] = 0;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, sblock, &sFN_send_to_dw);
	MPI_Type_commit(&sFN_send_to_dw);
	
	lsizes[0] = xdf-1;
	lsizes[1] = 1;
	array_start[0] = 0;
	array_start[1] = ydf-1;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, sblock, &sFN_recv_fr_up);
	MPI_Type_commit(&sFN_recv_fr_up);
	
	lsizes[0] = 1;
	lsizes[1] = ydf-1;
	array_start[0] = 0;
	array_start[1] = 0;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, sblock, &sFN_send_to_le);
	MPI_Type_commit(&sFN_send_to_le);
	
	lsizes[0] = 1;
	lsizes[1] = ydf-1;
	array_start[0] = xdf-1;
	array_start[1] = 0;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, sblock, &sFN_recv_fr_ri);
	MPI_Type_commit(&sFN_recv_fr_ri);
	
	lsizes[0] = 1;
	lsizes[1] = 1;
	array_start[0] = 0;
	array_start[1] = 0;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, sblock, &sFN_send_to_dl);
	MPI_Type_commit(&sFN_send_to_dl);
	
	lsizes[0] = 1;
	lsizes[1] = 1;
	array_start[0] = xdf-1;
	array_start[1] = ydf-1;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, sblock, &sFN_recv_fr_ur);
	MPI_Type_commit(&sFN_recv_fr_ur);
	
/*****************************************************************************************************************************************/
  
	sprintf(fn,"lbf.%.4d",cntr); 
	MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
	MPI_File_set_view(fh,0,MPI_DOUBLE,sarray,"native",MPI_INFO_NULL);
	
// read V.ftemp	
	MPI_File_read_all(fh, adr, 1, larray, &tstatus); 
// exchange the boundaries
/*****************************************************************************************************************************************/

	MPI_Sendrecv(adr,1,send_to_ri,pd.right,11,adr,1,recv_fr_le,pd.left,11,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(adr,1,send_to_le,pd.left,12,adr,1,recv_fr_ri,pd.right,12,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(adr,1,send_to_dw,pd.dw,13,adr,1,recv_fr_up,pd.up,13,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(adr,1,send_to_up,pd.up,14,adr,1,recv_fr_dw,pd.dw,14,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(adr,1,send_to_ur,pd.ur,15,adr,1,recv_fr_dl,pd.dl,15,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(adr,1,send_to_ul,pd.ul,16,adr,1,recv_fr_dr,pd.dr,16,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(adr,1,send_to_dl,pd.dl,17,adr,1,recv_fr_ur,pd.ur,17,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(adr,1,send_to_dr,pd.dr,18,adr,1,recv_fr_ul,pd.ul,18,MPI_COMM_WORLD,&tstatus);
	
/*******************************************************************************************************************************************/
	
	
// do interpolation
// starting from lower fine grid
	for (i=0;i<(xdf-1);i++)
	{
		for (j=0;j<(ydf-1);j++)
		{
			for (k=1;k<zdf;k++)
			{
				iref = (int)floor((i+0.5)*(dxf/dxr) - 0.5);
				jref = (int)floor((j+0.5)*(dxf/dxr) - 0.5);
				kref = (int)floor((k-0.5)*(dxf/dxr) + 0.5);
				
				ind[0] = iref - 1;
				ind[1] = iref;
				ind[2] = iref + 1;
				
				dx = (i+0.5)*dxf - (iref+0.5)*dxr;
				dy = (j+0.5)*dxf - (jref+0.5)*dxr;
				
				jnd[0] = jref - 1;
				jnd[1] = jref;
				jnd[2] = jref + 1;
				
				if (kref == 0)
				{
				  knd[0] = kref;
				  knd[1] = kref + 1;
				  knd[2] = kref + 2;
				  
				  dz = -dxr + (k-0.5)*dxf + 0.5*dxr - kref*dxr; 
				}
				else
				{
				  knd[0] = kref-1;
				  knd[1] = kref;
				  knd[2] = kref+1;
				  
				  dz = (k-0.5)*dxf + 0.5*dxr - kref*dxr; 
				}
				
		//		dz = (k-0.5)*dxf + 0.5*dxr - knd[1]*dxr;
				dx = dx/dxr; dy =dy/dxr; dz=dz/dxr;
				for (a=0;a<19;a++)
				  fintpl[a] = 0.;
				
				sp_interpolate_qd(ind,jnd,knd,(extx+6),(exty+6),zd,dx,dy,dz,adr,fintpl);
				
				v = 0.;
				ux=0.;uy=0.;uz=0.;dens=0.;
				for (a=0;a<3;a++)
				  ueq[a] = 0.;
				
				for (a=0;a<19;a++)
				{
					ux += E0[a]*fintpl[a];
					uy += E1[a]*fintpl[a];
					uz += E2[a]*fintpl[a];
					dens += fintpl[a];
				}

				ueq[0]=ux/dens;
				ueq[1]=uy/dens;
				ueq[2]=uz/dens;
				
				v = ueq[0]*ueq[0]+ueq[1]*ueq[1]+ueq[2]*ueq[2];
				for (a=0;a<19;a++)
				{
					u=E0[a]*ueq[0]+E1[a]*ueq[1]+E2[a]*ueq[2];
					feq=w[a]*dens*(1.+3.*u+4.5*u*u-1.5*v);
					fneq=fintpl[a]-feq;
					Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19) = feq + (dxr/dxf)*(tau/tauf)*fneq ;
				}
				
//				if (k==2||k==3)
//				  printf("k=%d dx=%f dy=%f dz=%f i0=%d i1=%d i2=%d j0=%d j1=%d j2=%d k0=%d k1=%d k2=%d i=%d j=%d u=%f\n",k,dx,dy,dz,ind[0],ind[1],ind[2],jnd[0],jnd[1],jnd[2],knd[0],knd[1],knd[2],i,j,ueq[0]);
			}
		}
	}
// // Interpolation on to coarse grid
	for (i=0;i<(xdc-1);i++)
	{
		for (j=0;j<(ydc-1);j++)
		{
			for (k=0;k<zdc;k++)
			{
				iref = (int)floor((i+0.5)*dxc/dxr - 0.5);
				jref = (int)floor((j+0.5)*dxc/dxr - 0.5);
				kref = (int)floor(0.5 + (k*GR+(zdf-GR-1-0.5))*dxf/dxr);
				
				ind[0] = iref - 1;
				ind[1] = iref;
				ind[2] = iref + 1;
				
				dx = (i+0.5)*dxc - (iref+0.5)*dxr;
				dy = (j+0.5)*dxc - (jref+0.5)*dxr;
				
				jnd[0] = jref - 1;
				jnd[1] = jref;
				jnd[2] = jref + 1;
				
				knd[0] = kref-1;
				knd[1] = kref;
				knd[2] = kref+1;
				  
				dz = (k*GR+(zdf-GR-1-0.5))*dxf - (kref-0.5)*dxr; 
				
		//		dz = (k*dxc+(zdf-1-0.5-GR)*dxf) + 0.5*dxr - knd[1]*dxr; 
				
				for (a=0;a<19;a++)
				  fintpl[a] = 0.;
				
				dx = dx/dxr; dy =dy/dxr; dz=dz/dxr;
				sp_interpolate_qd(ind,jnd,knd,(extx+6),(exty+6),zd,dx,dy,dz,adr,fintpl);
				
		//		for (a=0;a<19;a++)
		//		  Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19) = fintpl[a];
				
				v = 0.;
				ux=0.;uy=0.;uz=0.;dens=0.;
				for (a=0;a<3;a++)
				  ueq[a] = 0.;
				
				for (a=0;a<19;a++)
				{
					ux += E0[a]*fintpl[a];
					uy += E1[a]*fintpl[a];
					uz += E2[a]*fintpl[a];
					dens += fintpl[a];
				}

				ueq[0]=ux/dens;
				ueq[1]=uy/dens;
				ueq[2]=uz/dens;
				
				v = ueq[0]*ueq[0]+ueq[1]*ueq[1]+ueq[2]*ueq[2];
				for (a=0;a<19;a++)
				{
					u=E0[a]*ueq[0]+E1[a]*ueq[1]+E2[a]*ueq[2];
					feq=w[a]*dens*(1.+3.*u+4.5*u*u-1.5*v);
					fneq=fintpl[a]-feq;
					Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19) = feq + (dxr/dxc)*(tau/tauc)*fneq ;
				}
			}
		}
	}
// Interpolation on to the upper fine grid
	for (i=0;i<(xdf-1);i++)
	{
		for (j=0;j<(ydf-1);j++)
		{
			for (k=0;k<(zdf-1);k++)
			{
				iref = (int)floor((i+0.5)*dxf/dxr - 0.5);
				jref = (int)floor((j+0.5)*dxf/dxr - 0.5);
				kref = (int)floor(0.5+ (k+GR*(zdc-3)+zdf-1-0.5)*dxf/dxr);
				
				ind[0] = iref - 1;
				ind[1] = iref;
				ind[2] = iref + 1;
				
				dx = (i+0.5)*dxf - (iref+0.5)*dxr;
				dy = (j+0.5)*dxf - (jref+0.5)*dxr;
				
				jnd[0] = jref - 1;
				jnd[1] = jref;
				jnd[2] = jref + 1;
				
				
				knd[0] = kref-1;
				knd[1] = kref;
				knd[2] = kref+1;
				
				dz = (k+GR*(zdc-3)+zdf-1-0.5)*dxf - (kref-0.5)*dxr; 
				
		//		dz = (k+GR*(zdc-3)+zdf-1-0.5)*dxf + 0.5*dxr - knd[1]*dxr; 
				
				for (a=0;a<19;a++)
				  fintpl[a] = 0.;
				
				dx = dx/dxr; dy =dy/dxr; dz=dz/dxr;
				sp_interpolate_qd(ind,jnd,knd,(extx+6),(exty+6),zd,dx,dy,dz,adr,fintpl);

				v = 0.;
				ux=0.;uy=0.;uz=0.;dens=0.;
				for (a=0;a<3;a++)
				  ueq[a] = 0.;
				
				for (a=0;a<19;a++)
				{
					ux += E0[a]*fintpl[a];
					uy += E1[a]*fintpl[a];
					uz += E2[a]*fintpl[a];
					dens += fintpl[a];
				}

				ueq[0]=ux/dens;
				ueq[1]=uy/dens;
				ueq[2]=uz/dens;
				
				v = ueq[0]*ueq[0]+ueq[1]*ueq[1]+ueq[2]*ueq[2];
				for (a=0;a<19;a++)
				{
					u=E0[a]*ueq[0]+E1[a]*ueq[1]+E2[a]*ueq[2];
					feq=w[a]*dens*(1.+3.*u+4.5*u*u-1.5*v);
					fneq=fintpl[a]-feq;
					Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19) = feq + (dxr/dxf)*(tau/tauf)*fneq ;
				}				
			}
		}
	}
	MPI_Sendrecv(V.suftemp,1,sFN_send_to_dl,pd.dl,28,V.suftemp,1,sFN_recv_fr_ur,pd.ur,28,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(V.suftemp,1,sFN_send_to_dw,pd.dw,29,V.suftemp,1,sFN_recv_fr_up,pd.up,29,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(V.suftemp,1,sFN_send_to_le,pd.left,30,V.suftemp,1,sFN_recv_fr_ri,pd.right,30,MPI_COMM_WORLD,&tstatus);
	
	MPI_Sendrecv(V.ftemp,1,FN_send_to_dl,pd.dl,31,V.ftemp,1,FN_recv_fr_ur,pd.ur,31,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(V.ftemp,1,FN_send_to_dw,pd.dw,32,V.ftemp,1,FN_recv_fr_up,pd.up,32,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(V.ftemp,1,FN_send_to_le,pd.left,33,V.ftemp,1,FN_recv_fr_ri,pd.right,33,MPI_COMM_WORLD,&tstatus);
	
	MPI_Sendrecv(V.slftemp,1,sFN_send_to_dl,pd.dl,34,V.slftemp,1,sFN_recv_fr_ur,pd.ur,34,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(V.slftemp,1,sFN_send_to_dw,pd.dw,35,V.slftemp,1,sFN_recv_fr_up,pd.up,35,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(V.slftemp,1,sFN_send_to_le,pd.left,36,V.slftemp,1,sFN_recv_fr_ri,pd.right,36,MPI_COMM_WORLD,&tstatus);
	
	
	MPI_File_set_view(fh,(xd-1)*yd*zd*19*pd.numproc*sizeof(double),MPI_DOUBLE,sarray,"native",MPI_INFO_NULL);
// read V.f	
	MPI_File_read_all(fh, adr, 1, larray, &tstatus); 
// exchange the boundaries
/*****************************************************************************************************************************************/

	MPI_Sendrecv(adr,1,send_to_ri,pd.right,11,adr,1,recv_fr_le,pd.left,11,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(adr,1,send_to_le,pd.left,12,adr,1,recv_fr_ri,pd.right,12,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(adr,1,send_to_dw,pd.dw,13,adr,1,recv_fr_up,pd.up,13,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(adr,1,send_to_up,pd.up,14,adr,1,recv_fr_dw,pd.dw,14,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(adr,1,send_to_ur,pd.ur,15,adr,1,recv_fr_dl,pd.dl,15,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(adr,1,send_to_ul,pd.ul,16,adr,1,recv_fr_dr,pd.dr,16,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(adr,1,send_to_dl,pd.dl,17,adr,1,recv_fr_ur,pd.ur,17,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(adr,1,send_to_dr,pd.dr,18,adr,1,recv_fr_ul,pd.ul,18,MPI_COMM_WORLD,&tstatus);
	
/*******************************************************************************************************************************************/
// Close the files
	MPI_File_close(&fh);
// do interpolation
// On lower fine grid
	for (i=0;i<(xdf-1);i++)
	{
		for (j=0;j<(ydf-1);j++)
		{
			for (k=1;k<zdf;k++)
			{
				iref = (int)floor((i+0.5)*dxf/dxr - 0.5);
				jref = (int)floor((j+0.5)*dxf/dxr - 0.5);
				kref = (int)floor((k-0.5)*(dxf/dxr) + 0.5);
				
				ind[0] = iref - 1;
				ind[1] = iref;
				ind[2] = iref + 1;
				
				dx = (i+0.5)*dxf - (iref+0.5)*dxr;
				dy = (j+0.5)*dxf - (jref+0.5)*dxr;
				
				jnd[0] = jref - 1;
				jnd[1] = jref;
				jnd[2] = jref + 1;
				
				if (kref == 0)
				{
				  knd[0] = kref;
				  knd[1] = kref + 1;
				  knd[2] = kref + 2;
				  
				  dz = -dxr + (k-0.5)*dxf + 0.5*dxr - kref*dxr; 
				}
				else
				{
				  knd[0] = kref-1;
				  knd[1] = kref;
				  knd[2] = kref+1;
				  
				  dz = (k-0.5)*dxf + 0.5*dxr - kref*dxr; 
				}
		//		dz = (k-0.5)*dxf + 0.5*dxr - knd[1]*dxr; 
				
				dx = dx/dxr; dy =dy/dxr; dz=dz/dxr;
				for (a=0;a<19;a++)
				  fintpl[a] = 0.;
				
				sp_interpolate_qd(ind,jnd,knd,(extx+6),(exty+6),zd,dx,dy,dz,adr,fintpl);
				
				v = 0.;
				ux=0.;uy=0.;uz=0.;dens=0.;
				for (a=0;a<3;a++)
				  ueq[a] = 0.;
				
				for (a=0;a<19;a++)
				{
					ux += E0[a]*fintpl[a];
					uy += E1[a]*fintpl[a];
					uz += E2[a]*fintpl[a];
					dens += fintpl[a];
				}

				ueq[0]=ux/dens;
				ueq[1]=uy/dens;
				ueq[2]=uz/dens;
				
				v = ueq[0]*ueq[0]+ueq[1]*ueq[1]+ueq[2]*ueq[2];
				for (a=0;a<19;a++)
				{
					u=E0[a]*ueq[0]+E1[a]*ueq[1]+E2[a]*ueq[2];
					feq=w[a]*dens*(1.+3.*u+4.5*u*u-1.5*v);
					fneq=fintpl[a]-feq;
					Fb(V.slf,i,j,k,a,xdf,ydf,zdf,19) = feq + (dxr/dxf)*(tau/tauf)*fneq ;
				}
			}
		}
	}
	
// Interpolation on to coarse grid
	for (i=0;i<(xdc-1);i++)
	{
		for (j=0;j<(ydc-1);j++)
		{
			for (k=0;k<zdc;k++)
			{
				iref = (int)floor((i+0.5)*dxc/dxr - 0.5);
				jref = (int)floor((j+0.5)*dxc/dxr - 0.5);
				kref = (int)floor(0.5 + (k*GR+(zdf-GR-1-0.5))*dxf/dxr);
				
				ind[0] = iref - 1;
				ind[1] = iref;
				ind[2] = iref + 1;
				
				dx = (i+0.5)*dxc - (iref+0.5)*dxr;
				dy = (j+0.5)*dxc - (jref+0.5)*dxr;
				
				jnd[0] = jref - 1;
				jnd[1] = jref;
				jnd[2] = jref + 1;
				
				knd[0] = kref-1;
				knd[1] = kref;
				knd[2] = kref+1;
				  
				dz = (k*GR+(zdf-GR-1-0.5))*dxf - (kref-0.5)*dxr; 
				
		//		dz = (k*dxc+(zdf-1-0.5-GR)*dxf) + 0.5*dxr - knd[1]*dxr;
				
				for (a=0;a<19;a++)
				  fintpl[a] = 0.;
				
				dx = dx/dxr; dy =dy/dxr; dz=dz/dxr;
				sp_interpolate_qd(ind,jnd,knd,(extx+6),(exty+6),zd,dx,dy,dz,adr,fintpl);
				
				v = 0.;
				ux=0.;uy=0.;uz=0.;dens=0.;
				for (a=0;a<3;a++)
				  ueq[a] = 0.;
				
				for (a=0;a<19;a++)
				{
					ux += E0[a]*fintpl[a];
					uy += E1[a]*fintpl[a];
					uz += E2[a]*fintpl[a];
					dens += fintpl[a];
				}

				ueq[0]=ux/dens;
				ueq[1]=uy/dens;
				ueq[2]=uz/dens;
				
				v = ueq[0]*ueq[0]+ueq[1]*ueq[1]+ueq[2]*ueq[2];
				for (a=0;a<19;a++)
				{
					u=E0[a]*ueq[0]+E1[a]*ueq[1]+E2[a]*ueq[2];
					feq=w[a]*dens*(1.+3.*u+4.5*u*u-1.5*v);
					fneq=fintpl[a]-feq;
					Fb(V.f,i,j,k,a,xdc,ydc,zdc,19) = feq + (dxr/dxc)*(tau/tauc)*fneq ;
				}
			}
		}
	}
// Interpolation on to the upper fine grid
	for (i=0;i<(xdf-1);i++)
	{
		for (j=0;j<(ydf-1);j++)
		{
			for (k=0;k<(zdf-1);k++)
			{
				iref = (int)floor((i+0.5)*dxf/dxr - 0.5);
				jref = (int)floor((j+0.5)*dxf/dxr - 0.5);
				kref = (int)floor(0.5+ (k+GR*(zdc-3)+zdf-1-0.5)*dxf/dxr);
				
				ind[0] = iref - 1;
				ind[1] = iref;
				ind[2] = iref + 1;
				
				dx = (i+0.5)*dxf - (iref+0.5)*dxr;
				dy = (j+0.5)*dxf - (jref+0.5)*dxr;
				
				jnd[0] = jref - 1;
				jnd[1] = jref;
				jnd[2] = jref + 1;
				
				knd[0] = kref-1;
				knd[1] = kref;
				knd[2] = kref+1;
				
				dz = (k+GR*(zdc-3)+zdf-1-0.5)*dxf + 0.5*dxr - kref*dxr; 
				
		//		dz = (k+GR*(zdc-3)+zdf-1-0.5)*dxf + 0.5*dxr - knd[1]*dxr; 
				
				for (a=0;a<19;a++)
				  fintpl[a] = 0.;
				
				dx = dx/dxr; dy =dy/dxr; dz=dz/dxr;
				sp_interpolate_qd(ind,jnd,knd,(extx+6),(exty+6),zd,dx,dy,dz,adr,fintpl);
				
				v = 0.;
				ux=0.;uy=0.;uz=0.;dens=0.;
				for (a=0;a<3;a++)
				  ueq[a] = 0.;
				
				for (a=0;a<19;a++)
				{
					ux += E0[a]*fintpl[a];
					uy += E1[a]*fintpl[a];
					uz += E2[a]*fintpl[a];
					dens += fintpl[a];
				}

				ueq[0]=ux/dens;
				ueq[1]=uy/dens;
				ueq[2]=uz/dens;
				
				v = ueq[0]*ueq[0]+ueq[1]*ueq[1]+ueq[2]*ueq[2];
				for (a=0;a<19;a++)
				{
					u=E0[a]*ueq[0]+E1[a]*ueq[1]+E2[a]*ueq[2];
					feq=w[a]*dens*(1.+3.*u+4.5*u*u-1.5*v);
					fneq=fintpl[a]-feq;
					Fb(V.suf,i,j,k,a,xdf,ydf,zdf,19) = feq + (dxr/dxf)*(tau/tauf)*fneq ;
				}
			}
		}
	}
	
	MPI_Sendrecv(V.suf,1,sFN_send_to_dl,pd.dl,37,V.suf,1,sFN_recv_fr_ur,pd.ur,37,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(V.suf,1,sFN_send_to_dw,pd.dw,38,V.suf,1,sFN_recv_fr_up,pd.up,38,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(V.suf,1,sFN_send_to_le,pd.left,39,V.suf,1,sFN_recv_fr_ri,pd.right,39,MPI_COMM_WORLD,&tstatus);
	
	MPI_Sendrecv(V.f,1,FN_send_to_dl,pd.dl,40,V.f,1,FN_recv_fr_ur,pd.ur,40,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(V.f,1,FN_send_to_dw,pd.dw,41,V.f,1,FN_recv_fr_up,pd.up,41,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(V.f,1,FN_send_to_le,pd.left,42,V.f,1,FN_recv_fr_ri,pd.right,42,MPI_COMM_WORLD,&tstatus);
	
	MPI_Sendrecv(V.slf,1,sFN_send_to_dl,pd.dl,43,V.slf,1,sFN_recv_fr_ur,pd.ur,43,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(V.slf,1,sFN_send_to_dw,pd.dw,44,V.slf,1,sFN_recv_fr_up,pd.up,44,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(V.slf,1,sFN_send_to_le,pd.left,45,V.slf,1,sFN_recv_fr_ri,pd.right,45,MPI_COMM_WORLD,&tstatus);
  
	free(adr);
	
	MPI_Type_free(&distf);
	MPI_Type_free(&block);
	MPI_Type_free(&send_to_dw);
	MPI_Type_free(&recv_fr_dw);
	MPI_Type_free(&send_to_up);
	MPI_Type_free(&recv_fr_up);
	
	MPI_Type_free(&send_to_le);
	MPI_Type_free(&recv_fr_le);
	MPI_Type_free(&send_to_ri);
	MPI_Type_free(&recv_fr_ri);
	
	MPI_Type_free(&send_to_dl);
	MPI_Type_free(&recv_fr_dl);
	MPI_Type_free(&send_to_dr);
	MPI_Type_free(&recv_fr_dr);
	MPI_Type_free(&send_to_ul);
	MPI_Type_free(&recv_fr_ul);
	MPI_Type_free(&send_to_ur);
	MPI_Type_free(&recv_fr_ur);
	
	MPI_Type_free(&FN_send_to_dl);
	MPI_Type_free(&FN_send_to_dw);
	MPI_Type_free(&FN_send_to_le);
	MPI_Type_free(&FN_recv_fr_ur);
	MPI_Type_free(&FN_recv_fr_ri);
	MPI_Type_free(&FN_recv_fr_up);
	
	MPI_Type_free(&sFN_send_to_dl);
	MPI_Type_free(&sFN_send_to_dw);
	MPI_Type_free(&sFN_send_to_le);
	MPI_Type_free(&sFN_recv_fr_ur);
	MPI_Type_free(&sFN_recv_fr_ri);
	MPI_Type_free(&sFN_recv_fr_up);
	
  return V;
}