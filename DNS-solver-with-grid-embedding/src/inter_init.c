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

POINTER inter_init(int xdc, int ydc, int zdc, POINTER V, PDATA pd)
{
	int i,j,k,a;
	DP ux,uy,uz,dens;
	DP feq,fneq,ueq[3],v,u;
	MPI_Status tstatus;
	MPI_Datatype block,send_to_dw,recv_fr_up;
	MPI_Datatype send_to_up,recv_fr_dw,send_to_ri,recv_fr_ri,send_to_le,recv_fr_le;
	int gsizes[2],lsizes[2],psizes[2],array_start[2],memsizes[2];
	MPI_Datatype send_to_ul,send_to_dr,send_to_ur,send_to_dl;
	MPI_Datatype recv_fr_ul,recv_fr_dr,recv_fr_ur,recv_fr_dl;
	
	/****************************************************************************************************************/
	MPI_Type_contiguous(23,MPI_DOUBLE,&block);
	MPI_Type_commit(&block);
	
	memsizes[0] = xdc+5;
	memsizes[1] = ydc+5;
/****************************************************************************************************************/
	
	lsizes[0] = xdc-1;
	lsizes[1] = 2;
	array_start[0] = 2;
	array_start[1] = ydc-1-2 +2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &send_to_up);
	MPI_Type_commit(&send_to_up);
	
	lsizes[0] = xdc-1;
	lsizes[1] = 2;
	array_start[0] = 2;
	array_start[1] = 0;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &recv_fr_dw);
	MPI_Type_commit(&recv_fr_dw);
	
	lsizes[0] = xdc-1;
	lsizes[1] = 4;
	array_start[0] = 2;
	array_start[1] = 2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &send_to_dw);
	MPI_Type_commit(&send_to_dw);
	
	lsizes[0] = xdc-1;
	lsizes[1] = 4;
	array_start[0] = 2;
	array_start[1] = ydc-1+2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &recv_fr_up);
	MPI_Type_commit(&recv_fr_up);
	
/********************************************************************************************************************************/
	lsizes[0] = 2;
	lsizes[1] = ydc-1;
	array_start[0] = xdc-1-2 +2;
	array_start[1] = 2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &send_to_ri);
	MPI_Type_commit(&send_to_ri);
	
	lsizes[0] = 2;
	lsizes[1] = ydc-1;
	array_start[0] = 0;
	array_start[1] = 2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &recv_fr_le);
	MPI_Type_commit(&recv_fr_le);
	
	lsizes[0] = 4;
	lsizes[1] = ydc-1;
	array_start[0] = 2;
	array_start[1] = 2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &send_to_le);
	MPI_Type_commit(&send_to_le);
	
	lsizes[0] = 4;
	lsizes[1] = ydc-1;
	array_start[0] = xdc-1 + 2;
	array_start[1] = 2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &recv_fr_ri);
	MPI_Type_commit(&recv_fr_ri);
	
/*********************************************************************************************************************************/
	lsizes[0] = 2;
	lsizes[1] = 2;
	array_start[0] = xdc-1-2 +2;
	array_start[1] = ydc-1-2 +2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &send_to_ur);
	MPI_Type_commit(&send_to_ur);
	
	lsizes[0] = 4;
	lsizes[1] = 4;
	array_start[0] = xdc-1+2;
	array_start[1] = ydc-1+2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &recv_fr_ur);
	MPI_Type_commit(&recv_fr_ur);
	
	lsizes[0] = 4;
	lsizes[1] = 2;
	array_start[0] = 2;
	array_start[1] = ydc -1 -2 +2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &send_to_ul);
	MPI_Type_commit(&send_to_ul);
	
	lsizes[0] = 2;
	lsizes[1] = 4;
	array_start[0] = 0;
	array_start[1] = ydc-1+2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &recv_fr_ul);
	MPI_Type_commit(&recv_fr_ul);
	
	lsizes[0] = 2;
	lsizes[1] = 4;
	array_start[0] = xdc-1-2 +2;
	array_start[1] = 2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &send_to_dr);
	MPI_Type_commit(&send_to_dr);
	
	lsizes[0] = 4;
	lsizes[1] = 2;
	array_start[0] = xdc -1 +2;
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
	k=1;
	for (i=0;i<(xdc-1);i++)
	{
		for (j=0;j<(ydc-1);j++)
		{
			v = 0.;
			ux=0.;uy=0.;uz=0.;dens=0.;
			for (a=0;a<19;a++)
			{
				ux += E0[a]*Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);
				uy += E1[a]*Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);
				uz += E2[a]*Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);
				dens += Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);
			}
			ueq[0]=ux/dens;
			ueq[1]=uy/dens;
			ueq[2]=uz/dens;
			v = ueq[0]*ueq[0]+ueq[1]*ueq[1]+ueq[2]*ueq[2];
			
			Fbi(V.lin,i,j,0,0,xdc,ydc,4,23) = dens;
			Fbi(V.lin,i,j,0,1,xdc,ydc,4,23) = ueq[0];
			Fbi(V.lin,i,j,0,2,xdc,ydc,4,23) = ueq[1];
			Fbi(V.lin,i,j,0,3,xdc,ydc,4,23) = ueq[2];
			
			for (a=0;a<19;a++)
			{
				u=E0[a]*ueq[0]+E1[a]*ueq[1]+E2[a]*ueq[2];
				feq=w[a]*dens*(1.+3.*u+4.5*u*u-1.5*v);
				fneq = Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19) - feq;
				Fbi(V.lin,i,j,0,(a+4),xdc,ydc,4,23) = fneq;
			}
		}
	}
	k=zdc-2;
	for (i=0;i<(xdc-1);i++)
	{
		for (j=0;j<(ydc-1);j++)
		{
			v = 0.;
			ux=0.;uy=0.;uz=0.;dens=0.;
			for (a=0;a<19;a++)
			{
				ux += E0[a]*Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);
				uy += E1[a]*Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);
				uz += E2[a]*Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);
				dens += Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);
			}
			ueq[0]=ux/dens;
			ueq[1]=uy/dens;
			ueq[2]=uz/dens;
			v = ueq[0]*ueq[0]+ueq[1]*ueq[1]+ueq[2]*ueq[2];
			
			Fbi(V.lin,i,j,0,0,xdc,ydc,4,23) = dens;
			Fbi(V.lin,i,j,0,1,xdc,ydc,4,23) = ueq[0];
			Fbi(V.lin,i,j,0,2,xdc,ydc,4,23) = ueq[1];
			Fbi(V.lin,i,j,0,3,xdc,ydc,4,23) = ueq[2];
			
			for (a=0;a<19;a++)
			{
				u=E0[a]*ueq[0]+E1[a]*ueq[1]+E2[a]*ueq[2];
				feq=w[a]*dens*(1.+3.*u+4.5*u*u-1.5*v);
				fneq = Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19) - feq;
				Fbi(V.uin,i,j,0,(a+4),xdc,ydc,4,23) = fneq;
			}
		}
	}
/*********************************************************************************************************************/
	
	MPI_Sendrecv(V.uin,1,send_to_ri,pd.right,11,V.uin,1,recv_fr_le,pd.left,11,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(V.uin,1,send_to_le,pd.left,12,V.uin,1,recv_fr_ri,pd.right,12,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(V.uin,1,send_to_dw,pd.dw,13,V.uin,1,recv_fr_up,pd.up,13,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(V.uin,1,send_to_up,pd.up,14,V.uin,1,recv_fr_dw,pd.dw,14,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(V.uin,1,send_to_ur,pd.ur,15,V.uin,1,recv_fr_dl,pd.dl,15,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(V.uin,1,send_to_ul,pd.ul,16,V.uin,1,recv_fr_dr,pd.dr,16,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(V.uin,1,send_to_dl,pd.dl,17,V.uin,1,recv_fr_ur,pd.ur,17,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(V.uin,1,send_to_dr,pd.dr,18,V.uin,1,recv_fr_ul,pd.ul,18,MPI_COMM_WORLD,&tstatus);
	
	MPI_Sendrecv(V.lin,1,send_to_ri,pd.right,11,V.lin,1,recv_fr_le,pd.left,11,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(V.lin,1,send_to_le,pd.left,12,V.lin,1,recv_fr_ri,pd.right,12,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(V.lin,1,send_to_dw,pd.dw,13,V.lin,1,recv_fr_up,pd.up,13,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(V.lin,1,send_to_up,pd.up,14,V.lin,1,recv_fr_dw,pd.dw,14,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(V.lin,1,send_to_ur,pd.ur,15,V.lin,1,recv_fr_dl,pd.dl,15,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(V.lin,1,send_to_ul,pd.ul,16,V.lin,1,recv_fr_dr,pd.dr,16,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(V.lin,1,send_to_dl,pd.dl,17,V.lin,1,recv_fr_ur,pd.ur,17,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(V.lin,1,send_to_dr,pd.dr,18,V.lin,1,recv_fr_ul,pd.ul,18,MPI_COMM_WORLD,&tstatus);
	
/*******************************************************************************************************************************************/
	for (i=-2;i<(xdc+3);i++)
	{
		for (j=-2;j<(ydc+3);j++)
		{
			for (a=0;a<19;a++)
			{
				Fbi(V.linp,i,j,0,a,xdc,ydc,4,19) = Fbi(V.lin,i,j,0,a,xdc,ydc,4,19);
				Fbi(V.linpp,i,j,0,a,xdc,ydc,4,19) = Fbi(V.lin,i,j,0,a,xdc,ydc,4,19);
				Fbi(V.uinp,i,j,0,a,xdc,ydc,4,19) = Fbi(V.uin,i,j,0,a,xdc,ydc,4,19);
				Fbi(V.uinpp,i,j,0,a,xdc,ydc,4,19) = Fbi(V.uin,i,j,0,a,xdc,ydc,4,19);
			}
		}
	}
	
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
	
	return V;
}