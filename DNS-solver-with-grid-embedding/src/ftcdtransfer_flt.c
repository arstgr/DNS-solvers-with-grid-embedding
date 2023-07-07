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

/* Data transfere from fine to coarse grid with volume averaging (approximate box filter) */
POINTER ftcdtransfer_flt(int xdf, int ydf, int zdf, int xdc, int ydc, int zdc, int itr, double tauc, double tauf, int GR, POINTER V, PDATA pd)
{
	int i,j,k,kk,a,q,m,n;
	double feq,u,v,ueq[3],fneq,fneqt[19];
	double ux,uy,uz,dens;
	double fneqf[xdf+2][ydf+2][3][23]; /* Note: must be initialized */
	int ic, jc, kc;
	
	double sigmas=GR/2.;
	
/***********************************************************************/
	MPI_Datatype block,send_to_dw,recv_fr_up;
	MPI_Datatype send_to_up,recv_fr_dw,send_to_ri,recv_fr_ri,send_to_le,recv_fr_le;
	int gsizes[2],lsizes[2],psizes[2],array_start[2],memsizes[2];
	MPI_Datatype send_to_ul,send_to_dr,send_to_ur,send_to_dl;
	MPI_Datatype recv_fr_ul,recv_fr_dr,recv_fr_ur,recv_fr_dl;
	int tag,flag;
	MPI_Status tstatus;
/***********************************************************************/

/****************************************************************************************************************/
	MPI_Type_contiguous(23*3,MPI_DOUBLE,&block);
	MPI_Type_commit(&block);
	
	memsizes[0] = xdf+2;
	memsizes[1] = ydf+2;
/****************************************************************************************************************/
	lsizes[0] = xdf;
	lsizes[1] = 1;
	array_start[0] = 1;
	array_start[1] = ydf-1;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &send_to_up);
	MPI_Type_commit(&send_to_up);
	
	lsizes[0] = xdf;
	lsizes[1] = 1;
	array_start[0] = 1;
	array_start[1] = 0;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &recv_fr_dw);
	MPI_Type_commit(&recv_fr_dw);
	
	lsizes[0] = xdf;
	lsizes[1] = 1;
	array_start[0] = 1;
	array_start[1] = 2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &send_to_dw);
	MPI_Type_commit(&send_to_dw);
	
	lsizes[0] = xdf;
	lsizes[1] = 1;
	array_start[0] = 1;
	array_start[1] = ydf+1;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &recv_fr_up);
	MPI_Type_commit(&recv_fr_up);
	
/********************************************************************************************************************************/
	lsizes[0] = 1;
	lsizes[1] = ydf;
	array_start[0] = xdf-1;
	array_start[1] = 1;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &send_to_ri);
	MPI_Type_commit(&send_to_ri);
	
	lsizes[0] = 1;
	lsizes[1] = ydf;
	array_start[0] = 0;
	array_start[1] = 1;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &recv_fr_le);
	MPI_Type_commit(&recv_fr_le);
	
	lsizes[0] = 1;
	lsizes[1] = ydf;
	array_start[0] = 2;
	array_start[1] = 1;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &send_to_le);
	MPI_Type_commit(&send_to_le);
	
	lsizes[0] = 1;
	lsizes[1] = ydf;
	array_start[0] = xdf+1;
	array_start[1] = 1;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &recv_fr_ri);
	MPI_Type_commit(&recv_fr_ri);
	
/*********************************************************************************************************************************/
	lsizes[0] = 1;
	lsizes[1] = 1;
	array_start[0] = xdf-1;
	array_start[1] = ydf-1;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &send_to_ur);
	MPI_Type_commit(&send_to_ur);
	
	lsizes[0] = 1;
	lsizes[1] = 1;
	array_start[0] = xdf+1;
	array_start[1] = ydf+1;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &recv_fr_ur);
	MPI_Type_commit(&recv_fr_ur);
	
	lsizes[0] = 1;
	lsizes[1] = 1;
	array_start[0] = 2;
	array_start[1] = ydf-1;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &send_to_ul);
	MPI_Type_commit(&send_to_ul);
	
	lsizes[0] = 1;
	lsizes[1] = 1;
	array_start[0] = 0;
	array_start[1] = ydf+1;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &recv_fr_ul);
	MPI_Type_commit(&recv_fr_ul);
	
	lsizes[0] = 1;
	lsizes[1] = 1;
	array_start[0] = xdf-1;
	array_start[1] = 2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &send_to_dr);
	MPI_Type_commit(&send_to_dr);
	
	lsizes[0] = 1;
	lsizes[1] = 1;
	array_start[0] = xdf+1;
	array_start[1] = 0;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &recv_fr_dr);
	MPI_Type_commit(&recv_fr_dr);
	
	lsizes[0] = 1;
	lsizes[1] = 1;
	array_start[0] = 2;
	array_start[1] = 2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &send_to_dl);
	MPI_Type_commit(&send_to_dl);
	
	lsizes[0] = 1;
	lsizes[1] = 1;
	array_start[0] = 0;
	array_start[1] = 0;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &recv_fr_dl);
	MPI_Type_commit(&recv_fr_dl);
	
/*****************************************************************************************************************************************/
	/* calculating the fneq on the lower grid */
	for (i=0;i<xdf;i++)
	{
		for (j=0;j<ydf;j++)
		{
			for (kk=-1;kk<2;kk++)
			{
				k= zdf-1-GR + kk;
				
				v = 0.;
				ux=0.;uy=0.;uz=0.;dens=0.;
				for (a=0;a<19;a++)
				{
					ux += E0[a]*Fb(V.slf,i,j,k,a,xdf,ydf,zdf,19);
					uy += E1[a]*Fb(V.slf,i,j,k,a,xdf,ydf,zdf,19);
					uz += E2[a]*Fb(V.slf,i,j,k,a,xdf,ydf,zdf,19);
					dens += Fb(V.slf,i,j,k,a,xdf,ydf,zdf,19);
				}
				ueq[0]=ux/dens;
				ueq[1]=uy/dens;
				ueq[2]=uz/dens;
				v = ueq[0]*ueq[0]+ueq[1]*ueq[1]+ueq[2]*ueq[2];
				
				fneqf[i+1][j+1][kk+1][0] = dens;
				fneqf[i+1][j+1][kk+1][1] = ux;
				fneqf[i+1][j+1][kk+1][2] = uy;
				fneqf[i+1][j+1][kk+1][3] = uz;
				
				for (a=0;a<19;a++)
				{
					u=E0[a]*ueq[0]+E1[a]*ueq[1]+E2[a]*ueq[2];
					feq=w[a]*dens*(1.+3.*u+4.5*u*u-1.5*v);
					fneqt[a]=(Fb(V.slf,i,j,k,a,xdf,ydf,zdf,19)-feq);
					//fneqf[i+1][j+1][kk+1][a+4]=fneqt;
				}
				memcpy(&fneqf[i+1][j+1][kk+1][4],&fneqt[0],19*sizeof(double));
			}
		}
	}

	mpi_errors = MPI_Sendrecv(&fneqf[0][0][0][0],1,send_to_ri,pd.right,11,&fneqf[0][0][0][0],1,recv_fr_le,pd.left,11,MPI_COMM_WORLD,&tstatus);
	mpi_errors = MPI_Sendrecv(&fneqf[0][0][0][0],1,send_to_le,pd.left,12,&fneqf[0][0][0][0],1,recv_fr_ri,pd.right,12,MPI_COMM_WORLD,&tstatus);
	mpi_errors = MPI_Sendrecv(&fneqf[0][0][0][0],1,send_to_dw,pd.dw,13,&fneqf[0][0][0][0],1,recv_fr_up,pd.up,13,MPI_COMM_WORLD,&tstatus);
	mpi_errors = MPI_Sendrecv(&fneqf[0][0][0][0],1,send_to_up,pd.up,14,&fneqf[0][0][0][0],1,recv_fr_dw,pd.dw,14,MPI_COMM_WORLD,&tstatus);
	mpi_errors = MPI_Sendrecv(&fneqf[0][0][0][0],1,send_to_ur,pd.ur,15,&fneqf[0][0][0][0],1,recv_fr_dl,pd.dl,15,MPI_COMM_WORLD,&tstatus);
	mpi_errors = MPI_Sendrecv(&fneqf[0][0][0][0],1,send_to_ul,pd.ul,16,&fneqf[0][0][0][0],1,recv_fr_dr,pd.dr,16,MPI_COMM_WORLD,&tstatus);
	mpi_errors = MPI_Sendrecv(&fneqf[0][0][0][0],1,send_to_dl,pd.dl,17,&fneqf[0][0][0][0],1,recv_fr_ur,pd.ur,17,MPI_COMM_WORLD,&tstatus);
	mpi_errors = MPI_Sendrecv(&fneqf[0][0][0][0],1,send_to_dr,pd.dr,18,&fneqf[0][0][0][0],1,recv_fr_ul,pd.ul,18,MPI_COMM_WORLD,&tstatus);
	
	for (ic=0;ic<xdc;ic++)
	{
		for (jc=0;jc<ydc;jc++)
		{
			i=ic*GR;
			j=jc*GR;
			k= zdf-1-GR;
			
			if (ic==(xdc-1))
			  i=xdf-1;
			if (jc==(ydc-1))
			  j=ydf-1;
			
			dens = fneqf[i+1][j+1][1][0];
			ueq[0] = fneqf[i+1][j+1][1][1];
			ueq[1] = fneqf[i+1][j+1][1][2];
			ueq[2] = fneqf[i+1][j+1][1][3];
			
			v = 0.; u=0.;
			v = ueq[0]*ueq[0]+ueq[1]*ueq[1]+ueq[2]*ueq[2];
					
			for (a=0;a<19;a++)
			{
				u=E0[a]*ueq[0]+E1[a]*ueq[1]+E2[a]*ueq[2];
				feq=w[a]*dens*(1.+3.*u+4.5*u*u-1.5*v);
				fneq=0.;
				
#ifdef BFLTR
				fneq = (1./7.)*fneqf[i+1][j+1][1][a+4];
				
				for (q=1;q<7;q++)
				  fneq += (1./14.)*fneqf[i+1+e0[q]][j+1+e1[q]][1+e2[q]][a+4];
				
				for (q=7;q<19;q++)
				  fneq += (1./28.)*fneqf[i+1+e0[q]][j+1+e1[q]][1+e2[q]][a+4];
#endif
#ifdef GFLTR
				fneq = (1. - 2.*sigmas*sigmas)*fneqf[i+1][j+1][1][a+4];
				
				for (q=1;q<7;q++)
				  fneq += (sigmas*sigmas/6.)*fneqf[i+1+e0[q]][j+1+e1[q]][1+e2[q]][a+4];
				
				for (q=7;q<19;q++)
				  fneq += (sigmas*sigmas/12.)*fneqf[i+1+e0[q]][j+1+e1[q]][1+e2[q]][a+4];
#endif
				Fb(V.f,ic,jc,0,a,xdc,ydc,zdc,19)= feq + GR*(tauf/tauc)*fneq ;	
			}
		}
	}
	
	/* Calculating on the upper fine grid */
	for (i=0;i<xdf;i++)
	{
		for (j=0;j<ydf;j++)
		{
			for (kk=-1;kk<2;kk++)
			{
				k= GR + kk;
				
				v = 0.;
				ux=0.;uy=0.;uz=0.;dens=0.;
				for (a=0;a<19;a++)
				{
					ux += E0[a]*Fb(V.suf,i,j,k,a,xdf,ydf,zdf,19);
					uy += E1[a]*Fb(V.suf,i,j,k,a,xdf,ydf,zdf,19);
					uz += E2[a]*Fb(V.suf,i,j,k,a,xdf,ydf,zdf,19);
					dens += Fb(V.suf,i,j,k,a,xdf,ydf,zdf,19);
				}
				ueq[0]=ux/dens;
				ueq[1]=uy/dens;
				ueq[2]=uz/dens;
				v = ueq[0]*ueq[0]+ueq[1]*ueq[1]+ueq[2]*ueq[2];
				
				fneqf[i+1][j+1][kk+1][0] = dens;
				fneqf[i+1][j+1][kk+1][1] = ux;
				fneqf[i+1][j+1][kk+1][2] = uy;
				fneqf[i+1][j+1][kk+1][3] = uz;
						
				for (a=0;a<19;a++)
				{
					u=E0[a]*ueq[0]+E1[a]*ueq[1]+E2[a]*ueq[2];
					feq=w[a]*dens*(1.+3.*u+4.5*u*u-1.5*v);
					fneqt[a] = (Fb(V.suf,i,j,k,a,xdf,ydf,zdf,19)-feq);
				//	fneqf[i+1][j+1][kk+1][a+4]=(Fb(V.suf,i,j,k,a,xdf,ydf,zdf,19)-feq);
				}
				memcpy(&fneqf[i+1][j+1][kk+1][4],&fneqt[0],19*sizeof(double));
			}
		}
	}

	mpi_errors = MPI_Sendrecv(&fneqf[0][0][0][0],1,send_to_ri,pd.right,11,&fneqf[0][0][0][0],1,recv_fr_le,pd.left,11,MPI_COMM_WORLD,&tstatus);
	mpi_errors = MPI_Sendrecv(&fneqf[0][0][0][0],1,send_to_le,pd.left,12,&fneqf[0][0][0][0],1,recv_fr_ri,pd.right,12,MPI_COMM_WORLD,&tstatus);
	mpi_errors = MPI_Sendrecv(&fneqf[0][0][0][0],1,send_to_dw,pd.dw,13,&fneqf[0][0][0][0],1,recv_fr_up,pd.up,13,MPI_COMM_WORLD,&tstatus);
	mpi_errors = MPI_Sendrecv(&fneqf[0][0][0][0],1,send_to_up,pd.up,14,&fneqf[0][0][0][0],1,recv_fr_dw,pd.dw,14,MPI_COMM_WORLD,&tstatus);
	mpi_errors = MPI_Sendrecv(&fneqf[0][0][0][0],1,send_to_ur,pd.ur,15,&fneqf[0][0][0][0],1,recv_fr_dl,pd.dl,15,MPI_COMM_WORLD,&tstatus);
	mpi_errors = MPI_Sendrecv(&fneqf[0][0][0][0],1,send_to_ul,pd.ul,16,&fneqf[0][0][0][0],1,recv_fr_dr,pd.dr,16,MPI_COMM_WORLD,&tstatus);
	mpi_errors = MPI_Sendrecv(&fneqf[0][0][0][0],1,send_to_dl,pd.dl,17,&fneqf[0][0][0][0],1,recv_fr_ur,pd.ur,17,MPI_COMM_WORLD,&tstatus);
	mpi_errors = MPI_Sendrecv(&fneqf[0][0][0][0],1,send_to_dr,pd.dr,18,&fneqf[0][0][0][0],1,recv_fr_ul,pd.ul,18,MPI_COMM_WORLD,&tstatus);
	
	for (ic=0;ic<xdc;ic++)
	{
		for (jc=0;jc<ydc;jc++)
		{
			i=ic*GR;
			j=jc*GR;
			k= GR;
			
			if (ic==(xdc-1))
			  i=xdf-1;
			if (jc==(ydc-1))
			  j=ydf-1;
			
			dens = fneqf[i+1][j+1][1][0];
			ueq[0] = fneqf[i+1][j+1][1][1];
			ueq[1] = fneqf[i+1][j+1][1][2];
			ueq[2] = fneqf[i+1][j+1][1][3];
			
			v = 0.; u=0.;
			v = ueq[0]*ueq[0]+ueq[1]*ueq[1]+ueq[2]*ueq[2];
					
			for (a=0;a<19;a++)
			{
				u=E0[a]*ueq[0]+E1[a]*ueq[1]+E2[a]*ueq[2];
				feq=w[a]*dens*(1.+3.*u+4.5*u*u-1.5*v);
				fneq=0.;
				
#ifdef BFLTR
				fneq = (1./7.)*fneqf[i+1][j+1][1][a+4];
				
				for (q=1;q<7;q++)
				  fneq += (1./14.)*fneqf[i+1+e0[q]][j+1+e1[q]][1+e2[q]][a+4];
				
				for (q=7;q<19;q++)
				  fneq += (1./28.)*fneqf[i+1+e0[q]][j+1+e1[q]][1+e2[q]][a+4];
#endif			
#ifdef GFLTR
				fneq = (1. - 2.*sigmas*sigmas)*fneqf[i+1][j+1][1][a+4];
				
				for (q=1;q<7;q++)
				  fneq += (sigmas*sigmas/6.)*fneqf[i+1+e0[q]][j+1+e1[q]][1+e2[q]][a+4];
				
				for (q=7;q<19;q++)
				  fneq += (sigmas*sigmas/12.)*fneqf[i+1+e0[q]][j+1+e1[q]][1+e2[q]][a+4];
#endif
				
				Fb(V.f,ic,jc,(zdc-1),a,xdc,ydc,zdc,19) = feq + GR*(tauf/tauc)*fneq;
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