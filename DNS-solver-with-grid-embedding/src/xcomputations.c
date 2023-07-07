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

/*Performing Collision and streaming Parallel cut planes normal to x direction */
POINTER xcomputations(int xdc, int ydc, int zdc, int GR, POINTER V,PDATA pd, DP Gx, DP tauc, DP ubulk)
{
	int i,j,k,a,t;
	const int x[4]={0,1,(xdc-2),(xdc-1)},y[4]={0,1,(ydc-2),(ydc-1)};
	int ito,jto,kto;
	int zt[4]={0,1,(zdc-2),(zdc-1)},num=3,ct;
	const int L[5]={2,8,9,13,14};
	const int R[5]={1,7,10,11,12};
	const int U[5]={3,7,8,15,16};
	const int D[5]={4,9,10,17,18};
	const int UR[9]={1,3,5,6,7,11,12,15,16};
	const int DL[9]={2,4,5,6,9,13,14,17,18};
	const int UL[9]={2,3,5,6,8,13,14,15,16};
	const int DR[9]={1,4,5,6,10,11,12,17,18};
	DP v;
	DP feq,fneq,u,z,dens,*ptr,*adr;
	int solid,sld1,sld2,sld;
	DP ux,uy,uz,ueq[3];
	int ii,jj,kk,il,jl,kl,ic,jc,kc,xc,yc,zc;
	DP Fi,*taddrs;
	DP *tmp;
	DP tpsend[zdc][9],tprecv[zdc][9],tprecv2[zdc][9];
/***********************************************************************/
	MPI_Datatype block,send_to_dw,recv_fr_up;
	MPI_Datatype send_to_up,recv_fr_dw,send_to_ri,recv_fr_ri,send_to_le,recv_fr_le;
	int gsizes[2],lsizes[2],psizes[2],array_start[2],memsizes[2];
	MPI_Datatype send_to_ul,send_to_dr,send_to_ur,send_to_dl;
	MPI_Datatype recv_fr_ul,recv_fr_dr,recv_fr_ur,recv_fr_dl;
	int tag,flag;
	MPI_Status status[12]={0},tstatus;
/***********************************************************************/

	Gx=V.force[0];
	for (t=0;t<4;t++)
	{
		i=x[t];
		for (j=0;j<ydc;j++)
		{
			for (k=0;k<zdc;k++)
			{
				solid = Fs(V.s,i,j,k,xdc,ydc,zdc);
				if (solid>-1)
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
					
					for (a=0;a<19;a++)
					{
						u=E0[a]*ueq[0]+E1[a]*ueq[1]+E2[a]*ueq[2];
						feq=w[a]*dens*(1.+3.*u+4.5*u*u-1.5*v);
						//Fi=w[a]*(1-0.5*tauc)*(3.*(E0[a]-ueq[0])+9.*u*E0[a])*Gx;
						//Fi=BI[a]*dens*(E0[a]*Gx);
//						Fi=w[a]*(3.*(E0[a]-ueq[0])+9.*u*E0[a])*Gx;
//						Fi=w[a]*dens*(3.*(E0[a]-ueq[0])+9.*u*E0[a])*Gx;
						Fi=BI[a]*dens*(E0[a]*Gx);
						//Fi=E0[a]*Gx/12.;
						//Fi=BI[a]*E0[a]*Gx;
						//Fi=0.;
						Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19) += (Fi-tauc*(Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19)-feq));
					}
					
				}
			}
		}
	}
	for (i=2;i<(xdc-2);i++)
	{
		for (t=0;t<4;t++)
		{
			j=y[t];
			for (k=0;k<zdc;k++)
			{
				solid = Fs(V.s,i,j,k,xdc,ydc,zdc);
				if (solid>-1)
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
					
					for (a=0;a<19;a++)
					{
						u=E0[a]*ueq[0]+E1[a]*ueq[1]+E2[a]*ueq[2];
						feq=w[a]*dens*(1.+3.*u+4.5*u*u-1.5*v);
						//Fi=w[a]*(1-0.5*tauc)*(3.*(E0[a]-ueq[0])+9.*u*E0[a])*Gx;
						//Fi=BI[a]*dens*(E0[a]*Gx);
//						Fi=w[a]*(3.*(E0[a]-ueq[0])+9.*u*E0[a])*Gx;
//						Fi=w[a]*dens*(3.*(E0[a]-ueq[0])+9.*u*E0[a])*Gx;
						Fi=BI[a]*dens*(E0[a]*Gx);
						//Fi=E0[a]*Gx/12.;
						//Fi=BI[a]*E0[a]*Gx;
						//Fi=0.;
						Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19) += (Fi-tauc*(Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19)-feq));
					}
					
				}
			}
		}
	}
	for (t=0;t<4;t++)
	{
		i=x[t];
		for (j=0;j<ydc;j++)
		{
			for (k=0;k<zdc;k++)
			{
				solid = Fs(V.s,i,j,k,xdc,ydc,zdc);
				if (solid==0)
				{
					for (a=0;a<19;a++)
					{
						ito=i+e0[a];
						jto=j+e1[a];
						kto=k+e2[a];
						Fb(V.f,ito,jto,kto,a,xdc,ydc,zdc,19)=Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);
					}
				}
				else if (solid==1)
				{
					for (a=0;a<19;a++)
					{
						jto=j+e1[a];
						
						kto=k+e2[a];
						
						ito=i+e0[a];
						
						if ((ito>-1)&&(ito<xdc)&&(jto>-1)&&(jto<ydc)&&(kto>-1)&&(kto<zdc))
							Fb(V.f,ito,jto,kto,a,xdc,ydc,zdc,19)=Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);
					}
				}
			}
		}
	}
	for (i=2;i<(xdc-2);i++)
	{
		for (t=0;t<4;t++)
		{
			j=y[t];
			for (k=0;k<zdc;k++)
			{
				solid = Fs(V.s,i,j,k,xdc,ydc,zdc);
				if (solid==0)
				{
					for (a=0;a<19;a++)
					{
						ito=i+e0[a];
						jto=j+e1[a];
						kto=k+e2[a];
						Fb(V.f,ito,jto,kto,a,xdc,ydc,zdc,19)=Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);
					}
				}
				else if (solid==1)
				{
					for (a=0;a<19;a++)
					{
						jto=j+e1[a];
						
						kto=k+e2[a];
						
						ito=i+e0[a];
						
						if ((ito>-1)&&(ito<xdc)&&(jto>-1)&&(jto<ydc)&&(kto>-1)&&(kto<zdc))
							Fb(V.f,ito,jto,kto,a,xdc,ydc,zdc,19)=Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);
					}
				}
			}
		}
	}

	i=0;
	ptr=V.leftbufs;
	for (j=0;j<ydc;j++)
	{
		for (k=0;k<zdc;k++)
		{
			for (a=0;a<5;a++)
			{
				*(ptr++)=Fb(V.f,i,j,k,L[a],xdc,ydc,zdc,19);
			}
		}
	}

	i=xdc-1;
	adr=V.rightbufs;
	for (j=0;j<ydc;j++)
	{
		for (k=0;k<zdc;k++)
		{
			for (a=0;a<5;a++)
			{
				(*(adr++))=Fb(V.f,i,j,k,R[a],xdc,ydc,zdc,19);
			}
		}
	}
	
	j=0;
	ptr=V.dwbufs;
	for (i=0;i<xdc;i++)
	{
		for (k=0;k<zdc;k++)
		{
			for (a=0;a<5;a++)
			{
				*(ptr++)=Fb(V.f,i,j,k,D[a],xdc,ydc,zdc,19);
			}
		}
	}

	j=ydc-1;
	adr=V.upbufs;
	for (i=0;i<xdc;i++)
	{
		for (k=0;k<zdc;k++)
		{
			for (a=0;a<5;a++)
			{
				(*(adr++))=Fb(V.f,i,j,k,U[a],xdc,ydc,zdc,19);
			}
		}
	}
	i=0;j=0;
	adr=V.dlbufs;
	for (k=0;k<zdc;k++)
	{
		for (a=0;a<9;a++)
			(*(adr++))=Fb(V.f,i,j,k,DL[a],xdc,ydc,zdc,19);
	}
	i=xdc-1;j=ydc-1;
	adr=V.urbufs;
	for (k=0;k<zdc;k++)
	{
		for (a=0;a<9;a++)
			(*(adr++))=Fb(V.f,i,j,k,UR[a],xdc,ydc,zdc,19);
	}
	
	mpi_errors = MPI_Startall(12,&V.req[0]);
// 	if (mpi_errors != MPI_SUCCESS) 
// 	{
// 		MPI_Error_class(mpi_errors,&mpi_errors_class);
// 		MPI_Error_string(mpi_errors_class,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error class %s reported by rank %d\n",mpi_errors_message,pd.myrank);
// 		MPI_Error_string(mpi_errors,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error on rank %d: %s by MPI_Startall in line %d of function %s\n",pd.myrank,mpi_errors_message,__LINE__,__FUNCTION__);
// 		MPI_Abort(MPI_COMM_WORLD,mpi_errors);	/* abort*/
// 	}
	
	il=(xdc-4)/csc;
	jl=(ydc-4)/csc;
	kl=zdc/csc;

	xc=(xdc-4)%csc;
	yc=(ydc-4)%csc;
	zc=zdc%csc;
	
	for (ii=0;ii<=il;ii++)
	{
		for (jj=0;jj<=jl;jj++)
		{
			for (kk=0;kk<=kl;kk++)
			{
				for (ic=0;(ic<csc)&&((i=ic+csc*ii+2)<(xdc-2));ic++)
				{
					for (jc=0;(jc<csc)&&((j=jc+csc*jj+2)<(ydc-2));jc++)
					{
						for (kc=0;(kc<csc)&&((k=kc+csc*kk)<zdc);kc++)
						{
							solid = Fs(V.s,i,j,k,xdc,ydc,zdc);
							if (solid>-1)
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
								
								for (a=0;a<19;a++)
								{
									u=E0[a]*ueq[0]+E1[a]*ueq[1]+E2[a]*ueq[2];
									feq=w[a]*dens*(1.+3.*u+4.5*u*u-1.5*v);
									//Fi=w[a]*(1-0.5*tauc)*(3.*(E0[a]-ueq[0])+9.*u*E0[a])*Gx;
									//Fi=w[a]*dens*(1-0.5*tauc)*(3.*(E0[a]-ueq[0])+9.*u*E0[a])*Gx;
//									Fi=w[a]*(3.*(E0[a]-ueq[0])+9.*u*E0[a])*Gx;
//									Fi=w[a]*dens*(3.*(E0[a]-ueq[0])+9.*u*E0[a])*Gx;
									Fi=BI[a]*dens*(E0[a]*Gx);
									//Fi=BI[a]*dens*(E0[a]*Gx);
									//Fi=3.*dens*w[a]*E0[a]*Gx;
									//Fi=E0[a]*Gx/12.;
									//Fi=BI[a]*E0[a]*Gx;
									//Fi=0.;
									Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19) += (Fi-tauc*(Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19)-feq));
								}
							}
						}
					}
				}
			}
		}
	}
	mpi_errors = MPI_Testall(12,&V.req[0],&flag,&status[0]);
// 	if (mpi_errors != MPI_SUCCESS) 
// 	{
// 		MPI_Error_class(mpi_errors,&mpi_errors_class);
// 		MPI_Error_string(mpi_errors_class,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error class %s reported by rank %d\n",mpi_errors_message,pd.myrank);
// 		MPI_Error_string(mpi_errors,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error on rank %d: %s by MPI_Testall in line %d of function %s\n",pd.myrank,mpi_errors_message,__LINE__,__FUNCTION__);
// 		MPI_Abort(MPI_COMM_WORLD,mpi_errors);	/* abort*/
// 	}
	for (i=2;i<(xdc-2);i++)
	{
		for (j=2;j<(ydc-2);j++)
		{
			for (t=0;t<4;t++)
			{
				k=zt[t];
				solid = Fs(V.s,i,j,k,xdc,ydc,zdc);
				if (solid==0)
				{
					for (a=0;a<19;a++)
					{
						ito=i+e0[a];
						jto=j+e1[a];
						kto=k+e2[a];

						Fb(V.f,ito,jto,kto,a,xdc,ydc,zdc,19)=Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);
					}
				}
				else if (solid==1)
				{
					for (a=0;a<19;a++)
					{
						jto=j+e1[a];
						
						kto=k+e2[a];
						
						ito=i+e0[a];
						
						if ((kto>-1)&&(kto<zdc))
							Fb(V.f,ito,jto,kto,a,xdc,ydc,zdc,19)=Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);
					}
				}
			}
		}
	}
	for (ii=0;ii<=il;ii++)
	{
		for (jj=0;jj<=jl;jj++)
		{
			for (kk=0;kk<=kl;kk++)
			{
				for (ic=0;(ic<csc)&&((i=ic+csc*ii+2)<(xdc-2));ic++)
				{
					for (jc=0;(jc<csc)&&((j=jc+csc*jj+2)<(ydc-2));jc++)
					{
						for (kc=0;(kc<csc)&&((k=kc+csc*kk+2)<(zdc-2));kc++)
						{
							for (a=0;a<19;a++)
							{
								ito=i+e0[a];
								jto=j+e1[a];
								jto=j+e1[a];
								if (jto<0)
									jto=ydc-1;
								else if (jto>(ydc-1))
									jto=0;
								kto=k+e2[a];
								
								Fb(V.f,ito,jto,kto,a,xdc,ydc,zdc,19)=Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);
							}	
						}
					}
				}
			}
			MPI_Testall(12,&V.req[0],&flag,&status[0]);
		}
	}

	mpi_errors = MPI_Waitall(12,&V.req[0],&status[0]);
// 	if (mpi_errors != MPI_SUCCESS) 
// 	{
// 		MPI_Error_class(mpi_errors,&mpi_errors_class);
// 		MPI_Error_string(mpi_errors_class,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error class %s reported by rank %d\n",mpi_errors_message,pd.myrank);
// 		MPI_Error_string(mpi_errors,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error on rank %d: %s by MPI_Waitall in line %d of function %s\n",pd.myrank,mpi_errors_message,__LINE__,__FUNCTION__);
// 		MPI_Abort(MPI_COMM_WORLD,mpi_errors);	/* abort*/
// 	}

	i=0;
	adr=V.leftbufr;
	for (j=0;j<ydc;j++)
	{
		for (k=0;k<zdc;k++)
		{
			for (a=0;a<5;a++)
			{
				Fb(V.f,i,j,k,R[a],xdc,ydc,zdc,19)=(*(adr++));
			}
		}
	}

	ptr=V.rightbufr;
	i=xdc-1;
	for (j=0;j<ydc;j++)
	{
		for (k=0;k<zdc;k++)
		{
			for (a=0;a<5;a++)
			{
				Fb(V.f,i,j,k,L[a],xdc,ydc,zdc,19)=(*(ptr++));
			}
		}
	}
	
	j=0;
	adr=V.dwbufr;
	for (i=0;i<xdc;i++)
	{
		for (k=0;k<zdc;k++)
		{
			for (a=0;a<5;a++)
			{
				Fb(V.f,i,j,k,U[a],xdc,ydc,zdc,19)=(*(adr++));
			}
		}
	}

	ptr=V.upbufr;
	j=ydc-1;
	for (i=0;i<xdc;i++)
	{
		for (k=0;k<zdc;k++)
		{
			for (a=0;a<5;a++)
			{
				Fb(V.f,i,j,k,D[a],xdc,ydc,zdc,19)=(*(ptr++));
			}
		}
	}
	i=0;j=0;
	ptr=V.dlbufr;
	for (k=0;k<zdc;k++)
		for (a=0;a<9;a++)
			Fb(V.f,i,j,k,UR[a],xdc,ydc,zdc,19)=(*(ptr++));

	i=xdc-1;j=ydc-1;
	ptr=V.urbufr;
	for (k=0;k<zdc;k++)
	{
		for (a=0;a<9;a++)
			Fb(V.f,i,j,k,DL[a],xdc,ydc,zdc,19)=(*(ptr++));
	}
	
	i=0;j=0;
	for (k=0;k<zdc;k++)
	{
		for (a=0;a<9;a++)
			tpsend[k][a] = Fb(V.f,i,j,k,UR[a],xdc,ydc,zdc,19);
	}
	
	MPI_Sendrecv(tpsend,9*zdc,MPI_DOUBLE,pd.left,221,tprecv,9*zdc,MPI_DOUBLE,pd.right,221,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(tpsend,9*zdc,MPI_DOUBLE,pd.dw,222,tprecv2,9*zdc,MPI_DOUBLE,pd.up,222,MPI_COMM_WORLD,&tstatus);
	
	i=xdc-1;j=0;
	for (k=0;k<zdc;k++)
	{
		for (a=0;a<9;a++)
			Fb(V.f,i,j,k,UR[a],xdc,ydc,zdc,19) = tprecv[k][a];
	}
	i=0;j=ydc-1;
	for (k=0;k<zdc;k++)
	{
		for (a=0;a<9;a++)
			Fb(V.f,i,j,k,UR[a],xdc,ydc,zdc,19) = tprecv2[k][a];
	}
	i=0;j=ydc-1;
	for (k=0;k<zdc;k++)
	{
		for (a=0;a<9;a++)
			tpsend[k][a] = Fb(V.f,i,j,k,UL[a],xdc,ydc,zdc,19);
	}
	MPI_Sendrecv(tpsend,9*zdc,MPI_DOUBLE,pd.ul,223,tprecv,9*zdc,MPI_DOUBLE,pd.dr,223,MPI_COMM_WORLD,&tstatus);
	i=xdc-1;j=0;
	for (k=0;k<zdc;k++)
	{
		for (a=0;a<9;a++)
			Fb(V.f,i,j,k,UL[a],xdc,ydc,zdc,19) = tprecv[k][a];
	}
	
	i=xdc-1;j=0;
	for (k=0;k<zdc;k++)
	{
		for (a=0;a<9;a++)
			tpsend[k][a] = Fb(V.f,i,j,k,DR[a],xdc,ydc,zdc,19);
	}
	MPI_Sendrecv(tpsend,9*zdc,MPI_DOUBLE,pd.dr,224,tprecv,9*zdc,MPI_DOUBLE,pd.ul,224,MPI_COMM_WORLD,&tstatus);
	i=0;j=ydc-1;
	for (k=0;k<zdc;k++)
	{
		for (a=0;a<9;a++)
			Fb(V.f,i,j,k,DR[a],xdc,ydc,zdc,19) = tprecv[k][a];
	}
		
	tmp = V.linpp;
	V.linpp = V.linp;
	V.linp = V.lin;
	V.lin = tmp;
	
	tmp = V.uinpp;
	V.uinpp = V.uinp;
	V.uinp = V.uin;
	V.uin = tmp;
	
	k=1;
	for (i=0;i<(xdc-1);i++)
//	for (i=0;i<xdc;i++)
	{
		for (j=0;j<(ydc-1);j++)
		{
			v = 0.;
			ux=0.;uy=0.;uz=0.;dens=0.;
			for (a=0;a<19;a++)
			{
				ux += E0[a]*Fb(V.f,i,j,k,a,xdc,ydc,zdc,19);
				uy += E1[a]*Fb(V.f,i,j,k,a,xdc,ydc,zdc,19);
				uz += E2[a]*Fb(V.f,i,j,k,a,xdc,ydc,zdc,19);
				dens += Fb(V.f,i,j,k,a,xdc,ydc,zdc,19);
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
				fneq = Fb(V.f,i,j,k,a,xdc,ydc,zdc,19) - feq;
				Fbi(V.lin,i,j,0,(a+4),xdc,ydc,4,23) = fneq;
			}
		}
	}
	k=zdc-2;
	for (i=0;i<(xdc-1);i++)
//	for (i=0;i<xdc;i++)
	{
		for (j=0;j<(ydc-1);j++)
		{
			v = 0.;
			ux=0.;uy=0.;uz=0.;dens=0.;
			for (a=0;a<19;a++)
			{
				ux += E0[a]*Fb(V.f,i,j,k,a,xdc,ydc,zdc,19);
				uy += E1[a]*Fb(V.f,i,j,k,a,xdc,ydc,zdc,19);
				uz += E2[a]*Fb(V.f,i,j,k,a,xdc,ydc,zdc,19);
				dens += Fb(V.f,i,j,k,a,xdc,ydc,zdc,19);
			}
			ueq[0]=ux/dens;
			ueq[1]=uy/dens;
			ueq[2]=uz/dens;
			v = ueq[0]*ueq[0]+ueq[1]*ueq[1]+ueq[2]*ueq[2];
			
			Fbi(V.uin,i,j,0,0,xdc,ydc,4,23) = dens;
			Fbi(V.uin,i,j,0,1,xdc,ydc,4,23) = ueq[0];
			Fbi(V.uin,i,j,0,2,xdc,ydc,4,23) = ueq[1];
			Fbi(V.uin,i,j,0,3,xdc,ydc,4,23) = ueq[2];
			
			for (a=0;a<19;a++)
			{
				u=E0[a]*ueq[0]+E1[a]*ueq[1]+E2[a]*ueq[2];
				feq=w[a]*dens*(1.+3.*u+4.5*u*u-1.5*v);
				fneq = Fb(V.f,i,j,k,a,xdc,ydc,zdc,19) - feq;
				Fbi(V.uin,i,j,0,(a+4),xdc,ydc,4,23) = fneq;
			}
		}
	}
	
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

	mpi_errors = MPI_Sendrecv(V.uin,1,send_to_ri,pd.right,11,V.uin,1,recv_fr_le,pd.left,11,MPI_COMM_WORLD,&tstatus);
// 	if (mpi_errors != MPI_SUCCESS) 
// 	{
// 		MPI_Error_class(mpi_errors,&mpi_errors_class);
// 		MPI_Error_string(mpi_errors_class,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error class %s reported by rank %d\n",mpi_errors_message,pd.myrank);
// 		MPI_Error_string(mpi_errors,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error on rank %d: %s by MPI_Sendrecv in line %d of function %s\n",pd.myrank,mpi_errors_message,__LINE__,__FUNCTION__);
// 		MPI_Abort(MPI_COMM_WORLD,mpi_errors);	/* abort*/
// 	}
	
	mpi_errors = MPI_Sendrecv(V.uin,1,send_to_le,pd.left,12,V.uin,1,recv_fr_ri,pd.right,12,MPI_COMM_WORLD,&tstatus);
// 	if (mpi_errors != MPI_SUCCESS) 
// 	{
// 		MPI_Error_class(mpi_errors,&mpi_errors_class);
// 		MPI_Error_string(mpi_errors_class,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error class %s reported by rank %d\n",mpi_errors_message,pd.myrank);
// 		MPI_Error_string(mpi_errors,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error on rank %d: %s by MPI_Sendrecv in line %d of function %s\n",pd.myrank,mpi_errors_message,__LINE__,__FUNCTION__);
// 		MPI_Abort(MPI_COMM_WORLD,mpi_errors);	/* abort*/
// 	}
	
	mpi_errors = MPI_Sendrecv(V.uin,1,send_to_dw,pd.dw,13,V.uin,1,recv_fr_up,pd.up,13,MPI_COMM_WORLD,&tstatus);
// 	if (mpi_errors != MPI_SUCCESS) 
// 	{
// 		MPI_Error_class(mpi_errors,&mpi_errors_class);
// 		MPI_Error_string(mpi_errors_class,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error class %s reported by rank %d\n",mpi_errors_message,pd.myrank);
// 		MPI_Error_string(mpi_errors,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error on rank %d: %s by MPI_Sendrecv in line %d of function %s\n",pd.myrank,mpi_errors_message,__LINE__,__FUNCTION__);
// 		MPI_Abort(MPI_COMM_WORLD,mpi_errors);	/* abort*/
// 	}
	
	mpi_errors = MPI_Sendrecv(V.uin,1,send_to_up,pd.up,14,V.uin,1,recv_fr_dw,pd.dw,14,MPI_COMM_WORLD,&tstatus);
// 	if (mpi_errors != MPI_SUCCESS) 
// 	{
// 		MPI_Error_class(mpi_errors,&mpi_errors_class);
// 		MPI_Error_string(mpi_errors_class,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error class %s reported by rank %d\n",mpi_errors_message,pd.myrank);
// 		MPI_Error_string(mpi_errors,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error on rank %d: %s by MPI_Sendrecv in line %d of function %s\n",pd.myrank,mpi_errors_message,__LINE__,__FUNCTION__);
// 		MPI_Abort(MPI_COMM_WORLD,mpi_errors);	/* abort*/
// 	}
	
	mpi_errors = MPI_Sendrecv(V.uin,1,send_to_ur,pd.ur,15,V.uin,1,recv_fr_dl,pd.dl,15,MPI_COMM_WORLD,&tstatus);
// 	if (mpi_errors != MPI_SUCCESS) 
// 	{
// 		MPI_Error_class(mpi_errors,&mpi_errors_class);
// 		MPI_Error_string(mpi_errors_class,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error class %s reported by rank %d\n",mpi_errors_message,pd.myrank);
// 		MPI_Error_string(mpi_errors,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error on rank %d: %s by MPI_Sendrecv in line %d of function %s\n",pd.myrank,mpi_errors_message,__LINE__,__FUNCTION__);
// 		MPI_Abort(MPI_COMM_WORLD,mpi_errors);	/* abort*/
// 	}
	
	mpi_errors = MPI_Sendrecv(V.uin,1,send_to_ul,pd.ul,16,V.uin,1,recv_fr_dr,pd.dr,16,MPI_COMM_WORLD,&tstatus);
// 	if (mpi_errors != MPI_SUCCESS) 
// 	{
// 		MPI_Error_class(mpi_errors,&mpi_errors_class);
// 		MPI_Error_string(mpi_errors_class,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error class %s reported by rank %d\n",mpi_errors_message,pd.myrank);
// 		MPI_Error_string(mpi_errors,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error on rank %d: %s by MPI_Sendrecv in line %d of function %s\n",pd.myrank,mpi_errors_message,__LINE__,__FUNCTION__);
// 		MPI_Abort(MPI_COMM_WORLD,mpi_errors);	/* abort*/
// 	}
	
	mpi_errors = MPI_Sendrecv(V.uin,1,send_to_dl,pd.dl,17,V.uin,1,recv_fr_ur,pd.ur,17,MPI_COMM_WORLD,&tstatus);
// 	if (mpi_errors != MPI_SUCCESS) 
// 	{
// 		MPI_Error_class(mpi_errors,&mpi_errors_class);
// 		MPI_Error_string(mpi_errors_class,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error class %s reported by rank %d\n",mpi_errors_message,pd.myrank);
// 		MPI_Error_string(mpi_errors,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error on rank %d: %s by MPI_Sendrecv in line %d of function %s\n",pd.myrank,mpi_errors_message,__LINE__,__FUNCTION__);
// 		MPI_Abort(MPI_COMM_WORLD,mpi_errors);	/* abort*/
// 	}
	
	mpi_errors = MPI_Sendrecv(V.uin,1,send_to_dr,pd.dr,18,V.uin,1,recv_fr_ul,pd.ul,18,MPI_COMM_WORLD,&tstatus);
// 	if (mpi_errors != MPI_SUCCESS) 
// 	{
// 		MPI_Error_class(mpi_errors,&mpi_errors_class);
// 		MPI_Error_string(mpi_errors_class,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error class %s reported by rank %d\n",mpi_errors_message,pd.myrank);
// 		MPI_Error_string(mpi_errors,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error on rank %d: %s by MPI_Sendrecv in line %d of function %s\n",pd.myrank,mpi_errors_message,__LINE__,__FUNCTION__);
// 		MPI_Abort(MPI_COMM_WORLD,mpi_errors);	/* abort*/
// 	}
	
	mpi_errors = MPI_Sendrecv(V.lin,1,send_to_ri,pd.right,11,V.lin,1,recv_fr_le,pd.left,11,MPI_COMM_WORLD,&tstatus);
// 	if (mpi_errors != MPI_SUCCESS) 
// 	{
// 		MPI_Error_class(mpi_errors,&mpi_errors_class);
// 		MPI_Error_string(mpi_errors_class,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error class %s reported by rank %d\n",mpi_errors_message,pd.myrank);
// 		MPI_Error_string(mpi_errors,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error on rank %d: %s by MPI_Sendrecv in line %d of function %s\n",pd.myrank,mpi_errors_message,__LINE__,__FUNCTION__);
// 		MPI_Abort(MPI_COMM_WORLD,mpi_errors);	/* abort*/
// 	}
	
	mpi_errors = MPI_Sendrecv(V.lin,1,send_to_le,pd.left,12,V.lin,1,recv_fr_ri,pd.right,12,MPI_COMM_WORLD,&tstatus);
// 	if (mpi_errors != MPI_SUCCESS) 
// 	{
// 		MPI_Error_class(mpi_errors,&mpi_errors_class);
// 		MPI_Error_string(mpi_errors_class,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error class %s reported by rank %d\n",mpi_errors_message,pd.myrank);
// 		MPI_Error_string(mpi_errors,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error on rank %d: %s by MPI_Sendrecv in line %d of function %s\n",pd.myrank,mpi_errors_message,__LINE__,__FUNCTION__);
// 		MPI_Abort(MPI_COMM_WORLD,mpi_errors);	/* abort*/
// 	}
	
	mpi_errors = MPI_Sendrecv(V.lin,1,send_to_dw,pd.dw,13,V.lin,1,recv_fr_up,pd.up,13,MPI_COMM_WORLD,&tstatus);
// 	if (mpi_errors != MPI_SUCCESS) 
// 	{
// 		MPI_Error_class(mpi_errors,&mpi_errors_class);
// 		MPI_Error_string(mpi_errors_class,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error class %s reported by rank %d\n",mpi_errors_message,pd.myrank);
// 		MPI_Error_string(mpi_errors,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error on rank %d: %s by MPI_Sendrecv in line %d of function %s\n",pd.myrank,mpi_errors_message,__LINE__,__FUNCTION__);
// 		MPI_Abort(MPI_COMM_WORLD,mpi_errors);	/* abort*/
// 	}
	
	mpi_errors = MPI_Sendrecv(V.lin,1,send_to_up,pd.up,14,V.lin,1,recv_fr_dw,pd.dw,14,MPI_COMM_WORLD,&tstatus);
// 	if (mpi_errors != MPI_SUCCESS) 
// 	{
// 		MPI_Error_class(mpi_errors,&mpi_errors_class);
// 		MPI_Error_string(mpi_errors_class,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error class %s reported by rank %d\n",mpi_errors_message,pd.myrank);
// 		MPI_Error_string(mpi_errors,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error on rank %d: %s by MPI_Sendrecv in line %d of function %s\n",pd.myrank,mpi_errors_message,__LINE__,__FUNCTION__);
// 		MPI_Abort(MPI_COMM_WORLD,mpi_errors);	/* abort*/
// 	}
	
	mpi_errors = MPI_Sendrecv(V.lin,1,send_to_ur,pd.ur,15,V.lin,1,recv_fr_dl,pd.dl,15,MPI_COMM_WORLD,&tstatus);
// 	if (mpi_errors != MPI_SUCCESS) 
// 	{
// 		MPI_Error_class(mpi_errors,&mpi_errors_class);
// 		MPI_Error_string(mpi_errors_class,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error class %s reported by rank %d\n",mpi_errors_message,pd.myrank);
// 		MPI_Error_string(mpi_errors,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error on rank %d: %s by MPI_Sendrecv in line %d of function %s\n",pd.myrank,mpi_errors_message,__LINE__,__FUNCTION__);
// 		MPI_Abort(MPI_COMM_WORLD,mpi_errors);	/* abort*/
// 	}
	
	mpi_errors = MPI_Sendrecv(V.lin,1,send_to_ul,pd.ul,16,V.lin,1,recv_fr_dr,pd.dr,16,MPI_COMM_WORLD,&tstatus);
// 	if (mpi_errors != MPI_SUCCESS) 
// 	{
// 		MPI_Error_class(mpi_errors,&mpi_errors_class);
// 		MPI_Error_string(mpi_errors_class,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error class %s reported by rank %d\n",mpi_errors_message,pd.myrank);
// 		MPI_Error_string(mpi_errors,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error on rank %d: %s by MPI_Sendrecv in line %d of function %s\n",pd.myrank,mpi_errors_message,__LINE__,__FUNCTION__);
// 		MPI_Abort(MPI_COMM_WORLD,mpi_errors);	/* abort*/
// 	}
	
	mpi_errors = MPI_Sendrecv(V.lin,1,send_to_dl,pd.dl,17,V.lin,1,recv_fr_ur,pd.ur,17,MPI_COMM_WORLD,&tstatus);
// 	if (mpi_errors != MPI_SUCCESS) 
// 	{
// 		MPI_Error_class(mpi_errors,&mpi_errors_class);
// 		MPI_Error_string(mpi_errors_class,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error class %s reported by rank %d\n",mpi_errors_message,pd.myrank);
// 		MPI_Error_string(mpi_errors,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error on rank %d: %s by MPI_Sendrecv in line %d of function %s\n",pd.myrank,mpi_errors_message,__LINE__,__FUNCTION__);
// 		MPI_Abort(MPI_COMM_WORLD,mpi_errors);	/* abort*/
// 	}
	
	mpi_errors = MPI_Sendrecv(V.lin,1,send_to_dr,pd.dr,18,V.lin,1,recv_fr_ul,pd.ul,18,MPI_COMM_WORLD,&tstatus);
// 	if (mpi_errors != MPI_SUCCESS) 
// 	{
// 		MPI_Error_class(mpi_errors,&mpi_errors_class);
// 		MPI_Error_string(mpi_errors_class,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error class %s reported by rank %d\n",mpi_errors_message,pd.myrank);
// 		MPI_Error_string(mpi_errors,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error on rank %d: %s by MPI_Sendrecv in line %d of function %s\n",pd.myrank,mpi_errors_message,__LINE__,__FUNCTION__);
// 		MPI_Abort(MPI_COMM_WORLD,mpi_errors);	/* abort*/
// 	}
	
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
	
//	fprintf(stderr,"in function xcomputations, c error handler is %s\n",strerror(errno));

return V;
}