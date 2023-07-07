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

/* Performs collision and streaming on the fine grids */
POINTER xfcomputations(int xdf, int ydf, int zdf, double GR, POINTER V, PDATA pd, DP tauf, DP ubulk)
{
	int i,j,k,a,t;
	const int xf[4]={0,1,(xdf-2),(xdf-1)},yf[4]={0,1,(ydf-2),(ydf-1)};
	const int zt[4]={0,1,(zdf-2),(zdf-1)};
	int ito,jto,kto;
	int ct;
	const int L[5]={2,8,9,13,14};
	const int R[5]={1,7,10,11,12};
	const int U[5]={3,7,8,15,16};
	const int D[5]={4,9,10,17,18};
	const int UR[9]={1,3,5,6,7,11,12,15,16};
	const int DL[9]={2,4,5,6,9,13,14,17,18};
	const int UL[9]={2,3,5,6,8,13,14,15,16};
	const int DR[9]={1,4,5,6,10,11,12,17,18};
	int tag,flag;
	MPI_Status status[12]={0},tstatus;
	DP v,Gx;
	DP feq,u,z,dens,*ptr,*adr;
	int solid,sld1,sld2,sld;
	DP ux,uy,uz,ueq[3];
	int ii,jj,kk,il,jl,kl,ic,jc,kc,xc,yc,zc;
	DP Fi,*taddrs,*tmp;
	DP rl,rr,dq;
	DP tpsend[zdf][9],tprecv[zdf][9],tprecv2[zdf][9];

	Gx=GR*V.force[0];
	for (t=0;t<4;t++)
	{
		i=xf[t];
		for (j=0;j<ydf;j++)
		{
			for (k=0;k<zdf;k++)
			{
				solid = Fs(V.sus,i,j,k,xdf,ydf,zdf);
				if (solid>-1)
				{	
					v = 0.;
					ux=0.;uy=0.;uz=0.;dens=0.;
					for (a=0;a<19;a++)
					{
						ux += E0[a]*Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
						uy += E1[a]*Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
						uz += E2[a]*Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
						dens += Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
					}
					ueq[0]=ux/dens;
					ueq[1]=uy/dens;
					ueq[2]=uz/dens;
					v = ueq[0]*ueq[0]+ueq[1]*ueq[1]+ueq[2]*ueq[2];
					
					for (a=0;a<19;a++)
					{
						u=E0[a]*ueq[0]+E1[a]*ueq[1]+E2[a]*ueq[2];
						feq=w[a]*dens*(1.+3.*u+4.5*u*u-1.5*v);
						//Fi=w[a]*(1-0.5*tauf)*(3.*(E0[a]-ueq[0])+9.*u*E0[a])*Gx;
						//Fi=BI[a]*dens*(E0[a]*Gx+E1[a]*Gy+E2[a]*Gz);
						//Fi=BI[a]*dens*E0[a]*Gx;
//						Fi=w[a]*(3.*(E0[a]-ueq[0])+9.*u*E0[a])*Gx;
//						Fi=w[a]*dens*(3.*(E0[a]-ueq[0])+9.*u*E0[a])*Gx;
						Fi=BI[a]*dens*(E0[a]*Gx);
						//Fi=E0[a]*Gx/12.;
						//Fi=BI[a]*E0[a]*Gx;
						//Fi=0.;
						Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19) += (Fi-tauf*(Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19)-feq));
					}
				}
			}
		}
	}
	for (i=2;i<(xdf-2);i++)
	{
		for (t=0;t<4;t++)
		{
			j=yf[t];
			for (k=0;k<zdf;k++)
			{
				solid = Fs(V.sus,i,j,k,xdf,ydf,zdf);
				if (solid>-1)
				{	
					v = 0.;
					ux=0.;uy=0.;uz=0.;dens=0.;
					for (a=0;a<19;a++)
					{
						ux += E0[a]*Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
						uy += E1[a]*Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
						uz += E2[a]*Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
						dens += Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
					}
					ueq[0]=ux/dens;
					ueq[1]=uy/dens;
					ueq[2]=uz/dens;
					v = ueq[0]*ueq[0]+ueq[1]*ueq[1]+ueq[2]*ueq[2];
					
					for (a=0;a<19;a++)
					{
						u=E0[a]*ueq[0]+E1[a]*ueq[1]+E2[a]*ueq[2];
						feq=w[a]*dens*(1.+3.*u+4.5*u*u-1.5*v);
						//Fi=w[a]*(1-0.5*tauf)*(3.*(E0[a]-ueq[0])+9.*u*E0[a])*Gx;
						//Fi=BI[a]*dens*(E0[a]*Gx+E1[a]*Gy+E2[a]*Gz);
						//Fi=BI[a]*dens*E0[a]*Gx;
//						Fi=w[a]*(3.*(E0[a]-ueq[0])+9.*u*E0[a])*Gx;
//						Fi=w[a]*dens*(3.*(E0[a]-ueq[0])+9.*u*E0[a])*Gx;
						Fi=BI[a]*dens*(E0[a]*Gx);
						//Fi=E0[a]*Gx/12.;
						//Fi=BI[a]*E0[a]*Gx;
						//Fi=0.;
						Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19) += (Fi-tauf*(Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19)-feq));
					}
				}
			}
		}
	}
	for (t=0;t<4;t++)
	{
		i=xf[t];
		for (j=0;j<ydf;j++)
		{
			for (k=0;k<zdf;k++)
			{
				solid = Fs(V.sls,i,j,k,xdf,ydf,zdf);
				if (solid>-1)
				{	
					v = 0.;
					ux=0.;uy=0.;uz=0.;dens=0.;
					for (a=0;a<19;a++)
					{
						ux += E0[a]*Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
						uy += E1[a]*Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
						uz += E2[a]*Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
						dens += Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
					}
					ueq[0]=ux/dens;
					ueq[1]=uy/dens;
					ueq[2]=uz/dens;
					v = ueq[0]*ueq[0]+ueq[1]*ueq[1]+ueq[2]*ueq[2];
					
					for (a=0;a<19;a++)
					{
						u=E0[a]*ueq[0]+E1[a]*ueq[1]+E2[a]*ueq[2];
						feq=w[a]*dens*(1.+3.*u+4.5*u*u-1.5*v);
						//Fi=w[a]*(1-0.5*tauf)*(3.*(E0[a]-ueq[0])+9.*u*E0[a])*Gx;
						//Fi=BI[a]*dens*(E0[a]*Gx+E1[a]*Gy+E2[a]*Gz);
						//Fi=BI[a]*dens*E0[a]*Gx;
//						Fi=w[a]*(3.*(E0[a]-ueq[0])+9.*u*E0[a])*Gx;
//						Fi=w[a]*dens*(3.*(E0[a]-ueq[0])+9.*u*E0[a])*Gx;
						Fi=BI[a]*dens*(E0[a]*Gx);
						//Fi=E0[a]*Gx/12.;
						//Fi=BI[a]*E0[a]*Gx;
						//Fi=0.;
						Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19) += (Fi-tauf*(Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19)-feq));
					}
				}
			}
		}
	}
	for (i=2;i<(xdf-2);i++)
	{
		for (t=0;t<4;t++)
		{
			j=yf[t];
			for (k=0;k<zdf;k++)
			{
				solid = Fs(V.sls,i,j,k,xdf,ydf,zdf);
				if (solid>-1)
				{	
					v = 0.;
					ux=0.;uy=0.;uz=0.;dens=0.;
					for (a=0;a<19;a++)
					{
						ux += E0[a]*Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
						uy += E1[a]*Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
						uz += E2[a]*Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
						dens += Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
					}
					ueq[0]=ux/dens;
					ueq[1]=uy/dens;
					ueq[2]=uz/dens;
					v = ueq[0]*ueq[0]+ueq[1]*ueq[1]+ueq[2]*ueq[2];
					
					for (a=0;a<19;a++)
					{
						u=E0[a]*ueq[0]+E1[a]*ueq[1]+E2[a]*ueq[2];
						feq=w[a]*dens*(1.+3.*u+4.5*u*u-1.5*v);
						//Fi=w[a]*(1-0.5*tauf)*(3.*(E0[a]-ueq[0])+9.*u*E0[a])*Gx;
						//Fi=BI[a]*dens*(E0[a]*Gx+E1[a]*Gy+E2[a]*Gz);
						//Fi=BI[a]*dens*E0[a]*Gx;
//						Fi=w[a]*(3.*(E0[a]-ueq[0])+9.*u*E0[a])*Gx;
//						Fi=w[a]*dens*(3.*(E0[a]-ueq[0])+9.*u*E0[a])*Gx;
						Fi=BI[a]*dens*(E0[a]*Gx);
						//Fi=E0[a]*Gx/12.;
						//Fi=BI[a]*E0[a]*Gx;
						//Fi=0.;
						Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19) += (Fi-tauf*(Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19)-feq));
					}
				}
			}
		}
	}
	for (t=0;t<4;t++)
	{
		i=xf[t];
		for (j=0;j<ydf;j++)
		{
			for (k=0;k<zdf;k++)
			{
				solid = Fs(V.sus,i,j,k,xdf,ydf,zdf);
				if (solid==0)
				{
					for (a=0;a<19;a++)
					{
						ito=i+e0[a];
						jto=j+e1[a];
						kto=k+e2[a];
						sld=Fs(V.sus,ito,jto,kto,xdf,ydf,zdf);
						sld1=Fs(V.sus,i,j,kto,xdf,ydf,zdf);
						if (sld>-1)
							Fb(V.suf,ito,jto,kto,a,xdf,ydf,zdf,19)=Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
						else if ((sld==-1)&&(sld1==-1))
							Fb(V.suf,i,j,k,opposite[a],xdf,ydf,zdf,19)=Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
						else if ((sld==-1)&&(sld1==-2))
							Fb(V.suf,ito,jto,k,mirror[a],xdf,ydf,zdf,19)=Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
						else if (sld==-2)
							Fb(V.suf,ito,jto,k,mirror[a],xdf,ydf,zdf,19)=Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
					}
				}
				else if (solid==1)
				{
					for (a=0;a<19;a++)
					{
						jto=j+e1[a];
						
						kto=k+e2[a];
						
						ito=i+e0[a];
						
						if ((ito>-1)&&(ito<xdf)&&(jto>-1)&&(jto<ydf)&&(kto>-1)&&(kto<zdf))
						{
							sld=Fs(V.sus,ito,jto,kto,xdf,ydf,zdf);
							sld1=Fs(V.sus,i,j,kto,xdf,ydf,zdf);
							if (sld>-1)
								Fb(V.suf,ito,jto,kto,a,xdf,ydf,zdf,19)=Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
							else if ((sld==-1)&&(sld1==-1))
								Fb(V.suf,i,j,k,opposite[a],xdf,ydf,zdf,19)=Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
							else if ((sld==-1)&&(sld1==-2))
								Fb(V.suf,ito,jto,k,mirror[a],xdf,ydf,zdf,19)=Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
							else if (sld==-2)
								Fb(V.suf,ito,jto,k,mirror[a],xdf,ydf,zdf,19)=Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
						}
					}
				}
			}
		}
	}
	for (i=2;i<(xdf-2);i++)
	{
		for (t=0;t<4;t++)
		{
			j=yf[t];
			for (k=0;k<zdf;k++)
			{
				solid = Fs(V.sus,i,j,k,xdf,ydf,zdf);
				if (solid==0)
				{
					for (a=0;a<19;a++)
					{
						ito=i+e0[a];
						jto=j+e1[a];
						kto=k+e2[a];
						sld=Fs(V.sus,ito,jto,kto,xdf,ydf,zdf);
						sld1=Fs(V.sus,i,j,kto,xdf,ydf,zdf);
						if (sld>-1)
							Fb(V.suf,ito,jto,kto,a,xdf,ydf,zdf,19)=Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
						else if ((sld==-1)&&(sld1==-1))
							Fb(V.suf,i,j,k,opposite[a],xdf,ydf,zdf,19)=Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
						else if ((sld==-1)&&(sld1==-2))
							Fb(V.suf,ito,jto,k,mirror[a],xdf,ydf,zdf,19)=Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
						else if (sld==-2)
							Fb(V.suf,ito,jto,k,mirror[a],xdf,ydf,zdf,19)=Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
					}
				}
				else if (solid==1)
				{
					for (a=0;a<19;a++)
					{
						jto=j+e1[a];
						
						kto=k+e2[a];
						
						ito=i+e0[a];
						
						if ((ito>-1)&&(ito<xdf)&&(jto>-1)&&(jto<ydf)&&(kto>-1)&&(kto<zdf))
						{
							sld=Fs(V.sus,ito,jto,kto,xdf,ydf,zdf);
							sld1=Fs(V.sus,i,j,kto,xdf,ydf,zdf);
							if (sld>-1)
								Fb(V.suf,ito,jto,kto,a,xdf,ydf,zdf,19)=Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
							else if ((sld==-1)&&(sld1==-1))
								Fb(V.suf,i,j,k,opposite[a],xdf,ydf,zdf,19)=Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
							else if ((sld==-1)&&(sld1==-2))
								Fb(V.suf,ito,jto,k,mirror[a],xdf,ydf,zdf,19)=Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
							else if (sld==-2)
								Fb(V.suf,ito,jto,k,mirror[a],xdf,ydf,zdf,19)=Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
						}
					}
				}
			}
		}
	}
	for (t=0;t<4;t++)
	{
		i=xf[t];
		for (j=0;j<ydf;j++)
		{
			for (k=0;k<zdf;k++)
			{
				solid = Fs(V.sls,i,j,k,xdf,ydf,zdf);
				if (solid==0)
				{
					for (a=0;a<19;a++)
					{
						ito=i+e0[a];
						jto=j+e1[a];
						kto=k+e2[a];
						sld=Fs(V.sls,ito,jto,kto,xdf,ydf,zdf);
						sld1=Fs(V.sls,i,j,kto,xdf,ydf,zdf);
						if (sld>-1)
							Fb(V.slf,ito,jto,kto,a,xdf,ydf,zdf,19)=Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
						else if ((sld==-1)&&(sld1==-1))
							Fb(V.slf,i,j,k,opposite[a],xdf,ydf,zdf,19)=Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
						else if ((sld==-1)&&(sld1==-2))
							Fb(V.slf,ito,jto,k,mirror[a],xdf,ydf,zdf,19)=Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
						else if (sld==-2)
							Fb(V.slf,ito,jto,k,mirror[a],xdf,ydf,zdf,19)=Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
					}
				}
				else if (solid==1)
				{
					for (a=0;a<19;a++)
					{
						jto=j+e1[a];
						if (jto<0)
							jto=ydf-1;
						else if (jto>(ydf-1))
							jto=0;
						
						kto=k+e2[a];
						
						ito=i+e0[a];
						
						if ((ito>-1)&&(ito<xdf)&&(jto>-1)&&(jto<ydf)&&(kto>-1)&&(kto<zdf))
						{
							sld=Fs(V.sls,ito,jto,kto,xdf,ydf,zdf);
							sld1=Fs(V.sls,i,j,kto,xdf,ydf,zdf);
							if (sld>-1)
								Fb(V.slf,ito,jto,kto,a,xdf,ydf,zdf,19)=Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
							else if ((sld==-1)&&(sld1==-1))
								Fb(V.slf,i,j,k,opposite[a],xdf,ydf,zdf,19)=Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
							else if ((sld==-1)&&(sld1==-2))
								Fb(V.slf,ito,jto,k,mirror[a],xdf,ydf,zdf,19)=Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
							else if (sld==-2)
								Fb(V.slf,ito,jto,k,mirror[a],xdf,ydf,zdf,19)=Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
						}
					}
				}
			}
		}
	}
	for (i=2;i<(xdf-2);i++)
	{
		for (t=0;t<4;t++)
		{
			j=yf[t];
			for (k=0;k<zdf;k++)
			{
				solid = Fs(V.sls,i,j,k,xdf,ydf,zdf);
				if (solid==0)
				{
					for (a=0;a<19;a++)
					{
						ito=i+e0[a];
						jto=j+e1[a];
						kto=k+e2[a];
						sld=Fs(V.sls,ito,jto,kto,xdf,ydf,zdf);
						sld1=Fs(V.sls,i,j,kto,xdf,ydf,zdf);
						if (sld>-1)
							Fb(V.slf,ito,jto,kto,a,xdf,ydf,zdf,19)=Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
						else if ((sld==-1)&&(sld1==-1))
							Fb(V.slf,i,j,k,opposite[a],xdf,ydf,zdf,19)=Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
						else if ((sld==-1)&&(sld1==-2))
							Fb(V.slf,ito,jto,k,mirror[a],xdf,ydf,zdf,19)=Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
						else if (sld==-2)
							Fb(V.slf,ito,jto,k,mirror[a],xdf,ydf,zdf,19)=Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
					}
				}
				else if (solid==1)
				{
					for (a=0;a<19;a++)
					{
						jto=j+e1[a];
						if (jto<0)
							jto=ydf-1;
						else if (jto>(ydf-1))
							jto=0;
						
						kto=k+e2[a];
						
						ito=i+e0[a];
						
						if ((ito>-1)&&(ito<xdf)&&(jto>-1)&&(jto<ydf)&&(kto>-1)&&(kto<zdf))
						{
							sld=Fs(V.sls,ito,jto,kto,xdf,ydf,zdf);
							sld1=Fs(V.sls,i,j,kto,xdf,ydf,zdf);
							if (sld>-1)
								Fb(V.slf,ito,jto,kto,a,xdf,ydf,zdf,19)=Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
							else if ((sld==-1)&&(sld1==-1))
								Fb(V.slf,i,j,k,opposite[a],xdf,ydf,zdf,19)=Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
							else if ((sld==-1)&&(sld1==-2))
								Fb(V.slf,ito,jto,k,mirror[a],xdf,ydf,zdf,19)=Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
							else if (sld==-2)
								Fb(V.slf,ito,jto,k,mirror[a],xdf,ydf,zdf,19)=Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
						}
					}
				}
			}
		}
	}
	
	i=0;
	ptr=V.suleftbufs;
	for (j=0;j<ydf;j++)
	{
		for (k=0;k<zdf;k++)
		{
			for (a=0;a<5;a++)
			{
				*(ptr++)=Fb(V.suf,i,j,k,L[a],xdf,ydf,zdf,19);
			}
		}
	}
	i=0;
	ptr=V.slleftbufs;
	for (j=0;j<ydf;j++)
	{
		for (k=0;k<zdf;k++)
		{
			for (a=0;a<5;a++)
			{
				*(ptr++)=Fb(V.slf,i,j,k,L[a],xdf,ydf,zdf,19);
			}
		}
	}
	j=0;
	ptr=V.sudwbufs;
	for (i=0;i<xdf;i++)
	{
		for (k=0;k<zdf;k++)
		{
			for (a=0;a<5;a++)
			{
				*(ptr++)=Fb(V.suf,i,j,k,D[a],xdf,ydf,zdf,19);
			}
		}
	}
	j=0;
	ptr=V.sldwbufs;
	for (i=0;i<xdf;i++)
	{
		for (k=0;k<zdf;k++)
		{
			for (a=0;a<5;a++)
			{
				*(ptr++)=Fb(V.slf,i,j,k,D[a],xdf,ydf,zdf,19);
			}
		}
	}

	i=xdf-1;
	adr=V.surightbufs;
	for (j=0;j<ydf;j++)
	{
		for (k=0;k<zdf;k++)
		{
			for (a=0;a<5;a++)
			{
				(*(adr++))=Fb(V.suf,i,j,k,R[a],xdf,ydf,zdf,19);
			}
		}
	}
	i=xdf-1;
	adr=V.slrightbufs;
	for (j=0;j<ydf;j++)
	{
		for (k=0;k<zdf;k++)
		{
			for (a=0;a<5;a++)
			{
				(*(adr++))=Fb(V.slf,i,j,k,R[a],xdf,ydf,zdf,19);
			}
		}
	}
	j=ydf-1;
	adr=V.suupbufs;
	for (i=0;i<xdf;i++)
	{
		for (k=0;k<zdf;k++)
		{
			for (a=0;a<5;a++)
			{
				(*(adr++))=Fb(V.suf,i,j,k,U[a],xdf,ydf,zdf,19);
			}
		}
	}
	j=ydf-1;
	adr=V.slupbufs;
	for (i=0;i<xdf;i++)
	{
		for (k=0;k<zdf;k++)
		{
			for (a=0;a<5;a++)
			{
				(*(adr++))=Fb(V.slf,i,j,k,U[a],xdf,ydf,zdf,19);
			}
		}
	}
		
	i=0;j=0;
	adr=V.sldlbufs;
	for (k=0;k<zdf;k++)
	{
		for (a=0;a<9;a++)
			(*(adr++))=Fb(V.slf,i,j,k,DL[a],xdf,ydf,zdf,19);
	}
	adr=V.sudlbufs;
	for (k=0;k<zdf;k++)
	{
		for (a=0;a<9;a++)
			(*(adr++))=Fb(V.suf,i,j,k,DL[a],xdf,ydf,zdf,19);
	}
	
	i=xdf-1;j=ydf-1;
	adr=V.slurbufs;
	for (k=0;k<zdf;k++)
	{
		for (a=0;a<9;a++)
			(*(adr++))=Fb(V.slf,i,j,k,UR[a],xdf,ydf,zdf,19);
	}
	adr=V.suurbufs;
	for (k=0;k<zdf;k++)
	{
		for (a=0;a<9;a++)
			(*(adr++))=Fb(V.suf,i,j,k,UR[a],xdf,ydf,zdf,19);
	}
	
	mpi_errors = MPI_Startall(12,&V.sureq[0]);
// 	if (mpi_errors != MPI_SUCCESS) 
// 	{
// 		MPI_Error_class(mpi_errors,&mpi_errors_class);
// 		MPI_Error_string(mpi_errors_class,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error class %s reported by rank %d\n",mpi_errors_message,pd.myrank);
// 		MPI_Error_string(mpi_errors,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error on rank %d: %s by MPI_Startall in xcomputations, in line %d of function %s\n",pd.myrank,mpi_errors_message,__LINE__,__FUNCTION__);
// 		MPI_Abort(MPI_COMM_WORLD,mpi_errors);	/* abort*/
// 	}
	mpi_errors = MPI_Startall(12,&V.slreq[0]);
// 	if (mpi_errors != MPI_SUCCESS) 
// 	{
// 		MPI_Error_class(mpi_errors,&mpi_errors_class);
// 		MPI_Error_string(mpi_errors_class,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error class %s reported by rank %d\n",mpi_errors_message,pd.myrank);
// 		MPI_Error_string(mpi_errors,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error on rank %d: %s by MPI_Startall in xcomputations, in line %d of function %s\n",pd.myrank,mpi_errors_message,__LINE__,__FUNCTION__);
// 		MPI_Abort(MPI_COMM_WORLD,mpi_errors);	/* abort*/
// 	}
	
	il=(xdf-4)/csf;
	jl=(ydf-4)/csf;
	kl=zdf/csf;

	xc=(xdf-4)%csf;
	yc=(ydf-4)%csf;
	zc=zdf%csf;

	for (ii=0;ii<=il;ii++)
	{
		for (jj=0;jj<=jl;jj++)
		{
			for (kk=0;kk<=kl;kk++)
			{
				for (ic=0;(ic<csf)&&((i=ic+csf*ii+2)<(xdf-2));ic++)
				{
					for (jc=0;(jc<csf)&&((j=jc+csf*jj+2)<(ydf-2));jc++)
					{
						for (kc=0;(kc<csf)&&((k=kc+csf*kk)<zdf);kc++)
						{
							solid = Fs(V.sus,i,j,k,xdf,ydf,zdf);
							if (solid>-1)
							{	
								v = 0.;
								ux=0.;uy=0.;uz=0.;dens=0.;
								for (a=0;a<19;a++)
								{
									ux += E0[a]*Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
									uy += E1[a]*Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
									uz += E2[a]*Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
									dens += Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
								}
								ueq[0]=ux/dens;
								ueq[1]=uy/dens;
								ueq[2]=uz/dens;
								v = ueq[0]*ueq[0]+ueq[1]*ueq[1]+ueq[2]*ueq[2];
								
								for (a=0;a<19;a++)
								{
									u=E0[a]*ueq[0]+E1[a]*ueq[1]+E2[a]*ueq[2];
									feq=w[a]*dens*(1.+3.*u+4.5*u*u-1.5*v);
									//Fi=w[a]*(1-0.5*tauf)*(3.*(E0[a]-ueq[0])+9.*u*E0[a])*Gx;
									//Fi=BI[a]*dens*(E0[a]*Gx+E1[a]*Gy+E2[a]*Gz);
									//Fi=BI[a]*dens*E0[a]*Gx;
									//Fi=3.*dens*w[a]*E0[a]*Gx;
//									Fi=w[a]*(3.*(E0[a]-ueq[0])+9.*u*E0[a])*Gx;
//									Fi=w[a]*dens*(3.*(E0[a]-ueq[0])+9.*u*E0[a])*Gx;
									Fi=BI[a]*dens*(E0[a]*Gx);
									//Fi=E0[a]*Gx/12.;
									//Fi=BI[a]*E0[a]*Gx;
									//Fi=0.;
									Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19) += (Fi-tauf*(Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19)-feq));
								}
							}
						}
					}
				}
			}
		}
	}
	mpi_errors = MPI_Testall(12,&V.sureq[0],&flag,&status[0]);
// 	if (mpi_errors != MPI_SUCCESS) 
// 	{
// 		MPI_Error_class(mpi_errors,&mpi_errors_class);
// 		MPI_Error_string(mpi_errors_class,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error class %s reported by rank %d\n",mpi_errors_message,pd.myrank);
// 		MPI_Error_string(mpi_errors,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error on rank %d: %s by MPI_Testall in xcomputations, in line %d of function %s\n",pd.myrank,mpi_errors_message,__LINE__,__FUNCTION__);
// 		MPI_Abort(MPI_COMM_WORLD,mpi_errors);	/* abort*/
// 	}
	for (ii=0;ii<=il;ii++)
	{
		for (jj=0;jj<=jl;jj++)
		{
			for (kk=0;kk<=kl;kk++)
			{
				for (ic=0;(ic<csf)&&((i=ic+csf*ii+2)<(xdf-2));ic++)
				{
					for (jc=0;(jc<csf)&&((j=jc+csf*jj+2)<(ydf-2));jc++)
					{
						for (kc=0;(kc<csf)&&((k=kc+csf*kk)<zdf);kc++)
						{
							solid = Fs(V.sls,i,j,k,xdf,ydf,zdf);
							if (solid>-1)
							{	
								v = 0.;
								ux=0.;uy=0.;uz=0.;dens=0.;
								for (a=0;a<19;a++)
								{
									ux += E0[a]*Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
									uy += E1[a]*Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
									uz += E2[a]*Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
									dens += Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
								}
								ueq[0]=ux/dens;
								ueq[1]=uy/dens;
								ueq[2]=uz/dens;
								v = ueq[0]*ueq[0]+ueq[1]*ueq[1]+ueq[2]*ueq[2];
								
								for (a=0;a<19;a++)
								{
									u=E0[a]*ueq[0]+E1[a]*ueq[1]+E2[a]*ueq[2];
									feq=w[a]*dens*(1.+3.*u+4.5*u*u-1.5*v);
									//Fi=w[a]*(1-0.5*tauf)*(3.*(E0[a]-ueq[0])+9.*u*E0[a])*Gx;
									//Fi=BI[a]*dens*(E0[a]*Gx+E1[a]*Gy+E2[a]*Gz);
									//Fi=BI[a]*dens*E0[a]*Gx;
//									Fi=w[a]*(3.*(E0[a]-ueq[0])+9.*u*E0[a])*Gx;
//									Fi=w[a]*dens*(3.*(E0[a]-ueq[0])+9.*u*E0[a])*Gx;
									Fi=BI[a]*dens*(E0[a]*Gx);
									//Fi=3.*dens*w[a]*E0[a]*Gx;
									//Fi=E0[a]*Gx/12.;
									//Fi=BI[a]*E0[a]*Gx;
									//Fi=0.;
									Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19) += (Fi-tauf*(Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19)-feq));
								}
							}
						}
					}
				}
			}
		}
	}
	mpi_errors = MPI_Testall(12,&V.slreq[0],&flag,&status[0]);
// 	if (mpi_errors != MPI_SUCCESS) 
// 	{
// 		MPI_Error_class(mpi_errors,&mpi_errors_class);
// 		MPI_Error_string(mpi_errors_class,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error class %s reported by rank %d\n",mpi_errors_message,pd.myrank);
// 		MPI_Error_string(mpi_errors,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error on rank %d: %s by MPI_Testall in xcomputations, in li//ne %d of function %s\n",pd.myrank,mpi_errors_message,__LINE__,__FUNCTION__);
// 		MPI_Abort(MPI_COMM_WORLD,mpi_errors);	/* abort*/
// 	}
	
	for (i=2;i<(xdf-2);i++)
	{
		for (j=2;j<(ydf-2);j++)
		{
			for (t=0;t<4;t++)
			{
				k=zt[t];
				solid = Fs(V.sus,i,j,k,xdf,ydf,zdf);
				if (solid==0)
				{
					for (a=0;a<19;a++)
					{
						ito=i+e0[a];
						jto=j+e1[a];
						kto=k+e2[a];
						sld=Fs(V.sus,ito,jto,kto,xdf,ydf,zdf);
						sld1=Fs(V.sus,i,j,kto,xdf,ydf,zdf);
						if (sld>-1)
							Fb(V.suf,ito,jto,kto,a,xdf,ydf,zdf,19)=Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
						else if ((sld==-1)&&(sld1==-1))
							Fb(V.suf,i,j,k,opposite[a],xdf,ydf,zdf,19)=Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
						else if ((sld==-1)&&(sld1==-2))
							Fb(V.suf,ito,jto,k,mirror[a],xdf,ydf,zdf,19)=Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
						else if (sld==-2)
							Fb(V.suf,ito,jto,k,mirror[a],xdf,ydf,zdf,19)=Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
					}
				}
				else if (solid==1)
				{
					for (a=0;a<19;a++)
					{
						jto=j+e1[a];
						
						kto=k+e2[a];
						
						ito=i+e0[a];
						
						if ((kto>-1)&&(kto<zdf))
						{
							sld=Fs(V.sus,ito,jto,kto,xdf,ydf,zdf);
							sld1=Fs(V.sus,i,j,kto,xdf,ydf,zdf);
							if (sld>-1)
								Fb(V.suf,ito,jto,kto,a,xdf,ydf,zdf,19)=Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
							else if ((sld==-1)&&(sld1==-1))
								Fb(V.suf,i,j,k,opposite[a],xdf,ydf,zdf,19)=Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
							else if ((sld==-1)&&(sld1==-2))
								Fb(V.suf,ito,jto,k,mirror[a],xdf,ydf,zdf,19)=Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
							else if (sld==-2)
								Fb(V.suf,ito,jto,k,mirror[a],xdf,ydf,zdf,19)=Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
						}
					}
				}
			}
		}
	}
	
	for (i=2;i<(xdf-2);i++)
	{
		for (j=2;j<(ydf-2);j++)
		{
			for (t=0;t<4;t++)
			{
				k=zt[t];
				solid = Fs(V.sls,i,j,k,xdf,ydf,zdf);
				if (solid==0)
				{
					for (a=0;a<19;a++)
					{
						ito=i+e0[a];
						jto=j+e1[a];
						kto=k+e2[a];
						sld=Fs(V.sls,ito,jto,kto,xdf,ydf,zdf);
						sld1=Fs(V.sls,i,j,kto,xdf,ydf,zdf);
						if (sld>-1)
							Fb(V.slf,ito,jto,kto,a,xdf,ydf,zdf,19)=Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
						else if ((sld==-1)&&(sld1==-1))
							Fb(V.slf,i,j,k,opposite[a],xdf,ydf,zdf,19)=Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
						else if ((sld==-1)&&(sld1==-2))
							Fb(V.slf,ito,jto,k,mirror[a],xdf,ydf,zdf,19)=Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
						else if (sld==-2)
							Fb(V.slf,ito,jto,k,mirror[a],xdf,ydf,zdf,19)=Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
					}
				}
				else if (solid==1)
				{
					for (a=0;a<19;a++)
					{
						jto=j+e1[a];
						
						kto=k+e2[a];
						
						ito=i+e0[a];
						
						if ((kto>-1)&&(kto<zdf))
						{
							sld=Fs(V.sls,ito,jto,kto,xdf,ydf,zdf);
							sld1=Fs(V.sls,i,j,kto,xdf,ydf,zdf);
							if (sld>-1)
								Fb(V.slf,ito,jto,kto,a,xdf,ydf,zdf,19)=Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
							else if ((sld==-1)&&(sld1==-1))
								Fb(V.slf,i,j,k,opposite[a],xdf,ydf,zdf,19)=Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
							else if ((sld==-1)&&(sld1==-2))
								Fb(V.slf,ito,jto,k,mirror[a],xdf,ydf,zdf,19)=Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
							else if (sld==-2)
								Fb(V.slf,ito,jto,k,mirror[a],xdf,ydf,zdf,19)=Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
						}
					}
				}
			}
		}
	}

	il=(xdf-4)/csf;
	jl=(ydf-4)/csf;
	kl=zdf/csf;

	xc=(xdf-4)%csf;
	yc=(ydf-4)%csf;
	zc=zdf%csf;
	for (ii=0;ii<=il;ii++)
	{
		for (jj=0;jj<=jl;jj++)
		{
			for (kk=0;kk<=kl;kk++)
			{
				for (ic=0;(ic<csf)&&((i=ic+csf*ii+2)<(xdf-2));ic++)
				{
					for (jc=0;(jc<csf)&&((j=jc+csf*jj+2)<(ydf-2));jc++)
					{
						for (kc=0;(kc<csf)&&((k=kc+csf*kk+2)<(zdf-2));kc++)
						{
							for (a=0;a<19;a++)
							{
								ito=i+e0[a];
								jto=j+e1[a];
								kto=k+e2[a];
								
								Fb(V.suf,ito,jto,kto,a,xdf,ydf,zdf,19)=Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
							}
						}
					}
				}
			}
			MPI_Testall(12,&V.sureq[0],&flag,&status[0]);
		}
	}

	for (ii=0;ii<=il;ii++)
	{
		for (jj=0;jj<=jl;jj++)
		{
			for (kk=0;kk<=kl;kk++)
			{
				for (ic=0;(ic<csf)&&((i=ic+csf*ii+2)<(xdf-2));ic++)
				{
					for (jc=0;(jc<csf)&&((j=jc+csf*jj+2)<(ydf-2));jc++)
					{
						for (kc=0;(kc<csf)&&((k=kc+csf*kk+2)<(zdf-2));kc++)
						{
							for (a=0;a<19;a++)
							{
								ito=i+e0[a];
								jto=j+e1[a];
								kto=k+e2[a];
								
								Fb(V.slf,ito,jto,kto,a,xdf,ydf,zdf,19)=Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
							}
						}
					}
				}
			}
			MPI_Testall(12,&V.slreq[0],&flag,&status[0]);
		}
	}
	
	mpi_errors = MPI_Waitall(12,&V.sureq[0],&status[0]);
// 	if (mpi_errors != MPI_SUCCESS) 
// 	{
// 		MPI_Error_class(mpi_errors,&mpi_errors_class);
// 		MPI_Error_string(mpi_errors_class,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error class %s reported by rank %d\n",mpi_errors_message,pd.myrank);
// 		MPI_Error_string(mpi_errors,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error on rank %d: %s by MPI_Waitall in xcomputations, in line %d of function %s\n",pd.myrank,mpi_errors_message,__LINE__,__FUNCTION__);
// 		MPI_Abort(MPI_COMM_WORLD,mpi_errors);	/* abort*/
// 	}
	mpi_errors = MPI_Waitall(12,&V.slreq[0],&status[0]);
// 	if (mpi_errors != MPI_SUCCESS) 
// 	{
// 		MPI_Error_class(mpi_errors,&mpi_errors_class);
// 		MPI_Error_string(mpi_errors_class,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error class %s reported by rank %d\n",mpi_errors_message,pd.myrank);
// 		MPI_Error_string(mpi_errors,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error on rank %d: %s by MPI_Testall in xcomputations, in line %d of function %s\n",pd.myrank,mpi_errors_message,__LINE__,__FUNCTION__);
// 		MPI_Abort(MPI_COMM_WORLD,mpi_errors);	/* abort*/
// 	}

	i=0;
	adr=V.suleftbufr;
	for (j=0;j<ydf;j++)
	{
		for (k=0;k<zdf;k++)
		{
			for (a=0;a<5;a++)
			{
				Fb(V.suf,i,j,k,R[a],xdf,ydf,zdf,19)=(*(adr++));
			}
		}
	}
	i=0;
	adr=V.slleftbufr;
	for (j=0;j<ydf;j++)
	{
		for (k=0;k<zdf;k++)
		{
			for (a=0;a<5;a++)
			{
				Fb(V.slf,i,j,k,R[a],xdf,ydf,zdf,19)=(*(adr++));
			}
		}
	}
	j=0;
	adr=V.sudwbufr;
	for (i=0;i<xdf;i++)
	{
		for (k=0;k<zdf;k++)
		{
			for (a=0;a<5;a++)
			{
				Fb(V.suf,i,j,k,U[a],xdf,ydf,zdf,19)=(*(adr++));
			}
		}
	}
	j=0;
	adr=V.sldwbufr;
	for (i=0;i<xdf;i++)
	{
		for (k=0;k<zdf;k++)
		{
			for (a=0;a<5;a++)
			{
				Fb(V.slf,i,j,k,U[a],xdf,ydf,zdf,19)=(*(adr++));
			}
		}
	}

	ptr=V.surightbufr;
	i=xdf-1;
	for (j=0;j<ydf;j++)
	{
		for (k=0;k<zdf;k++)
		{
			for (a=0;a<5;a++)
			{
				Fb(V.suf,i,j,k,L[a],xdf,ydf,zdf,19)=(*(ptr++));
			}
		}
	}
	ptr=V.slrightbufr;
	i=xdf-1;
	for (j=0;j<ydf;j++)
	{
		for (k=0;k<zdf;k++)
		{
			for (a=0;a<5;a++)
			{
				Fb(V.slf,i,j,k,L[a],xdf,ydf,zdf,19)=(*(ptr++));
			}
		}
	}
	ptr=V.suupbufr;
	j=ydf-1;
	for (i=0;i<xdf;i++)
	{
		for (k=0;k<zdf;k++)
		{
			for (a=0;a<5;a++)
			{
				Fb(V.suf,i,j,k,D[a],xdf,ydf,zdf,19)=(*(ptr++));
			}
		}
	}
	ptr=V.slupbufr;
	j=ydf-1;
	for (i=0;i<xdf;i++)
	{
		for (k=0;k<zdf;k++)
		{
			for (a=0;a<5;a++)
			{
				Fb(V.slf,i,j,k,D[a],xdf,ydf,zdf,19)=(*(ptr++));
			}
		}
	}
	i=0;j=0;
	ptr=V.sldlbufr;
	for (k=0;k<zdf;k++)
	{
		for (a=0;a<9;a++)
			Fb(V.slf,i,j,k,UR[a],xdf,ydf,zdf,19)=(*(ptr++));
	}
	ptr=V.sudlbufr;
	for (k=0;k<zdf;k++)
	{
		for (a=0;a<9;a++)
			Fb(V.suf,i,j,k,UR[a],xdf,ydf,zdf,19)=(*(ptr++));
	}
	
	i=xdf-1;j=ydf-1;
	ptr=V.slurbufr;
	for (k=0;k<zdf;k++)
	{
		for (a=0;a<9;a++)
			Fb(V.slf,i,j,k,DL[a],xdf,ydf,zdf,19)=(*(ptr++));
	}
	ptr=V.suurbufr;
	for (k=0;k<zdf;k++)
	{
		for (a=0;a<9;a++)
			Fb(V.suf,i,j,k,DL[a],xdf,ydf,zdf,19)=(*(ptr++));
	}
	
	i=0;j=0;
	for (k=0;k<zdf;k++)
	{
		for (a=0;a<9;a++)
			tpsend[k][a] = Fb(V.slf,i,j,k,UR[a],xdf,ydf,zdf,19);
	}
	MPI_Sendrecv(tpsend,9*zdf,MPI_DOUBLE,pd.left,221,tprecv,9*zdf,MPI_DOUBLE,pd.right,221,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(tpsend,9*zdf,MPI_DOUBLE,pd.dw,222,tprecv2,9*zdf,MPI_DOUBLE,pd.up,222,MPI_COMM_WORLD,&tstatus);
	i=xdf-1;j=0;
	for (k=0;k<zdf;k++)
	{
		for (a=0;a<9;a++)
			Fb(V.slf,i,j,k,UR[a],xdf,ydf,zdf,19) = tprecv[k][a];
	}
	i=0;j=ydf-1;
	for (k=0;k<zdf;k++)
	{
		for (a=0;a<9;a++)
			Fb(V.slf,i,j,k,UR[a],xdf,ydf,zdf,19) = tprecv2[k][a];
	}
	
	i=0;j=0;
	for (k=0;k<zdf;k++)
	{
		for (a=0;a<9;a++)
			tpsend[k][a] = Fb(V.suf,i,j,k,UR[a],xdf,ydf,zdf,19);
	}
	MPI_Sendrecv(tpsend,9*zdf,MPI_DOUBLE,pd.left,221,tprecv,9*zdf,MPI_DOUBLE,pd.right,221,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(tpsend,9*zdf,MPI_DOUBLE,pd.dw,222,tprecv2,9*zdf,MPI_DOUBLE,pd.up,222,MPI_COMM_WORLD,&tstatus);
	i=xdf-1;j=0;
	for (k=0;k<zdf;k++)
	{
		for (a=0;a<9;a++)
			Fb(V.suf,i,j,k,UR[a],xdf,ydf,zdf,19) = tprecv[k][a];
	}
	i=0;j=ydf-1;
	for (k=0;k<zdf;k++)
	{
		for (a=0;a<9;a++)
			Fb(V.suf,i,j,k,UR[a],xdf,ydf,zdf,19) = tprecv2[k][a];
	}
	
	i=0;j=ydf-1;
	for (k=0;k<zdf;k++)
	{
		for (a=0;a<9;a++)
			tpsend[k][a] = Fb(V.slf,i,j,k,UL[a],xdf,ydf,zdf,19);
	}
	MPI_Sendrecv(tpsend,9*zdf,MPI_DOUBLE,pd.ul,223,tprecv,9*zdf,MPI_DOUBLE,pd.dr,223,MPI_COMM_WORLD,&tstatus);
	i=xdf-1;j=0;
	for (k=0;k<zdf;k++)
	{
		for (a=0;a<9;a++)
			Fb(V.slf,i,j,k,UL[a],xdf,ydf,zdf,19) = tprecv[k][a];
	}
	i=0;j=ydf-1;
	for (k=0;k<zdf;k++)
	{
		for (a=0;a<9;a++)
			tpsend[k][a] = Fb(V.suf,i,j,k,UL[a],xdf,ydf,zdf,19);
	}
	MPI_Sendrecv(tpsend,9*zdf,MPI_DOUBLE,pd.ul,223,tprecv,9*zdf,MPI_DOUBLE,pd.dr,223,MPI_COMM_WORLD,&tstatus);
	i=xdf-1;j=0;
	for (k=0;k<zdf;k++)
	{
		for (a=0;a<9;a++)
			Fb(V.suf,i,j,k,UL[a],xdf,ydf,zdf,19) = tprecv[k][a];
	}
	
	i=xdf-1;j=0;
	for (k=0;k<zdf;k++)
	{
		for (a=0;a<9;a++)
			tpsend[k][a] = Fb(V.slf,i,j,k,DR[a],xdf,ydf,zdf,19);
	}
	MPI_Sendrecv(tpsend,9*zdf,MPI_DOUBLE,pd.dr,224,tprecv,9*zdf,MPI_DOUBLE,pd.ul,224,MPI_COMM_WORLD,&tstatus);
	i=0;j=ydf-1;
	for (k=0;k<zdf;k++)
	{
		for (a=0;a<9;a++)
			Fb(V.slf,i,j,k,DR[a],xdf,ydf,zdf,19) = tprecv[k][a];
	}
	i=xdf-1;j=0;
	for (k=0;k<zdf;k++)
	{
		for (a=0;a<9;a++)
			tpsend[k][a] = Fb(V.suf,i,j,k,DR[a],xdf,ydf,zdf,19);
	}
	MPI_Sendrecv(tpsend,9*zdf,MPI_DOUBLE,pd.dr,224,tprecv,9*zdf,MPI_DOUBLE,pd.ul,224,MPI_COMM_WORLD,&tstatus);
	i=0;j=ydf-1;
	for (k=0;k<zdf;k++)
	{
		for (a=0;a<9;a++)
			Fb(V.suf,i,j,k,DR[a],xdf,ydf,zdf,19) = tprecv[k][a];
	}
	
//	fprintf(stderr,"in function xfcomputations, c error handler is %s\n",strerror(errno));
	
	return V;
}
