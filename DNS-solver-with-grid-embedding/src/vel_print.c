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

POINTER vel_print(int xdf, int ydf, int zdf, int xdc, int ydc, int zdc, POINTER V, long cntr)
{
	int i,j,k,a,q;
	char fl[30];
	double ux,uy,uz,dens,*tmp,*tmpd;
	MPI_File fh,fhd;
	MPI_Offset offset,offsetd;
	MPI_Status status,statusd;
	int NZ=2*zdf+(zdc-4);
	
	int coord[2], rank, size;
	MPI_Datatype block,cblock,fblock,sarray,fine_sarray,d_sarray,fine_d_sarray;
	MPI_Datatype cdblock,fdblock;
	MPI_Info info;
	int gsizes[2],lsizes[2],psizes[2],array_start[2],memsizes[2];

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Cart_coords(cart_grid, rank, 2, coord);
	
	MPI_Info_create(&info);
	MPI_Info_set(info, "striping_factor", STRIPE_COUNT);
//	MPI_Info_set(info, "striping_unit", STRIPE_SIZE);
	MPI_Info_set(info, "romio_cb_write", "enable");
//	MPI_Info_set(info, "romio_ds_write", "disable");
	
	MPI_Type_contiguous(3,MPI_DOUBLE,&block);
	MPI_Type_commit(&block);
	MPI_Type_contiguous(zdc,MPI_DOUBLE,&cdblock);
	MPI_Type_contiguous(zdf,MPI_DOUBLE,&fdblock);
	MPI_Type_contiguous(zdc,block,&cblock);
	MPI_Type_contiguous(zdf,block,&fblock);
	
	MPI_Type_commit(&cdblock);
	MPI_Type_commit(&fdblock);
	MPI_Type_commit(&cblock);
	MPI_Type_commit(&fblock);
	
	gsizes[0] = (xdc-1)*XDIM;    gsizes[1] = (ydc-1)*YDIM;     
	psizes[0] = XDIM;    psizes[1] = YDIM;
	lsizes[0] = xdc-1; lsizes[1] = ydc-1;
	array_start[0] = coord[0] * lsizes[0];
	array_start[1] = coord[1] * lsizes[1];
	MPI_Type_create_subarray(2, gsizes, lsizes, array_start, MPI_ORDER_C, cblock, &sarray);
	MPI_Type_commit(&sarray);
	MPI_Type_create_subarray(2, gsizes, lsizes, array_start, MPI_ORDER_C, cdblock, &d_sarray);
	MPI_Type_commit(&d_sarray);
	
	gsizes[0] = (xdf-1)*XDIM;    gsizes[1] = (ydf-1)*YDIM;     
	psizes[0] = XDIM;    psizes[1] = YDIM;
	lsizes[0] = xdf-1; lsizes[1] = ydf-1;
	array_start[0] = coord[0] * lsizes[0];
	array_start[1] = coord[1] * lsizes[1];
	MPI_Type_create_subarray(2, gsizes, lsizes, array_start, MPI_ORDER_C, fblock, &fine_sarray);
	MPI_Type_commit(&fine_sarray);
	MPI_Type_create_subarray(2, gsizes, lsizes, array_start, MPI_ORDER_C, fdblock, &fine_d_sarray);
	MPI_Type_commit(&fine_d_sarray);

	tmp=(double *)calloc(3*zdc*(ydc-1)*(xdc-1),sizeof(double));
	tmpd=(double *)calloc(zdc*(ydc-1)*(xdc-1),sizeof(double));// it should be xdc-1
	if (tmp==NULL)
		fprintf(stderr,"Unable to allocte memory for vel_print function tmp\n error code: %s\n",strerror(errno));

	if (tmpd==NULL)
		fprintf(stderr,"Unable to allocte memory for vel_print function tmpd\n error code: %s\n",strerror(errno));

	sprintf(fl,"vel-c.%.4d",((int)cntr));
	MPI_File_open(MPI_COMM_WORLD,fl,MPI_MODE_CREATE | MPI_MODE_WRONLY,info,&fh);
	MPI_File_set_view(fh, 0, MPI_DOUBLE, sarray, "native", info);

	sprintf(fl,"den-c.%.4d",((int)cntr));
	MPI_File_open(MPI_COMM_WORLD,fl,MPI_MODE_CREATE | MPI_MODE_WRONLY,info,&fhd);
	MPI_File_set_view(fhd, 0, MPI_DOUBLE, d_sarray, "native", info);


	for (i=0;i<(xdc-1);i++)
	{
		for (j=0;j<(ydc-1);j++)
		{
			for (k=0;k<zdc;k++)
			{
				ux=0.; uy=0.; uz=0.; dens=0.;
				for (a=0;a<19;a++)
				{
					ux += E0[a]*Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);
					uy += E1[a]*Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);
					uz += E2[a]*Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);
					dens += Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);
				}
				if (dens>0.)
				{
					ux /= dens;
					uy /= dens;
					uz /= dens;
				}
				else
				{
					ux =0.;
					uy =0.;
					uz =0.;
					dens=0.;
				}
// 				*(tmp+j*zdc*3+k*3+0)=ux;
// 				*(tmp+j*zdc*3+k*3+1)=uy;
// 				*(tmp+j*zdc*3+k*3+2)=uz;
// 				*(tmpd+j*zdc+k)=dens;

				*(tmp+i*(ydc-1)*zdc*3+j*zdc*3+k*3+0)=ux;
				*(tmp+i*(ydc-1)*zdc*3+j*zdc*3+k*3+1)=uy;
				*(tmp+i*(ydc-1)*zdc*3+j*zdc*3+k*3+2)=uz;
 				*(tmpd+i*(ydc-1)*zdc+j*zdc+k)=dens;
			}
		}
// 		MPI_File_write_all(fh,tmp,(ydc-1)*zdc*3,MPI_DOUBLE,&status);
// 		MPI_File_write_all(fhd,tmpd,(ydc-1)*zdc,MPI_DOUBLE,&statusd);
	}
 	MPI_File_write_all(fh,tmp,(xdc-1)*(ydc-1)*zdc*3,MPI_DOUBLE,&status);
	MPI_File_write_all(fhd,tmpd,(xdc-1)*(ydc-1)*zdc,MPI_DOUBLE,&statusd);
	free(tmp);
	free(tmpd);
	MPI_File_close(&fh);
	MPI_File_close(&fhd);

	tmp=(double *)calloc(3*zdf*(ydf-1)*(xdf-1),sizeof(double));
	tmpd=(double *)calloc(zdf*(ydf-1)*(xdf-1),sizeof(double));
	if (tmp==NULL)
		fprintf(stderr,"Unable to allocte memory for vel_print function tmp\n error code: %s\n",strerror(errno));

	if (tmpd==NULL)
		fprintf(stderr,"Unable to allocte memory for vel_print function tmpd\n error code: %s\n",strerror(errno));

	sprintf(fl,"vel-suf.%.4d",((int)cntr));
	MPI_File_open(MPI_COMM_WORLD,fl,MPI_MODE_CREATE | MPI_MODE_WRONLY,info,&fh);
	MPI_File_set_view(fh,0,MPI_DOUBLE,fine_sarray,"native",info);

	sprintf(fl,"den-suf.%.4d",((int)cntr));
	MPI_File_open(MPI_COMM_WORLD,fl,MPI_MODE_CREATE | MPI_MODE_WRONLY,info,&fhd);
	MPI_File_set_view(fhd,0,MPI_DOUBLE,fine_d_sarray,"native",info);


	for (i=0;i<(xdf-1);i++)
	{
		for (j=0;j<(ydf-1);j++)
		{
			for (k=0;k<zdf;k++)
			{
				ux=0.; uy=0.; uz=0.; dens=0.;
				for (a=0;a<19;a++)
				{
					ux += E0[a]*Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
					uy += E1[a]*Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
					uz += E2[a]*Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
					dens += Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
				}
				if (dens>0.)
				{
					ux /= dens;
					uy /= dens;
					uz /= dens;
				}
				else
				{
					ux =0.;
					uy =0.;
					uz =0.;
					dens=0.;
				}
// 				*(tmp+j*zdf*3+k*3+0)=ux;
// 				*(tmp+j*zdf*3+k*3+1)=uy;
// 				*(tmp+j*zdf*3+k*3+2)=uz;
// 				*(tmpd+j*zdf+k)=dens;
				
				*(tmp+i*(ydf-1)*zdf*3+j*zdf*3+k*3+0)=ux;
				*(tmp+i*(ydf-1)*zdf*3+j*zdf*3+k*3+1)=uy;
				*(tmp+i*(ydf-1)*zdf*3+j*zdf*3+k*3+2)=uz;
				*(tmpd+i*(ydf-1)*zdf+j*zdf+k)=dens;
			}
		}
// 		MPI_File_write_all(fh,tmp,(ydf-1)*zdf*3,MPI_DOUBLE,&status);
// 		MPI_File_write_all(fhd,tmpd,(ydf-1)*zdf,MPI_DOUBLE,&statusd);
	}
	MPI_File_write_all(fh,tmp,(xdf-1)*(ydf-1)*zdf*3,MPI_DOUBLE,&status);
	MPI_File_write_all(fhd,tmpd,(xdf-1)*(ydf-1)*zdf,MPI_DOUBLE,&statusd);
	
	MPI_File_close(&fh);
	MPI_File_close(&fhd);
	
	sprintf(fl,"vel-slf.%.4d",((int)cntr));
	MPI_File_open(MPI_COMM_WORLD,fl,MPI_MODE_CREATE | MPI_MODE_WRONLY,info,&fh);
	MPI_File_set_view(fh,0,MPI_DOUBLE,fine_sarray,"native",info);

	sprintf(fl,"den-slf.%.4d",((int)cntr));
	MPI_File_open(MPI_COMM_WORLD,fl,MPI_MODE_CREATE | MPI_MODE_WRONLY,info,&fhd);
	MPI_File_set_view(fhd,0,MPI_DOUBLE,fine_d_sarray,"native",info);


	for (i=0;i<(xdf-1);i++)
	{
		for (j=0;j<(ydf-1);j++)
		{
			for (k=0;k<zdf;k++)
			{
				ux=0.; uy=0.; uz=0.; dens=0.;
				for (a=0;a<19;a++)
				{
					ux += E0[a]*Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
					uy += E1[a]*Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
					uz += E2[a]*Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
					dens += Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);
				}
				if (dens>0.)
				{
					ux /= dens;
					uy /= dens;
					uz /= dens;
				}
				else
				{
					ux =0.;
					uy =0.;
					uz =0.;
					dens=0.;
				}
// 				*(tmp+j*zdf*3+k*3+0)=ux;
// 				*(tmp+j*zdf*3+k*3+1)=uy;
// 				*(tmp+j*zdf*3+k*3+2)=uz;
// 				*(tmpd+j*zdf+k)=dens;

				*(tmp+i*(ydf-1)*zdf*3+j*zdf*3+k*3+0)=ux;
				*(tmp+i*(ydf-1)*zdf*3+j*zdf*3+k*3+1)=uy;
				*(tmp+i*(ydf-1)*zdf*3+j*zdf*3+k*3+2)=uz;
				*(tmpd+i*(ydf-1)*zdf+j*zdf+k)=dens;
			}
		}
// 		MPI_File_write_all(fh,tmp,(ydf-1)*zdf*3,MPI_DOUBLE,&status);
// 		MPI_File_write_all(fhd,tmpd,(ydf-1)*zdf,MPI_DOUBLE,&statusd);
	}
	MPI_File_write_all(fh,tmp,(xdf-1)*(ydf-1)*zdf*3,MPI_DOUBLE,&status);
	MPI_File_write_all(fhd,tmpd,(xdf-1)*(ydf-1)*zdf,MPI_DOUBLE,&statusd);
	
	MPI_File_close(&fh);
	MPI_File_close(&fhd);
	free(tmp);
	free(tmpd);
	
	MPI_Type_free(&block);
	MPI_Type_free(&cdblock);
	MPI_Type_free(&fdblock);
	MPI_Type_free(&cblock);
	MPI_Type_free(&fblock);
	
	MPI_Type_free(&sarray);
	MPI_Type_free(&d_sarray);
	
	MPI_Type_free(&fine_sarray);
	MPI_Type_free(&fine_d_sarray);
	
	MPI_Info_free(&info);
		
	return V;
}