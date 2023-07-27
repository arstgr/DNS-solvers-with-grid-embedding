/* ****************************************************************************
 * by Amirreza Rastegari                                                      *
 * arstgri@gmail.com                                                          *
 *                                                                            *
 * To be used with postproc-spectra.c                                         *
 *                                                                            *
 * last time tested: 11/20/2014                                               *
 *                                                                            *
 * Computes 1-dimensional energy spectra from DNS of turbulent channel flow   *
 *                                                                            *
 * Inputs:                                                                    *
 *    xd: dimension of velocity array in streamwise direction                 *
 *    yd: dimension of velocity array in spanwise direction                   *
 *    zd: dimension of velocity array in wall-normal direction                *
 *    TSTR: starting time for the calculations                                *
 *    TEND: final time for the calculations                                   *
 *                                                                            *
 * Input files:                                                               *
 *    vel.time: velocity file with array indices running from                 *
 *    (u_x,u_y,u_z), then z index, then y index, then x index                 *
 *    example: vel.0001                                                       *
 *                                                                            *
 * Output files:                                                              *
 *    Euux.time, Evvx.time, Ewwx.time                                         *
 *    Euuy.time, Evvy.time, Ewwy.time                                         *
 *                                                                            *
 * This program uses the fast fourier transform, you need to load the         *
 * necessary modules before the compilation.                                  *
 * To compile on stampede:                                                    *
 * -> module load fftw3                                                       *
 * -> mpicc 1d-spectra.c -I$TACC_FFTW3_INC -L$TACC_FFTW3_LIB -lfftw3_mpi      *
 *    -lfftw3                                                                 *
 ******************************************************************************/
#include "definitions.h"
#include "fun_defs.h"

int main(int argc, char *argv[])
{
	int i,j,k,a;
	int nproc,rank;
	PDATA pd;
	POINTER V;
	int tst;
	
	char fn[80];
	FILE *fd;
	MPI_File fh;
	MPI_Offset offset;
	MPI_Status status;
	int xdf,ydf,zdf,xdc,ydc,zdc,GR;
	
	ptrdiff_t N0, N1;
	fftw_plan plan,inv_plan,plan_u,plan_v,plan_w;
	fftw_complex *data;
	double *rin_u,*rin_v,*rin_w;
	fftw_complex *cout_u,*cout_v,*cout_w;
	ptrdiff_t alloc_local, local_n0, local_n0_start, local_n1_aft, local_n1_start_aft;
	
	double euux,euuy,evvx,evvy,ewwx,ewwy;
	double *Euux,*Euuy,*Evvx,*Evvy,*Ewwx,*Ewwy;
	double *tempy;
	double fac,fac2;
	
	MPI_Init(&argc, &argv);
	fftw_mpi_init();
	
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	pd.myrank=rank;
	pd.numproc=nproc;
	
	if (pd.myrank != 0)
		pd.left = pd.myrank - 1;
	else
		pd.left = pd.numproc - 1;
	
	if (pd.myrank != (pd.numproc-1))
		pd.right = pd.myrank + 1;
	else
		pd.right = 0;

	GR=gratio;
	
	zdc=zdv;
	ydc=ydv;
	xdc=1+(xdv/pd.numproc);
	
	xdf = GR*(xdc-1)+1;
	ydf = GR*ydc;
	zdf = zds+1;	
	
	/***************************************************************************************************/
	
	N0 = xdv;
	N1 = ydv;
	
	local_n0 = N0/nproc;
	local_n0_start = rank*local_n0;

	local_n1_aft = (N1/2 +1)/nproc;
	local_n1_start_aft = rank*local_n1_aft;
	
	/* Allocating the velocity field */
	V.velc = (DP *)calloc((xdc-1)*ydc*zdc*3,sizeof(DP));
	
	/* Allocating the E variables */
	Euux = (DP *)calloc(local_n0*zdc,sizeof(DP));
	Evvx = (DP *)calloc(local_n0*zdc,sizeof(DP));
	Ewwx = (DP *)calloc(local_n0*zdc,sizeof(DP));
	
	Euuy = (DP *)calloc((N1/2 +1)*zdc,sizeof(DP));
	Evvy = (DP *)calloc((N1/2 +1)*zdc,sizeof(DP));
	Ewwy = (DP *)calloc((N1/2 +1)*zdc,sizeof(DP));
	
	tempy = (DP *)calloc((N1/2 +1)*zdc,sizeof(DP));
  
	/* get local data size and allocate */
	alloc_local = fftw_mpi_local_size_2d(N0, (N1/2 +1), MPI_COMM_WORLD, &local_n0, &local_n0_start);
	
	rin_u = fftw_alloc_real(2 * alloc_local);
	rin_v = fftw_alloc_real(2 * alloc_local);
	rin_w = fftw_alloc_real(2 * alloc_local);
	
	cout_u = fftw_alloc_complex(alloc_local);
	cout_v = fftw_alloc_complex(alloc_local);
	cout_w = fftw_alloc_complex(alloc_local);
     
        /* create plan for in-place forward DFT */
	plan_u = fftw_mpi_plan_dft_r2c_2d(N0, N1, rin_u, cout_u, MPI_COMM_WORLD, FFTW_ESTIMATE);
	plan_v = fftw_mpi_plan_dft_r2c_2d(N0, N1, rin_v, cout_v, MPI_COMM_WORLD, FFTW_ESTIMATE);
	plan_w = fftw_mpi_plan_dft_r2c_2d(N0, N1, rin_w, cout_w, MPI_COMM_WORLD, FFTW_ESTIMATE);
//	inv_plan = fftw_mpi_plan_dft_c2r_2d(N0, N1, cout,rin, MPI_COMM_WORLD, FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_IN);
	
	/* main loop */
	for (tst=TSTR;tst<=TEND;tst+=DT)
	{
		V = readingc(xdc,ydc,zdc,V,pd,tst);
		
		for (k=0;k<zdc;k++)
		{
			for (i=0;i<local_n0;i++)
			{
				Euux[i*zdc+k] = 0.;
				Evvx[i*zdc+k] = 0.;
				Ewwx[i*zdc+k] = 0.;
			}
			for (j=0;j<(N1/2+1);j++)
			{
				Euuy[j*zdc+k] = 0.;
				Evvy[j*zdc+k] = 0.;
				Ewwy[j*zdc+k] = 0.;
			}
		}
		
		/* moving the data into the fft plane */
		for (k=0;k<zdc;k++)
		{
			for (i=0;i<local_n0;i++)
			{
				for (j=0;j<N1;j++)
				{
					rin_u[i*(2*(N1/2 +1))+j] = Fb(V.velc,i,j,k,0,xdc,ydc,zdc,3);
					rin_v[i*(2*(N1/2 +1))+j] = Fb(V.velc,i,j,k,1,xdc,ydc,zdc,3);
					rin_w[i*(2*(N1/2 +1))+j] = Fb(V.velc,i,j,k,2,xdc,ydc,zdc,3);
				}
			}
			/* compute transforms, in-place, as many times as desired */
			fftw_execute(plan_u);
			fftw_execute(plan_v);
			fftw_execute(plan_w);
			
			/* Scale the R => F FFT */
			for (i=0; i<local_n0; i++)
			{
				for (j=0;j<(N1/2 +1); j++)
				{
					cout_u[i*(N1/2 +1) +j][0] *= 2.*(1./((double)(2.*N0*N1)));
					cout_u[i*(N1/2 +1) +j][1] *= 2.*(1./((double)(2.*N0*N1)));
					
					cout_v[i*(N1/2 +1) +j][0] *= 2.*(1./((double)(2.*N0*N1)));
					cout_v[i*(N1/2 +1) +j][1] *= 2.*(1./((double)(2.*N0*N1)));
					
					cout_w[i*(N1/2 +1) +j][0] *= 2.*(1./((double)(2.*N0*N1)));
					cout_w[i*(N1/2 +1) +j][1] *= 2.*(1./((double)(2.*N0*N1)));
				}
			}
			
			/* subtract the average */
			if (pd.myrank == 0)
			{
				cout_u[0][0] = 0.; cout_u[0][1] = 0.;
				cout_v[0][0] = 0.; cout_v[0][1] = 0.;
				cout_w[0][0] = 0.; cout_w[0][1] = 0.;
			}
			
			/* Calculation of the energy spectra */
			/* Y spectra */
			for (j=0;j<(N1/2 + 1 );j++)
			{				
				for (i=0; i< local_n0 ;i++)
				{
					if (j==0)
						fac = 1.;
					else
						fac = 2.;

					Euux[i*zdc+k] += fac*0.5*(cout_u[i*(N1/2 +1) +j][0]*cout_u[i*(N1/2 +1) +j][0] + cout_u[i*(N1/2 +1) +j][1]*cout_u[i*(N1/2 +1) +j][1]);
					Evvx[i*zdc+k] += fac*0.5*(cout_v[i*(N1/2 +1) +j][0]*cout_v[i*(N1/2 +1) +j][0] + cout_v[i*(N1/2 +1) +j][1]*cout_v[i*(N1/2 +1) +j][1]);
					Ewwx[i*zdc+k] += fac*0.5*(cout_w[i*(N1/2 +1) +j][0]*cout_w[i*(N1/2 +1) +j][0] + cout_w[i*(N1/2 +1) +j][1]*cout_w[i*(N1/2 +1) +j][1]);
					
					Euuy[j*zdc+k] += 0.5*fac*(cout_u[i*(N1/2 +1) +j][0]*cout_u[i*(N1/2 +1) +j][0] + cout_u[i*(N1/2 +1) +j][1]*cout_u[i*(N1/2 +1) +j][1]);
					Evvy[j*zdc+k] += 0.5*fac*(cout_v[i*(N1/2 +1) +j][0]*cout_v[i*(N1/2 +1) +j][0] + cout_v[i*(N1/2 +1) +j][1]*cout_v[i*(N1/2 +1) +j][1]);
					Ewwy[j*zdc+k] += 0.5*fac*(cout_w[i*(N1/2 +1) +j][0]*cout_w[i*(N1/2 +1) +j][0] + cout_w[i*(N1/2 +1) +j][1]*cout_w[i*(N1/2 +1) +j][1]);
				}
			  
			}
		}
			
		/* Mixing the y direction of all domains */
		MPI_Reduce(Euuy,tempy,(zdc*(N1/2+1)),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		memcpy(Euuy,tempy,(zdc*(N1/2+1))*sizeof(double));
		MPI_Reduce(Evvy,tempy,(zdc*(N1/2+1)),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		memcpy(Evvy,tempy,(zdc*(N1/2+1))*sizeof(double));
		MPI_Reduce(Ewwy,tempy,(zdc*(N1/2+1)),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		memcpy(Ewwy,tempy,(zdc*(N1/2+1))*sizeof(double));
		
		/* Exporting the y results: only by rank 0 */
		if (pd.myrank == 0)
		{
			sprintf(fn,"Euuy-c.%.4d",tst);
			fd=fopen(fn,"wb");
			fwrite(Euuy,sizeof(double),zdc*(N1/2+1),fd);
			fclose(fd);
			
			sprintf(fn,"Evvy-c.%.4d",tst);
			fd=fopen(fn,"wb");
			fwrite(Evvy,sizeof(double),zdc*(N1/2+1),fd);
			fclose(fd);
			
			sprintf(fn,"Ewwy-c.%.4d",tst);
			fd=fopen(fn,"wb");
			fwrite(Ewwy,sizeof(double),zdc*(N1/2+1),fd);
			fclose(fd);
		}
		
		MPI_Barrier(MPI_COMM_WORLD);
			
		/* Exporting the x results, by all in parallel */
		sprintf(fn,"Euux-c.%.4d",tst);
		MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
		MPI_File_set_view(fh,0,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
		MPI_File_write_at_all(fh,(pd.myrank*local_n0*zdc),Euux,(local_n0*zdc),MPI_DOUBLE,&status);
		MPI_File_close(&fh);
		
		sprintf(fn,"Evvx-c.%.4d",tst);
		MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
		MPI_File_set_view(fh,0,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
		MPI_File_write_at_all(fh,(pd.myrank*local_n0*zdc),Evvx,(local_n0*zdc),MPI_DOUBLE,&status);
		MPI_File_close(&fh);
		
		sprintf(fn,"Ewwx-c.%.4d",tst);
		MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
		MPI_File_set_view(fh,0,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
		MPI_File_write_at_all(fh,(pd.myrank*local_n0*zdc),Ewwx,(local_n0*zdc),MPI_DOUBLE,&status);
		MPI_File_close(&fh);
		
	}
     
	/* Destroying the plans */
	fftw_destroy_plan(plan_u);
	fftw_destroy_plan(plan_v);
	fftw_destroy_plan(plan_w);
//	fftw_destroy_plan(inv_plan);
	
	free(Euux);
	free(Euuy);
	free(Evvx);
	free(Evvy);
	free(Ewwx);
	free(Ewwy);
	free(tempy);
	
	/**************************************************************************************************************/
	
	N0 = xdv*GR;
	N1 = ydv*GR;
	
	local_n0 = N0/nproc;
	local_n0_start = rank*local_n0;

	local_n1_aft = (N1/2 +1)/nproc;
	local_n1_start_aft = rank*local_n1_aft;
	
	/* Allocating the velocity field */
	V.velfl = (DP *)calloc((xdf-1)*ydf*zdf*3,sizeof(DP));
	
	/* Allocating the E variables */
	Euux = (DP *)calloc(local_n0*zdf,sizeof(DP));
	Evvx = (DP *)calloc(local_n0*zdf,sizeof(DP));
	Ewwx = (DP *)calloc(local_n0*zdf,sizeof(DP));
	
	Euuy = (DP *)calloc((N1/2 +1)*zdf,sizeof(DP));
	Evvy = (DP *)calloc((N1/2 +1)*zdf,sizeof(DP));
	Ewwy = (DP *)calloc((N1/2 +1)*zdf,sizeof(DP));
	
	tempy = (DP *)calloc((N1/2 +1)*zdf,sizeof(DP));
  
	/* get local data size and allocate */
	alloc_local = fftw_mpi_local_size_2d(N0, (N1/2 +1), MPI_COMM_WORLD, &local_n0, &local_n0_start);
	
	rin_u = fftw_alloc_real(2 * alloc_local);
	rin_v = fftw_alloc_real(2 * alloc_local);
	rin_w = fftw_alloc_real(2 * alloc_local);
	
	cout_u = fftw_alloc_complex(alloc_local);
	cout_v = fftw_alloc_complex(alloc_local);
	cout_w = fftw_alloc_complex(alloc_local);
     
        /* create plan for in-place forward DFT */
	plan_u = fftw_mpi_plan_dft_r2c_2d(N0, N1, rin_u, cout_u, MPI_COMM_WORLD, FFTW_ESTIMATE);
	plan_v = fftw_mpi_plan_dft_r2c_2d(N0, N1, rin_v, cout_v, MPI_COMM_WORLD, FFTW_ESTIMATE);
	plan_w = fftw_mpi_plan_dft_r2c_2d(N0, N1, rin_w, cout_w, MPI_COMM_WORLD, FFTW_ESTIMATE);
//	inv_plan = fftw_mpi_plan_dft_c2r_2d(N0, N1, cout,rin, MPI_COMM_WORLD, FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_IN);
	
	/* main loop */
	for (tst=TSTR;tst<=TEND;tst+=DT)
	{
		V = readingl(xdf,ydf,zdf,V,pd,tst);
		
		for (k=0;k<zdf;k++)
		{
			for (i=0;i<local_n0;i++)
			{
				Euux[i*zdf+k] = 0.;
				Evvx[i*zdf+k] = 0.;
				Ewwx[i*zdf+k] = 0.;
			}
			for (j=0;j<(N1/2+1);j++)
			{
				Euuy[j*zdf+k] = 0.;
				Evvy[j*zdf+k] = 0.;
				Ewwy[j*zdf+k] = 0.;
			}
		}
		
		/* moving the data into the fft plane */
		for (k=0;k<zdf;k++)
		{
			for (i=0;i<local_n0;i++)
			{
				for (j=0;j<N1;j++)
				{
					rin_u[i*(2*(N1/2 +1))+j] = Fb(V.velfl,i,j,k,0,xdf,ydf,zdf,3);
					rin_v[i*(2*(N1/2 +1))+j] = Fb(V.velfl,i,j,k,1,xdf,ydf,zdf,3);
					rin_w[i*(2*(N1/2 +1))+j] = Fb(V.velfl,i,j,k,2,xdf,ydf,zdf,3);
				}
			}
			/* compute transforms, in-place, as many times as desired */
			fftw_execute(plan_u);
			fftw_execute(plan_v);
			fftw_execute(plan_w);
			
			/* Scale the R => F FFT */
			for (i=0; i<local_n0; i++)
			{
				for (j=0;j<(N1/2 +1); j++)
				{
					cout_u[i*(N1/2 +1) +j][0] *= 2.*(1./((double)(2.*N0*N1)));
					cout_u[i*(N1/2 +1) +j][1] *= 2.*(1./((double)(2.*N0*N1)));
					
					cout_v[i*(N1/2 +1) +j][0] *= 2.*(1./((double)(2.*N0*N1)));
					cout_v[i*(N1/2 +1) +j][1] *= 2.*(1./((double)(2.*N0*N1)));
					
					cout_w[i*(N1/2 +1) +j][0] *= 2.*(1./((double)(2.*N0*N1)));
					cout_w[i*(N1/2 +1) +j][1] *= 2.*(1./((double)(2.*N0*N1)));
				}
			}
			
			/* subtract the average */
			if (pd.myrank == 0)
			{
				cout_u[0][0] = 0.; cout_u[0][1] = 0.;
				cout_v[0][0] = 0.; cout_v[0][1] = 0.;
				cout_w[0][0] = 0.; cout_w[0][1] = 0.;
			}
			
			/* Calculation of the energy spectra */
			/* Y spectra */
			for (j=0;j<(N1/2 + 1 );j++)
			{				
				for (i=0; i< local_n0 ;i++)
				{
					if (j==0)
						fac = 1.;
					else
						fac = 2.;

					Euux[i*zdf+k] += fac*0.5*(cout_u[i*(N1/2 +1) +j][0]*cout_u[i*(N1/2 +1) +j][0] + cout_u[i*(N1/2 +1) +j][1]*cout_u[i*(N1/2 +1) +j][1]);
					Evvx[i*zdf+k] += fac*0.5*(cout_v[i*(N1/2 +1) +j][0]*cout_v[i*(N1/2 +1) +j][0] + cout_v[i*(N1/2 +1) +j][1]*cout_v[i*(N1/2 +1) +j][1]);
					Ewwx[i*zdf+k] += fac*0.5*(cout_w[i*(N1/2 +1) +j][0]*cout_w[i*(N1/2 +1) +j][0] + cout_w[i*(N1/2 +1) +j][1]*cout_w[i*(N1/2 +1) +j][1]);
					
					Euuy[j*zdf+k] += 0.5*fac*(cout_u[i*(N1/2 +1) +j][0]*cout_u[i*(N1/2 +1) +j][0] + cout_u[i*(N1/2 +1) +j][1]*cout_u[i*(N1/2 +1) +j][1]);
					Evvy[j*zdf+k] += 0.5*fac*(cout_v[i*(N1/2 +1) +j][0]*cout_v[i*(N1/2 +1) +j][0] + cout_v[i*(N1/2 +1) +j][1]*cout_v[i*(N1/2 +1) +j][1]);
					Ewwy[j*zdf+k] += 0.5*fac*(cout_w[i*(N1/2 +1) +j][0]*cout_w[i*(N1/2 +1) +j][0] + cout_w[i*(N1/2 +1) +j][1]*cout_w[i*(N1/2 +1) +j][1]);
				}
			  
			}
		}
			
		/* Mixing the y direction of all domains */
		MPI_Reduce(Euuy,tempy,(zdf*(N1/2+1)),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		memcpy(Euuy,tempy,(zdf*(N1/2+1))*sizeof(double));
		MPI_Reduce(Evvy,tempy,(zdf*(N1/2+1)),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		memcpy(Evvy,tempy,(zdf*(N1/2+1))*sizeof(double));
		MPI_Reduce(Ewwy,tempy,(zdf*(N1/2+1)),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		memcpy(Ewwy,tempy,(zdf*(N1/2+1))*sizeof(double));
		
		/* Exporting the y results: only by rank 0 */
		if (pd.myrank == 0)
		{
			sprintf(fn,"Euuy-fl.%.4d",tst);
			fd=fopen(fn,"wb");
			fwrite(Euuy,sizeof(double),zdf*(N1/2+1),fd);
			fclose(fd);
			
			sprintf(fn,"Evvy-fl.%.4d",tst);
			fd=fopen(fn,"wb");
			fwrite(Evvy,sizeof(double),zdf*(N1/2+1),fd);
			fclose(fd);
			
			sprintf(fn,"Ewwy-fl.%.4d",tst);
			fd=fopen(fn,"wb");
			fwrite(Ewwy,sizeof(double),zdf*(N1/2+1),fd);
			fclose(fd);
		}
		
		MPI_Barrier(MPI_COMM_WORLD);
			
		/* Exporting the x results, by all in parallel */
		sprintf(fn,"Euux-fl.%.4d",tst);
		MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
		MPI_File_set_view(fh,0,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
		MPI_File_write_at_all(fh,(pd.myrank*local_n0*zdf),Euux,(local_n0*zdf),MPI_DOUBLE,&status);
		MPI_File_close(&fh);
		
		sprintf(fn,"Evvx-fl.%.4d",tst);
		MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
		MPI_File_set_view(fh,0,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
		MPI_File_write_at_all(fh,(pd.myrank*local_n0*zdf),Evvx,(local_n0*zdf),MPI_DOUBLE,&status);
		MPI_File_close(&fh);
		
		sprintf(fn,"Ewwx-fl.%.4d",tst);
		MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
		MPI_File_set_view(fh,0,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
		MPI_File_write_at_all(fh,(pd.myrank*local_n0*zdf),Ewwx,(local_n0*zdf),MPI_DOUBLE,&status);
		MPI_File_close(&fh);
	}
     
	/* Destroying the plans */
	fftw_destroy_plan(plan_u);
	fftw_destroy_plan(plan_v);
	fftw_destroy_plan(plan_w);
//	fftw_destroy_plan(inv_plan);
	
	free(Euux);
	free(Euuy);
	free(Evvx);
	free(Evvy);
	free(Ewwx);
	free(Ewwy);
	free(tempy);
	
/**************************************************************************************************************/
	
	N0 = xdv*GR;
	N1 = ydv*GR;
	
	local_n0 = N0/nproc;
	local_n0_start = rank*local_n0;

	local_n1_aft = (N1/2 +1)/nproc;
	local_n1_start_aft = rank*local_n1_aft;
	
	/* Allocating the velocity field */
	V.velfu = (DP *)calloc((xdf-1)*ydf*zdf*3,sizeof(DP));
	
	/* Allocating the E variables */
	Euux = (DP *)calloc(local_n0*zdf,sizeof(DP));
	Evvx = (DP *)calloc(local_n0*zdf,sizeof(DP));
	Ewwx = (DP *)calloc(local_n0*zdf,sizeof(DP));
	
	Euuy = (DP *)calloc((N1/2 +1)*zdf,sizeof(DP));
	Evvy = (DP *)calloc((N1/2 +1)*zdf,sizeof(DP));
	Ewwy = (DP *)calloc((N1/2 +1)*zdf,sizeof(DP));
	
	tempy = (DP *)calloc((N1/2 +1)*zdf,sizeof(DP));
  
	/* get local data size and allocate */
	alloc_local = fftw_mpi_local_size_2d(N0, (N1/2 +1), MPI_COMM_WORLD, &local_n0, &local_n0_start);
	
	rin_u = fftw_alloc_real(2 * alloc_local);
	rin_v = fftw_alloc_real(2 * alloc_local);
	rin_w = fftw_alloc_real(2 * alloc_local);
	
	cout_u = fftw_alloc_complex(alloc_local);
	cout_v = fftw_alloc_complex(alloc_local);
	cout_w = fftw_alloc_complex(alloc_local);
     
        /* create plan for in-place forward DFT */
	plan_u = fftw_mpi_plan_dft_r2c_2d(N0, N1, rin_u, cout_u, MPI_COMM_WORLD, FFTW_ESTIMATE);
	plan_v = fftw_mpi_plan_dft_r2c_2d(N0, N1, rin_v, cout_v, MPI_COMM_WORLD, FFTW_ESTIMATE);
	plan_w = fftw_mpi_plan_dft_r2c_2d(N0, N1, rin_w, cout_w, MPI_COMM_WORLD, FFTW_ESTIMATE);
//	inv_plan = fftw_mpi_plan_dft_c2r_2d(N0, N1, cout,rin, MPI_COMM_WORLD, FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_IN);
	
	/* main loop */
	for (tst=TSTR;tst<=TEND;tst+=DT)
	{
		V = readingu(xdf,ydf,zdf,V,pd,tst);
		
		for (k=0;k<zdf;k++)
		{
			for (i=0;i<local_n0;i++)
			{
				Euux[i*zdf+k] = 0.;
				Evvx[i*zdf+k] = 0.;
				Ewwx[i*zdf+k] = 0.;
			}
			for (j=0;j<(N1/2+1);j++)
			{
				Euuy[j*zdf+k] = 0.;
				Evvy[j*zdf+k] = 0.;
				Ewwy[j*zdf+k] = 0.;
			}
		}
		
		/* moving the data into the fft plane */
		for (k=0;k<zdf;k++)
		{
			for (i=0;i<local_n0;i++)
			{
				for (j=0;j<N1;j++)
				{
					rin_u[i*(2*(N1/2 +1))+j] = Fb(V.velfu,i,j,k,0,xdf,ydf,zdf,3);
					rin_v[i*(2*(N1/2 +1))+j] = Fb(V.velfu,i,j,k,1,xdf,ydf,zdf,3);
					rin_w[i*(2*(N1/2 +1))+j] = Fb(V.velfu,i,j,k,2,xdf,ydf,zdf,3);
				}
			}
			/* compute transforms, in-place, as many times as desired */
			fftw_execute(plan_u);
			fftw_execute(plan_v);
			fftw_execute(plan_w);
			
			/* Scale the R => F FFT */
			for (i=0; i<local_n0; i++)
			{
				for (j=0;j<(N1/2 +1); j++)
				{
					cout_u[i*(N1/2 +1) +j][0] *= 2.*(1./((double)(2.*N0*N1)));
					cout_u[i*(N1/2 +1) +j][1] *= 2.*(1./((double)(2.*N0*N1)));
					
					cout_v[i*(N1/2 +1) +j][0] *= 2.*(1./((double)(2.*N0*N1)));
					cout_v[i*(N1/2 +1) +j][1] *= 2.*(1./((double)(2.*N0*N1)));
					
					cout_w[i*(N1/2 +1) +j][0] *= 2.*(1./((double)(2.*N0*N1)));
					cout_w[i*(N1/2 +1) +j][1] *= 2.*(1./((double)(2.*N0*N1)));
				}
			}
			
			/* subtract the average */
			if (pd.myrank == 0)
			{
				cout_u[0][0] = 0.; cout_u[0][1] = 0.;
				cout_v[0][0] = 0.; cout_v[0][1] = 0.;
				cout_w[0][0] = 0.; cout_w[0][1] = 0.;
			}
			
			/* Calculation of the energy spectra */
			/* Y spectra */
			for (j=0;j<(N1/2 + 1 );j++)
			{
				for (i=0; i< local_n0 ;i++)
				{
					if (j==0)
						fac = 1.;
					else
						fac = 2.;

					Euux[i*zdf+k] += fac*0.5*(cout_u[i*(N1/2 +1) +j][0]*cout_u[i*(N1/2 +1) +j][0] + cout_u[i*(N1/2 +1) +j][1]*cout_u[i*(N1/2 +1) +j][1]);
					Evvx[i*zdf+k] += fac*0.5*(cout_v[i*(N1/2 +1) +j][0]*cout_v[i*(N1/2 +1) +j][0] + cout_v[i*(N1/2 +1) +j][1]*cout_v[i*(N1/2 +1) +j][1]);
					Ewwx[i*zdf+k] += fac*0.5*(cout_w[i*(N1/2 +1) +j][0]*cout_w[i*(N1/2 +1) +j][0] + cout_w[i*(N1/2 +1) +j][1]*cout_w[i*(N1/2 +1) +j][1]);
					
					Euuy[j*zdf+k] += 0.5*fac*(cout_u[i*(N1/2 +1) +j][0]*cout_u[i*(N1/2 +1) +j][0] + cout_u[i*(N1/2 +1) +j][1]*cout_u[i*(N1/2 +1) +j][1]);
					Evvy[j*zdf+k] += 0.5*fac*(cout_v[i*(N1/2 +1) +j][0]*cout_v[i*(N1/2 +1) +j][0] + cout_v[i*(N1/2 +1) +j][1]*cout_v[i*(N1/2 +1) +j][1]);
					Ewwy[j*zdf+k] += 0.5*fac*(cout_w[i*(N1/2 +1) +j][0]*cout_w[i*(N1/2 +1) +j][0] + cout_w[i*(N1/2 +1) +j][1]*cout_w[i*(N1/2 +1) +j][1]);
				}
			  
			}
		}
			
		/* Mixing the y direction of all domains */
		MPI_Reduce(Euuy,tempy,(zdf*(N1/2+1)),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		memcpy(Euuy,tempy,(zdf*(N1/2+1))*sizeof(double));
		MPI_Reduce(Evvy,tempy,(zdf*(N1/2+1)),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		memcpy(Evvy,tempy,(zdf*(N1/2+1))*sizeof(double));
		MPI_Reduce(Ewwy,tempy,(zdf*(N1/2+1)),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		memcpy(Ewwy,tempy,(zdf*(N1/2+1))*sizeof(double));
		
		/* Exporting the y results: only by rank 0 */
		if (pd.myrank == 0)
		{
			sprintf(fn,"Euuy-fu.%.4d",tst);
			fd=fopen(fn,"wb");
			fwrite(Euuy,sizeof(double),zdf*(N1/2+1),fd);
			fclose(fd);
			
			sprintf(fn,"Evvy-fu.%.4d",tst);
			fd=fopen(fn,"wb");
			fwrite(Evvy,sizeof(double),zdf*(N1/2+1),fd);
			fclose(fd);
			
			sprintf(fn,"Ewwy-fu.%.4d",tst);
			fd=fopen(fn,"wb");
			fwrite(Ewwy,sizeof(double),zdf*(N1/2+1),fd);
			fclose(fd);
		}
		
		MPI_Barrier(MPI_COMM_WORLD);
			
		/* Exporting the x results, by all in parallel */
		sprintf(fn,"Euux-fu.%.4d",tst);
		MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
		MPI_File_set_view(fh,0,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
		MPI_File_write_at_all(fh,(pd.myrank*local_n0*zdf),Euux,(local_n0*zdf),MPI_DOUBLE,&status);
		MPI_File_close(&fh);
		
		sprintf(fn,"Evvx-fu.%.4d",tst);
		MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
		MPI_File_set_view(fh,0,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
		MPI_File_write_at_all(fh,(pd.myrank*local_n0*zdf),Evvx,(local_n0*zdf),MPI_DOUBLE,&status);
		MPI_File_close(&fh);
		
		sprintf(fn,"Ewwx-fu.%.4d",tst);
		MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
		MPI_File_set_view(fh,0,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
		MPI_File_write_at_all(fh,(pd.myrank*local_n0*zdf),Ewwx,(local_n0*zdf),MPI_DOUBLE,&status);
		MPI_File_close(&fh);
		
	}
     
	/* Destroying the plans */
	fftw_destroy_plan(plan_u);
	fftw_destroy_plan(plan_v);
	fftw_destroy_plan(plan_w);
//	fftw_destroy_plan(inv_plan);
	
	free(Euux);
	free(Euuy);
	free(Evvx);
	free(Evvy);
	free(Ewwx);
	free(Ewwy);
	free(tempy);
	
	MPI_Finalize();
  
 return 0; 
}
