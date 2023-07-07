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
#include "fun_defs.h"

int main(int argc, char *argv[])
{
	
	int xd,yd,zd;
	int NZs;
	int TM;
	int xdc,ydc,zdc;
	int xdf,ydf,zdf;
	int GR;
	int i,j,k,q,a,b,c,nut;
	int tagl=1,tagr=2;
	int cntr,code,iter;
	long ss,up;
	PDATA pd;
	POINTER V;
	DP coef;
	DP Fx=0.,Fy=0.,Fz=0.;
	DP Gx,Gy,Gz;
	DP tau=0.,tauc=0.,tauf=0.;
	DP rhozero=1.,rho;
	DP t1,t2,t3,t4,t5;
	FILE *tm,*sv;
	char fn[60],ch;
	DP utime,umax,dt;
	DP ux,uy,uz,dens;
	DP Ret, Reb, ut,dplus,Lplus,ubulk;

	/* MPI variables */
	MPI_File fh;
	MPI_Info info;
	MPI_Status status,tstatus;
	MPI_Offset offset;
	int dim[2], period[2], reorder;
	int coord[2], neigh[2], id;
	int rank, size, error;
	MPI_Datatype distf,block,sarray,larray,send_to_dw,recv_from_up,send_to_le,recv_from_ri;
	MPI_Datatype sfblock,sfsarray,sflarray,sfsend_to_dw,sfrecv_from_up,sfsend_to_le,sfrecv_from_ri;
	MPI_Datatype f_send_to_dl,f_recv_from_ur,f_sfsend_to_dl,f_sfrecv_from_ur;
	MPI_Datatype intplelmnt,intpl_send_to_dw,intpl_send_to_up,intpl_recv_fr_up,intpl_recv_fr_dw,intpllarray,intplsarray;
	MPI_Datatype intpl_send_to_ri,intpl_send_to_le,intpl_recv_fr_le,intpl_recv_fr_ri;
	int gsizes[2],lsizes[2],psizes[2],array_start[2],memsizes[2];
	MPI_Datatype intpl_send_to_ul,intpl_send_to_dr,intpl_send_to_ur,intpl_send_to_dl;
	MPI_Datatype intpl_recv_fr_ul,intpl_recv_fr_dr,intpl_recv_fr_ur,intpl_recv_fr_dl;
	char pcs_name[MPI_MAX_PROCESSOR_NAME];
	int pcs_name_length;
	
	/* MPI Initialization */
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Get_processor_name(pcs_name, &pcs_name_length);
	printf("%s: is running the process %d of total %d processors\n", pcs_name, rank, size);
	MPI_Barrier(MPI_COMM_WORLD);
	mpi_errors = MPI_Comm_set_errhandler(MPI_COMM_WORLD,MPI_ERRORS_RETURN);
	
	MPI_Info_create(&info);
	MPI_Info_set(info, "striping_factor", STRIPE_COUNT);
//	MPI_Info_set(info, "striping_unit", STRIPE_SIZE);
	MPI_Info_set(info, "romio_cb_write", "enable");
//	MPI_Info_set(info, "romio_ds_write", "disable");
	
	pd.myrank = rank;
	pd.numproc = size;
	
	if (XDIM * YDIM != pd.numproc)
	{
		printf("Incorrect dimensions for the topology\n");
		fflush(stdout);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	dim[0]=XDIM; dim[1]=YDIM;
	period[0]=1; period[1]=1;
	reorder=1;
	MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &cart_grid);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	pd.myrank = rank;
	MPI_Cart_coords(cart_grid, pd.myrank, 2, coord);
	
	neigh[0] = coord[0];
	neigh[1] = coord[1] + 1;
	MPI_Cart_rank(cart_grid, neigh, &rank);
	pd.up = rank;
	
	neigh[0] = coord[0];
	neigh[1] = coord[1] - 1;
	MPI_Cart_rank(cart_grid, neigh, &rank);
	pd.dw = rank;
	
	neigh[0] = coord[0] + 1;
	neigh[1] = coord[1];
	MPI_Cart_rank(cart_grid, neigh, &rank);
	pd.right = rank;
	
	neigh[0] = coord[0] - 1;
	neigh[1] = coord[1];
	MPI_Cart_rank(cart_grid, neigh, &rank);
	pd.left = rank;
	
	neigh[0] = coord[0] + 1;
	neigh[1] = coord[1] + 1;
	MPI_Cart_rank(cart_grid, neigh, &rank);
	pd.ur = rank;
	
	neigh[0] = coord[0] - 1;
	neigh[1] = coord[1] + 1;
	MPI_Cart_rank(cart_grid, neigh, &rank);
	pd.ul = rank;
	
	neigh[0] = coord[0] - 1;
	neigh[1] = coord[1] - 1;
	MPI_Cart_rank(cart_grid, neigh, &rank);
	pd.dl = rank;
	
	neigh[0] = coord[0] + 1;
	neigh[1] = coord[1] - 1;
	MPI_Cart_rank(cart_grid, neigh, &rank);
	pd.dr = rank;
	
	/* Grid Ratio */
	GR=gratio;
	
	/* Pi */
	PI = 4.*atan(1.);

	/* grid dimensions of the uniform grid for initialization */
	zd=zdu+2;
	yd=ydu;
	xd=1+(xdu/pd.numproc);

	/* Coarse grid's dimensions */
	xdc = 1+(xdv/XDIM);
	ydc = 1+(ydv/YDIM);
	zdc = zdv;

	/* Fine grid's dimensions */
	xdf = GR*(xdc-1)+1;
	ydf = GR*(ydc-1)+1;
	zdf = zds+1;	

	/* Re-defining zd: this is the total size of the grid in terms of fine grid's resolution */
	NZs = 2*zdf+GR*(zdc-4)+GR-1; // This number has the additional +2 value
	
	MPI_Type_contiguous(19,MPI_DOUBLE,&distf);
	MPI_Type_commit(&distf);
	MPI_Type_contiguous(23,MPI_DOUBLE,&intplelmnt);
	MPI_Type_contiguous(zdc,distf,&block);
	MPI_Type_contiguous(zdf,distf,&sfblock);
	MPI_Type_commit(&intplelmnt);
	MPI_Type_commit(&block);
	MPI_Type_commit(&sfblock);
	
	/* On coarse grid */
	/* no. of rows and columns in global array*/ 
	gsizes[0] = (xdc-1)*XDIM;    gsizes[1] = (ydc-1)*YDIM;     
	/* no. of processors in x and y directions */ 
	psizes[0] = XDIM;    psizes[1] = YDIM;
	/* size of the local array without the ghost cells */
	lsizes[0] = xdc-1; lsizes[1] = ydc-1;
	/* start of the array in terms of global array indices, ignoring the ghost cells */
	array_start[0] = coord[0] * lsizes[0];
	array_start[1] = coord[1] * lsizes[1];
	
	/* Creating the type for the arrays withou ghost cells
	 * the idea is to read them with a type of array with ghost cells */
	MPI_Type_create_subarray(2, gsizes, lsizes, array_start, MPI_ORDER_C, block, &sarray);
	MPI_Type_commit(&sarray);
	
	/* size of the local array with the ghost cells */
	memsizes[0] = lsizes[0] + 1;
	memsizes[1] = lsizes[1] + 1;
	
	/* start of the local array with ghost cells: note: ghost cells are the last row and columns */
	array_start[0] = 0;
	array_start[1] = 0;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &larray);
	MPI_Type_commit(&larray);
	
	/* data type for transfer in y direction */
	lsizes[0] = xdc-1;
	lsizes[1] = 1;
	array_start[0] = 0;
	array_start[1] = 0;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &send_to_dw);
	MPI_Type_commit(&send_to_dw);
	
	array_start[0] = 0;
	array_start[1] = ydc-1;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &recv_from_up);
	MPI_Type_commit(&recv_from_up);
	
	/* in the x direction */
	lsizes[0] = 1;
	lsizes[1] = ydc-1;
	array_start[0] = 0;
	array_start[1] = 0;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &send_to_le);
	MPI_Type_commit(&send_to_le);
	
	array_start[0] = xdc-1;
	array_start[1] = 0;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &recv_from_ri);
	MPI_Type_commit(&recv_from_ri);
	
	/* corner */
	lsizes[0] = 1;
	lsizes[1] = 1;
	array_start[0] = 0;
	array_start[1] = 0;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &f_send_to_dl);
	MPI_Type_commit(&f_send_to_dl);
	
	array_start[0] = xdc-1;
	array_start[1] = ydc-1;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &f_recv_from_ur);
	MPI_Type_commit(&f_recv_from_ur);
	
	/* On fine grid */
	/* no. of rows and columns in global array*/ 
	gsizes[0] = (xdf-1)*XDIM;    gsizes[1] = (ydf-1)*YDIM;     
	/* no. of processors in x and y directions */ 
	psizes[0] = XDIM;    psizes[1] = YDIM;
	/* size of the local array without the ghost cells */
	lsizes[0] = xdf-1; lsizes[1] = ydf-1;
	/* start of the array in terms of global array indices, ignoring the ghost cells */
	array_start[0] = coord[0] * lsizes[0];
	array_start[1] = coord[1] * lsizes[1];
	
	/* Creating the type for the arrays withou ghost cells
	 * the idea is to read them with a type of array with ghost cells */
	MPI_Type_create_subarray(2, gsizes, lsizes, array_start, MPI_ORDER_C, sfblock, &sfsarray);
	MPI_Type_commit(&sfsarray);
	
	/* size of the local array with the ghost cells */
	memsizes[0] = lsizes[0] + 1;
	memsizes[1] = lsizes[1] + 1;
	
	/* start of the local array with ghost cells: note: ghost cells are the last row and columns */
	array_start[0] = 0;
	array_start[1] = 0;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, sfblock, &sflarray);
	MPI_Type_commit(&sflarray);
	
	lsizes[0] = xdf-1;
	lsizes[1] = 1;
	array_start[0] = 0;
	array_start[1] = 0;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, sfblock, &sfsend_to_dw);
	MPI_Type_commit(&sfsend_to_dw);
	
	array_start[0] = 0;
	array_start[1] = ydf-1 ;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, sfblock, &sfrecv_from_up);
	MPI_Type_commit(&sfrecv_from_up);
	
	/* in the x direction */
	lsizes[0] = 1;
	lsizes[1] = ydf-1;
	array_start[0] = 0;
	array_start[1] = 0;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, sfblock, &sfsend_to_le);
	MPI_Type_commit(&sfsend_to_le);
	
	array_start[0] = xdf-1;
	array_start[1] = 0;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, sfblock, &sfrecv_from_ri);
	MPI_Type_commit(&sfrecv_from_ri);
	
	/* corner */
	lsizes[0] = 1;
	lsizes[1] = 1;
	array_start[0] = 0;
	array_start[1] = 0;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, sfblock, &f_sfsend_to_dl);
	MPI_Type_commit(&f_sfsend_to_dl);
	
	array_start[0] = xdf-1;
	array_start[1] = ydf-1;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, sfblock, &f_sfrecv_from_ur);
	MPI_Type_commit(&f_sfrecv_from_ur);
	
	/* for interpolation */
	/* no. of rows and columns in global array*/ 
	gsizes[0] = (xdc-1)*XDIM;    gsizes[1] = (ydc-1)*YDIM;     
	/* no. of processors in x and y directions */ 
	psizes[0] = XDIM;    psizes[1] = YDIM;
	/* size of the local array without the ghost cells */
	lsizes[0] = xdc-1; lsizes[1] = ydc-1;
	/* start of the array in terms of global array indices, ignoring the ghost cells */
	array_start[0] = coord[0] * lsizes[0];
	array_start[1] = coord[1] * lsizes[1];
	
	MPI_Type_create_subarray(2, gsizes, lsizes, array_start, MPI_ORDER_C, intplelmnt, &intplsarray);
	MPI_Type_commit(&intplsarray);
	
	/* size of the local array with the ghost cells */
	memsizes[0] = lsizes[0] + 6;
	memsizes[1] = lsizes[1] + 6;
	
	/* start of the local array with ghost cells: note: ghost cells are the last row and columns */
	array_start[0] = 2;
	array_start[1] = 2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, intplelmnt, &intpllarray);
	MPI_Type_commit(&intpllarray);
	
	/* data type for transfer in y direction */
	/* For the interpolation */
	lsizes[0] = xdc-1;
	lsizes[1] = 4;
	array_start[0] = 2;
	array_start[1] = 2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, intplelmnt, &intpl_send_to_dw);
	MPI_Type_commit(&intpl_send_to_dw);
	
	array_start[0] = 2;
	array_start[1] = ydc-1+2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, intplelmnt, &intpl_recv_fr_up);
	MPI_Type_commit(&intpl_recv_fr_up);
	
	lsizes[0] = xdc-1;
	lsizes[1] = 2;
	array_start[0] = 2;
	array_start[1] = ydc-1-2 +2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, intplelmnt, &intpl_send_to_up);
	MPI_Type_commit(&intpl_send_to_up);
	
	array_start[0] = 2;
	array_start[1] = 0;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, intplelmnt, &intpl_recv_fr_dw);
	MPI_Type_commit(&intpl_recv_fr_dw);
/*******************************************************************************************************************/
	
	lsizes[0] = 4;
	lsizes[1] = ydc-1;
	array_start[0] = 2;
	array_start[1] = 2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, intplelmnt, &intpl_send_to_le);
	MPI_Type_commit(&intpl_send_to_le);
	
	array_start[0] = xdc-1 + 2;
	array_start[1] = 2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, intplelmnt, &intpl_recv_fr_ri);
	MPI_Type_commit(&intpl_recv_fr_ri);
	
	lsizes[0] = 2;
	lsizes[1] = ydc-1;
	array_start[0] = xdc-1-2 +2;
	array_start[1] = 2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, intplelmnt, &intpl_send_to_ri);
	MPI_Type_commit(&intpl_send_to_ri);
	
	array_start[0] = 0;
	array_start[1] = 2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, intplelmnt, &intpl_recv_fr_le);
	MPI_Type_commit(&intpl_recv_fr_le);
	
	/* Temporary data layout for the interpolation at the corners */
	
	lsizes[0] = 2;
	lsizes[1] = 2;
	array_start[0] = xdc-1-2 +2;
	array_start[1] = ydc-1-2 +2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, intplelmnt, &intpl_send_to_ur);
	MPI_Type_commit(&intpl_send_to_ur);
	
	lsizes[0] = 4;
	lsizes[1] = 4;
	array_start[0] = xdc-1+2;
	array_start[1] = ydc-1+2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, intplelmnt, &intpl_recv_fr_ur);
	MPI_Type_commit(&intpl_recv_fr_ur);
	
	lsizes[0] = 4;
	lsizes[1] = 2;
	array_start[0] = 2;
	array_start[1] = ydc -1 -2 +2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, intplelmnt, &intpl_send_to_ul);
	MPI_Type_commit(&intpl_send_to_ul);
	
	lsizes[0] = 2;
	lsizes[1] = 4;
	array_start[0] = 0;
	array_start[1] = ydc-1+2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, intplelmnt, &intpl_recv_fr_ul);
	MPI_Type_commit(&intpl_recv_fr_ul);
	
	lsizes[0] = 2;
	lsizes[1] = 4;
	array_start[0] = xdc-1-2 +2;
	array_start[1] = 2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, intplelmnt, &intpl_send_to_dr);
	MPI_Type_commit(&intpl_send_to_dr);
	
	lsizes[0] = 4;
	lsizes[1] = 2;
	array_start[0] = xdc -1 +2;
	array_start[1] = 0;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, intplelmnt, &intpl_recv_fr_dr);
	MPI_Type_commit(&intpl_recv_fr_dr);
	
	lsizes[0] = 4;
	lsizes[1] = 4;
	array_start[0] = 2;
	array_start[1] = 2;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, intplelmnt, &intpl_send_to_dl);
	MPI_Type_commit(&intpl_send_to_dl);
	
	lsizes[0] = 2;
	lsizes[1] = 2;
	array_start[0] = 0;
	array_start[1] = 0;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, intplelmnt, &intpl_recv_fr_dl);
	MPI_Type_commit(&intpl_recv_fr_dl);

	t1=MPI_Wtime();

	/* Computing solution parameters */

	Reb=Rebs;
	ubulk=ub;

	/* Computing the Maximum velocity in lattice */
	umax=1.;//ut*(2.5*log(Ret)+5.5);

	/* Computing the required iterations for one non-dimensional time */
	utime=((double)(0.5*(NZs-2.)/CFL))/((double)(GR));
	dt=1./utime;
	utime=(int)utime;

	tau=(3./Reb)*ubulk*(zd-2.)+0.5;
	tau=1./tau;

	/* for fine and coarse grids */
	tauf = (3./Reb)*ubulk*(NZs-2.)+0.5;
	tauc = (1./((double)GR))*(3./Reb)*ubulk*(NZs-2.)+0.5;
	
	tauf = 1./tauf;
	tauc = 1./tauc;
	
//	Fx=CFL*CFL*ut*ut*2./((double)(zd-2));
	Fx=0.;
	Fy=0.;
	Fz=0.;

	/* defining the effective parameters in collision for performance issues */
	Gx=Fx;
	Gy=0.;
	Gz=0.;

	/*Maximum loop iterations */
		nut=TEND-TSTR;
		up=nut*utime;

#ifdef ITER
		up=ITER;
#endif

	/* Memory Alocation */
	/* For coarse grid */
	V.s=(int *)calloc(xdc*ydc*zdc,sizeof(int));
	if (V.s == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.s: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}

	V.f=(DP *)calloc(xdc*ydc*zdc*19,sizeof(DP));
	if (V.f == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.f: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}

	V.ftemp=(DP *)calloc(xdc*ydc*zdc*19,sizeof(DP));
	if (V.ftemp == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.ftemp: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	
	V.rightbufs=(DP *)calloc(5*ydc*zdc,sizeof(DP));
	if (V.rightbufs == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.rightbufs: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	
	V.rightbufr=(DP *)calloc(5*ydc*zdc,sizeof(DP));
	if (V.rightbufr == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.rightbufr: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	
	V.leftbufs=(DP *)calloc(5*ydc*zdc,sizeof(DP));
	if (V.leftbufs == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.leftbufs: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	
	V.leftbufr=(DP *)calloc(5*ydc*zdc,sizeof(DP));
	if (V.leftbufr == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.leftbufr: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
//
	V.upbufs=(DP *)calloc(5*xdc*zdc,sizeof(DP));
	if (V.upbufs == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.upbufs: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	
	V.upbufr=(DP *)calloc(5*xdc*zdc,sizeof(DP));
	if (V.upbufr == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.upbufr: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	
	V.dwbufs=(DP *)calloc(5*xdc*zdc,sizeof(DP));
	if (V.dwbufs == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.dwbufs: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	
	V.dwbufr=(DP *)calloc(5*xdc*zdc,sizeof(DP));
	if (V.dwbufr == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.dwbufr: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
//
	V.urbufs=(DP *)calloc(9*zdc,sizeof(DP));
	if (V.urbufs == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.urbufs: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	V.dlbufr=(DP *)calloc(9*zdc,sizeof(DP));
	if (V.dlbufr == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.dlbufr: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	V.dlbufs=(DP *)calloc(9*zdc,sizeof(DP));
	if (V.dlbufs == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.dlbufs: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	V.urbufr=(DP *)calloc(9*zdc,sizeof(DP));
	if (V.urbufr == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.urbufr: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}

	/* For the bottom fine grid */
	V.sls=(int *)calloc(xdf*ydf*zdf,sizeof(int));
	if (V.sls == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.sls: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}

	V.slf=(DP *)calloc(xdf*ydf*zdf*19,sizeof(DP));
	if (V.slf == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.slf: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}

	V.slftemp=(DP *)calloc(xdf*ydf*zdf*19,sizeof(DP));
	if (V.slftemp == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.slftemp: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	
	V.slrightbufs=(DP *)calloc(5*ydf*zdf,sizeof(DP));
	if (V.slrightbufs == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.slrightbufs: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	
	V.slrightbufr=(DP *)calloc(5*ydf*zdf,sizeof(DP));
	if (V.slrightbufr == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.slrightbufr: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	
	V.slleftbufs=(DP *)calloc(5*ydf*zdf,sizeof(DP));
	if (V.slleftbufs == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.slleftbufs: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	
	V.slleftbufr=(DP *)calloc(5*ydf*zdf,sizeof(DP));
	if (V.slleftbufr == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.slleftbufr: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
//	
	V.slupbufs=(DP *)calloc(5*xdf*zdf,sizeof(DP));
	if (V.slupbufs == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.slupbufs: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	
	V.slupbufr=(DP *)calloc(5*xdf*zdf,sizeof(DP));
	if (V.slupbufr == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.slupbufr: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	
	V.sldwbufs=(DP *)calloc(5*xdf*zdf,sizeof(DP));
	if (V.sldwbufs == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.sldwbufs: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	
	V.sldwbufr=(DP *)calloc(5*xdf*zdf,sizeof(DP));
	if (V.sldwbufr == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.sldwbufr: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
//
	V.slurbufs=(DP *)calloc(9*zdf,sizeof(DP));
	if (V.slurbufs == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.slurbufs: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	V.sldlbufr=(DP *)calloc(9*zdf,sizeof(DP));
	if (V.sldlbufr == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.sldlbufr: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	V.sldlbufs=(DP *)calloc(9*zdf,sizeof(DP));
	if (V.sldlbufs == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.sldlbufs: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	V.slurbufr=(DP *)calloc(9*zdf,sizeof(DP));
	if (V.slurbufr == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.slurbufr: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}

	/* For the upper fine grid */
	V.sus=(int *)calloc(xdf*ydf*zdf,sizeof(int));
	if (V.sus == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.sus: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}

	V.suf=(DP *)calloc(xdf*ydf*zdf*19,sizeof(DP));
	if (V.suf == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.suf: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}

	V.suftemp=(DP *)calloc(xdf*ydf*zdf*19,sizeof(DP));
	if (V.suftemp == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.suftemp: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	
	V.surightbufs=(DP *)calloc(5*ydf*zdf,sizeof(DP));
	if (V.surightbufs == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.surightbufs: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	
	V.surightbufr=(DP *)calloc(5*ydf*zdf,sizeof(DP));
	if (V.surightbufr == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.surightbufr: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	
	V.suleftbufs=(DP *)calloc(5*ydf*zdf,sizeof(DP));
	if (V.suleftbufs == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.suleftbufs: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	
	V.suleftbufr=(DP *)calloc(5*ydf*zdf,sizeof(DP));
	if (V.suleftbufr == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.suleftbufr: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
//
	V.suupbufs=(DP *)calloc(5*xdf*zdf,sizeof(DP));
	if (V.suupbufs == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.suupbufs: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	
	V.suupbufr=(DP *)calloc(5*xdf*zdf,sizeof(DP));
	if (V.suupbufr == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.suupbufr: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	
	V.sudwbufs=(DP *)calloc(5*xdf*zdf,sizeof(DP));
	if (V.sudwbufs == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.sudwbufs: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	
	V.sudwbufr=(DP *)calloc(5*xdf*zdf,sizeof(DP));
	if (V.sudwbufr == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.sudwbufr: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
//
	V.suurbufs=(DP *)calloc(9*zdf,sizeof(DP));
	if (V.suurbufs == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.suurbufs: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	V.sudlbufr=(DP *)calloc(9*zdf,sizeof(DP));
	if (V.sudlbufr == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.sudlbufr: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	V.sudlbufs=(DP *)calloc(9*zdf,sizeof(DP));
	if (V.sudlbufs == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.sudlbufs: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	V.suurbufr=(DP *)calloc(9*zdf,sizeof(DP));
	if (V.suurbufr == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.suurbufr: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}

	/* For interpolation */
	V.lin=(DP *)calloc(23*(ydc+5)*(xdc+5),sizeof(DP));
	if (V.lin == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.lin: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	V.linp=(DP *)calloc(23*(ydc+5)*(xdc+5),sizeof(DP));
	if (V.linp == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.linp: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	V.linpp=(DP *)calloc(23*(ydc+5)*(xdc+5),sizeof(DP));
	if (V.linpp == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.linpp: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	V.uin=(DP *)calloc(23*(ydc+5)*(xdc+5),sizeof(DP));
	if (V.uin == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.uin: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	V.uinp=(DP *)calloc(23*(ydc+5)*(xdc+5),sizeof(DP));
	if (V.uinp == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.uinp: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	V.uinpp=(DP *)calloc(23*(ydc+5)*(xdc+5),sizeof(DP));
	if (V.uinpp == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.uinpp: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	
	V.lint=(DP *)calloc(23*(ydc+5)*(xdc+5)*4,sizeof(DP));
	if (V.lint == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.lint: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	V.uint=(DP *)calloc(23*(ydc+5)*(xdc+5)*4,sizeof(DP));
	if (V.uint == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.uint: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	
	for (i=0;i<12;i++)
		V.req[i] = MPI_REQUEST_NULL;

	for (i=0;i<12;i++)
		V.slreq[i] = MPI_REQUEST_NULL;

	for (i=0;i<12;i++)
		V.sureq[i] = MPI_REQUEST_NULL;
	/*Two first req are for send (left-right) and the second two are for recv (left-right)
	so even numbers are for right and odd for left */

	/* Establishing a channel for communication on coarse grid */
	MPI_Send_init(V.leftbufs,(5*ydc*zdc),MPDP,pd.left,tagl,MPI_COMM_WORLD,&V.req[0]);
	MPI_Send_init(V.rightbufs,(5*ydc*zdc),MPDP,pd.right,tagr,MPI_COMM_WORLD,&V.req[1]);
	MPI_Recv_init(V.leftbufr,(5*ydc*zdc),MPDP,pd.left,tagr,MPI_COMM_WORLD,&V.req[2]);
	MPI_Recv_init(V.rightbufr,(5*ydc*zdc),MPDP,pd.right,tagl,MPI_COMM_WORLD,&V.req[3]);
	
	MPI_Send_init(V.dwbufs,(5*xdc*zdc),MPDP,pd.dw,tagl,MPI_COMM_WORLD,&V.req[4]);
	MPI_Send_init(V.upbufs,(5*xdc*zdc),MPDP,pd.up,tagr,MPI_COMM_WORLD,&V.req[5]);
	MPI_Recv_init(V.dwbufr,(5*xdc*zdc),MPDP,pd.dw,tagr,MPI_COMM_WORLD,&V.req[6]);
	MPI_Recv_init(V.upbufr,(5*xdc*zdc),MPDP,pd.up,tagl,MPI_COMM_WORLD,&V.req[7]);
	
	MPI_Send_init(V.urbufs,9*zdc,MPDP,pd.ur,tagl,MPI_COMM_WORLD,&V.req[8]);
	MPI_Recv_init(V.dlbufr,9*zdc,MPDP,pd.dl,tagl,MPI_COMM_WORLD,&V.req[9]);
	MPI_Send_init(V.dlbufs,9*zdc,MPDP,pd.dl,tagr,MPI_COMM_WORLD,&V.req[10]);
	MPI_Recv_init(V.urbufr,9*zdc,MPDP,pd.ur,tagr,MPI_COMM_WORLD,&V.req[11]);

	/* Establishing a channel for communication on the upper grid */
	MPI_Send_init(V.slleftbufs,(5*ydf*zdf),MPDP,pd.left,tagl,MPI_COMM_WORLD,&V.slreq[0]);
	MPI_Send_init(V.slrightbufs,(5*ydf*zdf),MPDP,pd.right,tagr,MPI_COMM_WORLD,&V.slreq[1]);
	MPI_Recv_init(V.slleftbufr,(5*ydf*zdf),MPDP,pd.left,tagr,MPI_COMM_WORLD,&V.slreq[2]);
	MPI_Recv_init(V.slrightbufr,(5*ydf*zdf),MPDP,pd.right,tagl,MPI_COMM_WORLD,&V.slreq[3]);
	
	MPI_Send_init(V.sldwbufs,(5*xdf*zdf),MPDP,pd.dw,tagl,MPI_COMM_WORLD,&V.slreq[4]);
	MPI_Send_init(V.slupbufs,(5*xdf*zdf),MPDP,pd.up,tagr,MPI_COMM_WORLD,&V.slreq[5]);
	MPI_Recv_init(V.sldwbufr,(5*xdf*zdf),MPDP,pd.dw,tagr,MPI_COMM_WORLD,&V.slreq[6]);
	MPI_Recv_init(V.slupbufr,(5*xdf*zdf),MPDP,pd.up,tagl,MPI_COMM_WORLD,&V.slreq[7]);
	
	MPI_Send_init(V.slurbufs,9*zdf,MPDP,pd.ur,tagl,MPI_COMM_WORLD,&V.slreq[8]);
	MPI_Recv_init(V.sldlbufr,9*zdf,MPDP,pd.dl,tagl,MPI_COMM_WORLD,&V.slreq[9]);
	MPI_Send_init(V.sldlbufs,9*zdf,MPDP,pd.dl,tagr,MPI_COMM_WORLD,&V.slreq[10]);
	MPI_Recv_init(V.slurbufr,9*zdf,MPDP,pd.ur,tagr,MPI_COMM_WORLD,&V.slreq[11]);

	/* Establishing a channel for communication on the bottom grid */
	MPI_Send_init(V.suleftbufs,(5*ydf*zdf),MPDP,pd.left,tagl,MPI_COMM_WORLD,&V.sureq[0]);
	MPI_Send_init(V.surightbufs,(5*ydf*zdf),MPDP,pd.right,tagr,MPI_COMM_WORLD,&V.sureq[1]);
	MPI_Recv_init(V.suleftbufr,(5*ydf*zdf),MPDP,pd.left,tagr,MPI_COMM_WORLD,&V.sureq[2]);
	MPI_Recv_init(V.surightbufr,(5*ydf*zdf),MPDP,pd.right,tagl,MPI_COMM_WORLD,&V.sureq[3]);
	
	MPI_Send_init(V.sudwbufs,(5*xdf*zdf),MPDP,pd.dw,tagl,MPI_COMM_WORLD,&V.sureq[4]);
	MPI_Send_init(V.suupbufs,(5*xdf*zdf),MPDP,pd.up,tagr,MPI_COMM_WORLD,&V.sureq[5]);
	MPI_Recv_init(V.sudwbufr,(5*xdf*zdf),MPDP,pd.dw,tagr,MPI_COMM_WORLD,&V.sureq[6]);
	MPI_Recv_init(V.suupbufr,(5*xdf*zdf),MPDP,pd.up,tagl,MPI_COMM_WORLD,&V.sureq[7]);
	
	MPI_Send_init(V.suurbufs,9*zdf,MPDP,pd.ur,tagl,MPI_COMM_WORLD,&V.sureq[8]);
	MPI_Recv_init(V.sudlbufr,9*zdf,MPDP,pd.dl,tagl,MPI_COMM_WORLD,&V.sureq[9]);
	MPI_Send_init(V.sudlbufs,9*zdf,MPDP,pd.dl,tagr,MPI_COMM_WORLD,&V.sureq[10]);
	MPI_Recv_init(V.suurbufr,9*zdf,MPDP,pd.ur,tagr,MPI_COMM_WORLD,&V.sureq[11]);

	ch=argv[1][1];
	mpi_errors=MPI_Bcast(&ch,1,MPI_CHAR,0,MPI_COMM_WORLD);
// 	if (mpi_errors != MPI_SUCCESS) 
// 	{
// 		MPI_Error_class(mpi_errors,&mpi_errors_class);
// 		MPI_Error_string(mpi_errors_class,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error class %s reported by rank %d\n",mpi_errors_message,pd.myrank);
// 		MPI_Error_string(mpi_errors,mpi_errors_message,&mpi_errors_length);
// 		fprintf(stderr,"MPI Error on rank %d: %s by MPI_Bcast, in line %d of function %s\n",pd.myrank,mpi_errors_message,__LINE__,__FUNCTION__);
// 		MPI_Abort(MPI_COMM_WORLD,mpi_errors);	/* abort*/
// 	}
	MPI_Barrier(MPI_COMM_WORLD);

	t3=MPI_Wtime();
	ss = 0;
	switch (ch)
	{
		case 's':
		{
			TM=TSTR;
			V=solid_init(xdc,ydc,zdc, xdf,ydf,zdf, V);
#ifdef RSLIP
			V=slip_ridges_init(xdf, ydf, zdf, V, pd);
#endif
#ifdef PSLIP
			V=slip_posts_init(xdf, ydf, zdf, V, pd);
#endif
#ifdef SSLIP
			V=slip_spanwise_init(xdf, ydf, zdf, V, pd);
#endif
			V=init(xdc,ydc,zdc, xdf,ydf,zdf,V);
			ss=0;
			break;
		}
		case 'i':
		{
			ss=0;
			TM=TSTR;
			cntr=TSTR;
			V=init(xdc,ydc,zdc, xdf,ydf,zdf,V);
			fprintf(stderr,"starting to read and interpolate the data %d\n",pd.myrank);
			V=solid_init(xdc,ydc,zdc,xdf,ydf,zdf,V);
#ifdef RSLIP
			V=slip_ridges_init(xdf, ydf, zdf, V, pd);
#endif
#ifdef PSLIP
			V=slip_posts_init(xdf, ydf, zdf, V, pd);
#endif
#ifdef SSLIP
			V=slip_spanwise_init(xdf, ydf, zdf, V, pd);
#endif
			V=initializer(xdf,ydf,zdf, tauf,xdc,ydc,zdc, tauc,xd,yd,zd, tau, GR, V, pd, cntr);
			V=inter_init(xdc, ydc, zdc, V, pd);
			break;
		}
		case 'p':
		{
			if (pd.myrank==0)
				fprintf(stderr,"TRYING TO READ DATA FROM FILES\n");
			cntr=TSTR;
			V=init(xdc,ydc,zdc, xdf,ydf,zdf,V);
			
			sprintf(fn,"lbf-c.%.4d",cntr); 
			printf("reading %s \n",fn);
			MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
			MPI_File_set_view(fh, 0, MPI_DOUBLE, sarray, "native", MPI_INFO_NULL);
			MPI_File_read_all(fh, V.ftemp, 1, larray, &status); 
			MPI_Get_count(&status,MPI_DOUBLE,&i);
			printf("V.ftemp %d\n",i);
			MPI_File_set_view(fh, (xdc-1)*XDIM*(ydc-1)*YDIM*zdc*19*sizeof(double), MPI_DOUBLE, sarray, "native", MPI_INFO_NULL);
			MPI_File_read_all(fh, V.f, 1, larray, &status); 
			MPI_Get_count(&status,MPI_DOUBLE,&i);
			printf("V.f %d\n",i);
			MPI_File_close(&fh);
			
			MPI_Sendrecv(V.ftemp,1,send_to_le,pd.left,1129,V.ftemp,1,recv_from_ri,pd.right,1129,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.f,1,send_to_le,pd.left,1128,V.f,1,recv_from_ri,pd.right,1128,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.ftemp,1,send_to_dw,pd.dw,4130,V.ftemp,1,recv_from_up,pd.up,4130,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.f,1,send_to_dw,pd.dw,4131,V.f,1,recv_from_up,pd.up,4131,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.ftemp,1,f_send_to_dl,pd.dl,4132,V.ftemp,1,f_recv_from_ur,pd.ur,4132,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.f,1,f_send_to_dl,pd.dl,4133,V.f,1,f_recv_from_ur,pd.ur,4133,MPI_COMM_WORLD,&status);
			
			sprintf(fn,"lbf-slf.%.4d",cntr); 
			MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
			MPI_File_set_view(fh, 0, MPI_DOUBLE, sfsarray, "native", MPI_INFO_NULL);
			MPI_File_read_all(fh, V.slftemp, 1, sflarray, &status); 
			MPI_Get_count(&status,MPI_DOUBLE,&i);
			printf("V.slftemp %d\n",i);
			MPI_File_set_view(fh, (xdf-1)*XDIM*(ydf-1)*YDIM*zdf*19*sizeof(double), MPI_DOUBLE, sfsarray, "native", MPI_INFO_NULL);
			MPI_File_read_all(fh, V.slf, 1, sflarray, &status); 
			MPI_Get_count(&status,MPI_DOUBLE,&i);
			printf("V.slf %d\n",i);
			MPI_File_close(&fh);

			MPI_Sendrecv(V.slftemp,1,sfsend_to_le,pd.left,1130,V.slftemp,1,sfrecv_from_ri,pd.right,1130,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.slf,1,sfsend_to_le,pd.left,1131,V.slf,1,sfrecv_from_ri,pd.right,1131,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.slftemp,1,sfsend_to_dw,pd.dw,4230,V.slftemp,1,sfrecv_from_up,pd.up,4230,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.slf,1,sfsend_to_dw,pd.dw,4231,V.slf,1,sfrecv_from_up,pd.up,4231,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.slftemp,1,f_sfsend_to_dl,pd.dl,4232,V.slftemp,1,f_sfrecv_from_ur,pd.ur,4232,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.slf,1,f_sfsend_to_dl,pd.dl,4233,V.slf,1,f_sfrecv_from_ur,pd.ur,4233,MPI_COMM_WORLD,&status);
			
			
			sprintf(fn,"lbf-suf.%.4d",cntr); 
			MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
			MPI_File_set_view(fh, 0, MPI_DOUBLE, sfsarray, "native", MPI_INFO_NULL);
			MPI_File_read_all(fh, V.suftemp, 1, sflarray, &status); 
			MPI_Get_count(&status,MPI_DOUBLE,&i);
			printf("V.suftemp %d\n",i);
			MPI_File_set_view(fh, (xdf-1)*XDIM*(ydf-1)*YDIM*zdf*19*sizeof(double), MPI_DOUBLE, sfsarray, "native", MPI_INFO_NULL);
			MPI_File_read_all(fh, V.suf, 1, sflarray, &status); 
			MPI_Get_count(&status,MPI_DOUBLE,&i);
			printf("V.suf %d\n",i);
			MPI_File_close(&fh);

			MPI_Sendrecv(V.suftemp,1,sfsend_to_le,pd.left,1930,V.suftemp,1,sfrecv_from_ri,pd.right,1930,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.suf,1,sfsend_to_le,pd.left,1931,V.suf,1,sfrecv_from_ri,pd.right,1931,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.suftemp,1,sfsend_to_dw,pd.dw,1932,V.suftemp,1,sfrecv_from_up,pd.up,1932,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.suf,1,sfsend_to_dw,pd.dw,1933,V.suf,1,sfrecv_from_up,pd.up,1933,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.suftemp,1,f_sfsend_to_dl,pd.dl,1934,V.suftemp,1,f_sfrecv_from_ur,pd.ur,1934,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.suf,1,f_sfsend_to_dl,pd.dl,1935,V.suf,1,f_sfrecv_from_ur,pd.ur,1935,MPI_COMM_WORLD,&status);

			if (pd.myrank==0)
			{
				sprintf(fn,"force.%.4d",cntr);
				sv=fopen(fn,"rb");
				i=fread(&Fx,sizeof(DP),1,sv);
				V.force[0] = Fx;
				fclose(sv);
				printf("force read %d\n",i);
			}
			MPI_Bcast(&Fx,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
			V.force[0] = Fx;
			printf("reading V.force=%.12f myrank=%d\n",V.force[0],pd.myrank);
			
			sprintf(fn,"lbf-uin.%.4d",cntr);
			MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
			MPI_File_set_view(fh, 0, MPI_DOUBLE, intplsarray, "native", MPI_INFO_NULL);
			MPI_File_read_all(fh, V.uin, 1, intpllarray, &status); 
			MPI_Get_count(&status,MPI_DOUBLE,&i);
			printf("V.uin %d\n",i);
			MPI_File_set_view(fh, (xdc-1)*(ydc-1)*23*YDIM*XDIM*sizeof(double), MPI_DOUBLE, intplsarray, "native", MPI_INFO_NULL);
			MPI_File_read_all(fh, V.uinp, 1, intpllarray, &status); 
			MPI_Get_count(&status,MPI_DOUBLE,&i);
			printf("V.uinp %d\n",i);
			MPI_File_set_view(fh, 2*(xdc-1)*(ydc-1)*23*YDIM*XDIM*sizeof(double), MPI_DOUBLE, intplsarray, "native", MPI_INFO_NULL);
			MPI_File_read_all(fh, V.uinpp, 1, intpllarray, &status); 
			MPI_Get_count(&status,MPI_DOUBLE,&i);
			printf("V.uinpp %d\n",i);
			MPI_File_close(&fh);

			MPI_Sendrecv(V.uin,1,intpl_send_to_ri,pd.right,1100,V.uin,1,intpl_recv_fr_le,pd.left,1100,MPI_COMM_WORLD,&tstatus);
			MPI_Sendrecv(V.uin,1,intpl_send_to_le,pd.left,1200,V.uin,1,intpl_recv_fr_ri,pd.right,1200,MPI_COMM_WORLD,&tstatus);
			MPI_Sendrecv(V.uin,1,intpl_send_to_dw,pd.dw,1300,V.uin,1,intpl_recv_fr_up,pd.up,1300,MPI_COMM_WORLD,&tstatus);
			MPI_Sendrecv(V.uin,1,intpl_send_to_up,pd.up,1400,V.uin,1,intpl_recv_fr_dw,pd.dw,1400,MPI_COMM_WORLD,&tstatus);
			MPI_Sendrecv(V.uin,1,intpl_send_to_ur,pd.ur,1500,V.uin,1,intpl_recv_fr_dl,pd.dl,1500,MPI_COMM_WORLD,&tstatus);
			MPI_Sendrecv(V.uin,1,intpl_send_to_ul,pd.ul,1600,V.uin,1,intpl_recv_fr_dr,pd.dr,1600,MPI_COMM_WORLD,&tstatus);
			MPI_Sendrecv(V.uin,1,intpl_send_to_dl,pd.dl,1700,V.uin,1,intpl_recv_fr_ur,pd.ur,1700,MPI_COMM_WORLD,&tstatus);
			MPI_Sendrecv(V.uin,1,intpl_send_to_dr,pd.dr,1800,V.uin,1,intpl_recv_fr_ul,pd.ul,1800,MPI_COMM_WORLD,&tstatus);

			MPI_Sendrecv(V.uinp,1,intpl_send_to_dw,pd.dw,5130,V.uinp,1,intpl_recv_fr_up,pd.up,5130,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.uinp,1,intpl_send_to_up,pd.up,5131,V.uinp,1,intpl_recv_fr_dw,pd.dw,5131,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.uinp,1,intpl_send_to_le,pd.left,1416,V.uinp,1,intpl_recv_fr_ri,pd.right,1416,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.uinp,1,intpl_send_to_ri,pd.right,1417,V.uinp,1,intpl_recv_fr_le,pd.left,1417,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.uinp,1,intpl_send_to_dl,pd.dl,6134,V.uinp,1,intpl_recv_fr_ur,pd.ur,6134,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.uinp,1,intpl_send_to_dr,pd.dr,6135,V.uinp,1,intpl_recv_fr_ul,pd.ul,6135,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.uinp,1,intpl_send_to_ul,pd.ul,6136,V.uinp,1,intpl_recv_fr_dr,pd.dr,6136,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.uinp,1,intpl_send_to_ur,pd.ur,6137,V.uinp,1,intpl_recv_fr_dl,pd.dl,6137,MPI_COMM_WORLD,&status);

			MPI_Sendrecv(V.uinpp,1,intpl_send_to_dw,pd.dw,5130,V.uinpp,1,intpl_recv_fr_up,pd.up,5130,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.uinpp,1,intpl_send_to_up,pd.up,5131,V.uinpp,1,intpl_recv_fr_dw,pd.dw,5131,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.uinpp,1,intpl_send_to_le,pd.left,1418,V.uinpp,1,intpl_recv_fr_ri,pd.right,1418,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.uinpp,1,intpl_send_to_ri,pd.right,1419,V.uinpp,1,intpl_recv_fr_le,pd.left,1419,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.uinpp,1,intpl_send_to_dl,pd.dl,6138,V.uinpp,1,intpl_recv_fr_ur,pd.ur,6138,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.uinpp,1,intpl_send_to_dr,pd.dr,6139,V.uinpp,1,intpl_recv_fr_ul,pd.ul,6139,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.uinpp,1,intpl_send_to_ul,pd.ul,6140,V.uinpp,1,intpl_recv_fr_dr,pd.dr,6140,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.uinpp,1,intpl_send_to_ur,pd.ur,6141,V.uinpp,1,intpl_recv_fr_dl,pd.dl,6141,MPI_COMM_WORLD,&status);
			
			sprintf(fn,"lbf-lin.%.4d",cntr);
			MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
			MPI_File_set_view(fh, 0, MPI_DOUBLE, intplsarray, "native", MPI_INFO_NULL);
			MPI_File_read_all(fh, V.lin, 1, intpllarray, &status); 
			MPI_Get_count(&status,MPI_DOUBLE,&i);
			printf("V.lin %d\n",i);
			MPI_File_set_view(fh, (xdc-1)*(ydc-1)*23*YDIM*XDIM*sizeof(double), MPI_DOUBLE, intplsarray, "native", MPI_INFO_NULL);
			MPI_File_read_all(fh, V.linp, 1, intpllarray, &status); 
			MPI_Get_count(&status,MPI_DOUBLE,&i);
			printf("V.linp %d\n",i);
			MPI_File_set_view(fh, 2*(xdc-1)*(ydc-1)*23*YDIM*XDIM*sizeof(double), MPI_DOUBLE, intplsarray, "native", MPI_INFO_NULL);
			MPI_File_read_all(fh, V.linpp, 1, intpllarray, &status); 
			MPI_Get_count(&status,MPI_DOUBLE,&i);
			printf("V.linpp %d\n",i);
			MPI_File_close(&fh);
			
			MPI_Sendrecv(V.lin,1,intpl_send_to_ri,pd.right,1101,V.lin,1,intpl_recv_fr_le,pd.left,1101,MPI_COMM_WORLD,&tstatus);
			MPI_Sendrecv(V.lin,1,intpl_send_to_le,pd.left,1202,V.lin,1,intpl_recv_fr_ri,pd.right,1202,MPI_COMM_WORLD,&tstatus);
			MPI_Sendrecv(V.lin,1,intpl_send_to_dw,pd.dw,1303,V.lin,1,intpl_recv_fr_up,pd.up,1303,MPI_COMM_WORLD,&tstatus);
			MPI_Sendrecv(V.lin,1,intpl_send_to_up,pd.up,1404,V.lin,1,intpl_recv_fr_dw,pd.dw,1404,MPI_COMM_WORLD,&tstatus);
			MPI_Sendrecv(V.lin,1,intpl_send_to_ur,pd.ur,1505,V.lin,1,intpl_recv_fr_dl,pd.dl,1505,MPI_COMM_WORLD,&tstatus);
			MPI_Sendrecv(V.lin,1,intpl_send_to_ul,pd.ul,1606,V.lin,1,intpl_recv_fr_dr,pd.dr,1606,MPI_COMM_WORLD,&tstatus);
			MPI_Sendrecv(V.lin,1,intpl_send_to_dl,pd.dl,1707,V.lin,1,intpl_recv_fr_ur,pd.ur,1707,MPI_COMM_WORLD,&tstatus);
			MPI_Sendrecv(V.lin,1,intpl_send_to_dr,pd.dr,1808,V.lin,1,intpl_recv_fr_ul,pd.ul,1808,MPI_COMM_WORLD,&tstatus);

			MPI_Sendrecv(V.linp,1,intpl_send_to_dw,pd.dw,6130,V.linp,1,intpl_recv_fr_up,pd.up,6130,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.linp,1,intpl_send_to_up,pd.up,6131,V.linp,1,intpl_recv_fr_dw,pd.dw,6131,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.linp,1,intpl_send_to_le,pd.left,1422,V.linp,1,intpl_recv_fr_ri,pd.right,1422,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.linp,1,intpl_send_to_ri,pd.right,1423,V.linp,1,intpl_recv_fr_le,pd.left,1423,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.linp,1,intpl_send_to_dl,pd.dl,7134,V.linp,1,intpl_recv_fr_ur,pd.ur,7134,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.linp,1,intpl_send_to_dr,pd.dr,7135,V.linp,1,intpl_recv_fr_ul,pd.ul,7135,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.linp,1,intpl_send_to_ul,pd.ul,7136,V.linp,1,intpl_recv_fr_dr,pd.dr,7136,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.linp,1,intpl_send_to_ur,pd.ur,7137,V.linp,1,intpl_recv_fr_dl,pd.dl,7137,MPI_COMM_WORLD,&status);

			MPI_Sendrecv(V.linpp,1,intpl_send_to_dw,pd.dw,6130,V.linpp,1,intpl_recv_fr_up,pd.up,6130,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.linpp,1,intpl_send_to_up,pd.up,6131,V.linpp,1,intpl_recv_fr_dw,pd.dw,6131,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.linpp,1,intpl_send_to_le,pd.left,1424,V.linpp,1,intpl_recv_fr_ri,pd.right,1424,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.linpp,1,intpl_send_to_ri,pd.right,1425,V.linpp,1,intpl_recv_fr_le,pd.left,1425,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.linpp,1,intpl_send_to_dl,pd.dl,7138,V.linpp,1,intpl_recv_fr_ur,pd.ur,7138,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.linpp,1,intpl_send_to_dr,pd.dr,7139,V.linpp,1,intpl_recv_fr_ul,pd.ul,7139,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.linpp,1,intpl_send_to_ul,pd.ul,7140,V.linpp,1,intpl_recv_fr_dr,pd.dr,7140,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.linpp,1,intpl_send_to_ur,pd.ur,7141,V.linpp,1,intpl_recv_fr_dl,pd.dl,7141,MPI_COMM_WORLD,&status);

			V=solid_init(xdc,ydc,zdc, xdf, ydf, zdf ,V);
#ifdef RSLIP
			V=slip_ridges_init(xdf, ydf, zdf, V, pd);
#endif
#ifdef PSLIP
			V=slip_posts_init(xdf, ydf, zdf, V, pd);
#endif
#ifdef SSLIP
			V=slip_spanwise_init(xdf, ydf, zdf, V, pd);
#endif
			if (pd.myrank==0)
			{
				fprintf(stderr,"END OF READING DATA FROM FILES\n");
				fprintf(stderr,"STARTING FROM %d\n",cntr);
			}
			TM=cntr;
			ss=0;
			//ss=ss*utime;
			break;
		}
		default :
		{
			perror("not a good file name");
			MPI_Finalize();
			return(-1);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	V.rho[0]=0.;rho=0.;rhozero=0.;
	for (i=0;i<(xdf-1);i++)
	{
		for (j=0;j<(ydf-1);j++)
		{
			for (k=0;k<zdf;k++)
			{
				rhozero=0.;
				for (a=0;a<19;a++)
					rhozero += Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);	
				if (k==(zdf-1))
					coef = 0.5;
				else
					coef = 1.;
				rho+= (coef*rhozero);
			}
		}
	}
	for (i=0;i<(xdc-1);i++)
	{
		for (j=0;j<(ydc-1);j++)
		{
			for (k=1;k<(zdc-1);k++)
			{
				rhozero=0.;
				for (a=0;a<19;a++)
					rhozero += Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);	
				if ((k==1)||(k==(zdc-2)))
					coef=0.5*GR*GR*GR;
				else
					coef=1.*GR*GR*GR;
				rho+=(coef*rhozero);
			}
		}
	}
	for (i=0;i<(xdf-1);i++)
	{
		for (j=0;j<(ydf-1);j++)
		{
			for (k=0;k<zdf;k++)
			{
				rhozero=0.;
				for (a=0;a<19;a++)
					rhozero += Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
				if (k==0)
					coef = 0.5;
				else
					coef = 1.;
				rho+= (coef*rhozero);
			}
		}
	}
	V.rho[0]=rho;
	printf("V.rho = %f xdf=%d ydf=%d zdf=%d NZs=%d GR=%d\n",V.rho[0],xdf,ydf,zdf,NZs,GR);
	rho /= ((double)((xdf-1.)*(ydf-1.)*XDIM*YDIM*(NZs-2.)));
	MPI_Allreduce(&rho,&rhozero,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    	V.rho[0]=rhozero;
	printf("Average Density in the channel at start: %.12f  myrank=%d\n",V.rho[0],pd.myrank);

	t4=MPI_Wtime();

 	for (ss;ss<(up+1);ss++)
	{
//		V=encalc(xdc,ydc,zdc, xdf,ydf,zdf, GR, V, pd, TM,ss,ubulk);
		V=encalc2(xdc,ydc,zdc, xdf,ydf,zdf, GR, V, pd, TM,ss,ubulk);
		V=wshearstr(xdf, ydf, zdf, zdc, GR, V, pd, tauf, TM, ss);

		if (!(ss%(int)utime))
		{
			if (!(ss%(50*(int)utime)))
			{
				sprintf(fn,"lbf-c.%.4ld",(TM+(ss/(int)utime))); 
				MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_CREATE | MPI_MODE_WRONLY,info,&fh);
				MPI_File_set_view(fh, 0, MPI_DOUBLE, sarray, "native", info);
				MPI_File_write_all(fh, V.ftemp, 1, larray, &status); 
				MPI_Get_count(&status,MPI_DOUBLE,&i);
				printf("V.ftemp %d\n",i);
				MPI_File_set_view(fh, (xdc-1)*XDIM*(ydc-1)*YDIM*zdc*19*sizeof(double), MPI_DOUBLE, sarray, "native", info);
				MPI_File_write_all(fh, V.f, 1, larray, &status); 
				MPI_Get_count(&status,MPI_DOUBLE,&i);
				printf("V.f %d\n",i);
				MPI_File_close(&fh);
				
				sprintf(fn,"lbf-slf.%.4ld",(TM+(ss/(int)utime))); 
				MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_CREATE | MPI_MODE_WRONLY,info,&fh);
				MPI_File_set_view(fh, 0, MPI_DOUBLE, sfsarray, "native", info);
				MPI_File_write_all(fh, V.slftemp, 1, sflarray, &status); 
				MPI_Get_count(&status,MPI_DOUBLE,&i);
				printf("V.slftemp %d\n",i);
				MPI_File_set_view(fh, (xdf-1)*XDIM*(ydf-1)*YDIM*zdf*19*sizeof(double), MPI_DOUBLE, sfsarray, "native", info);
				MPI_File_write_all(fh, V.slf, 1, sflarray, &status); 
				MPI_Get_count(&status,MPI_DOUBLE,&i);
				printf("V.slf %d\n",i);
				MPI_File_close(&fh);
				
				sprintf(fn,"lbf-suf.%.4ld",(TM+(ss/(int)utime))); 
				MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_CREATE | MPI_MODE_WRONLY,info,&fh);
				MPI_File_set_view(fh, 0, MPI_DOUBLE, sfsarray, "native", info);
				MPI_File_write_all(fh, V.suftemp, 1, sflarray, &status); 
				MPI_Get_count(&status,MPI_DOUBLE,&i);
				printf("V.suftemp %d\n",i);
				MPI_File_set_view(fh, (xdf-1)*XDIM*(ydf-1)*YDIM*zdf*19*sizeof(double), MPI_DOUBLE, sfsarray, "native", info);
				MPI_File_write_all(fh, V.suf, 1, sflarray, &status); 
				MPI_Get_count(&status,MPI_DOUBLE,&i);
				printf("V.suf %d\n",i);
				MPI_File_close(&fh);
				
				if (pd.myrank==0)
				{
					sprintf(fn,"force.%.4ld",(TM+(ss/(int)utime)));
					sv=fopen(fn,"wb");
					i=fwrite(&V.force[0],sizeof(DP),1,sv);
					fclose(sv);
				}
				MPI_Bcast(&i,1,MPI_INT,0,MPI_COMM_WORLD);
				printf("writing V.force=%.12f num written=%d myrank=%d   V.rho=%.12f\n",V.force[0],i,pd.myrank,V.rho[0]);
				
				sprintf(fn,"lbf-uin.%.4ld",(TM+(ss/(int)utime)));
				MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_CREATE | MPI_MODE_WRONLY,info,&fh);
				MPI_File_set_view(fh, 0, MPI_DOUBLE, intplsarray, "native", info);
				MPI_File_write_all(fh, V.uin, 1, intpllarray, &status); 
				MPI_Get_count(&status,MPI_DOUBLE,&i);
				printf("V.uin %d\n",i);
				MPI_File_set_view(fh, (xdc-1)*(ydc-1)*23*YDIM*XDIM*sizeof(double), MPI_DOUBLE, intplsarray, "native", info);
				MPI_File_write_all(fh, V.uinp, 1, intpllarray, &status); 
				MPI_Get_count(&status,MPI_DOUBLE,&i);
				printf("V.uinp %d\n",i);
				MPI_File_set_view(fh, 2*(xdc-1)*(ydc-1)*23*YDIM*XDIM*sizeof(double), MPI_DOUBLE, intplsarray, "native", info);
				MPI_File_write_all(fh, V.uinpp, 1, intpllarray, &status); 
				MPI_Get_count(&status,MPI_DOUBLE,&i);
				printf("V.uinpp %d\n",i);
				MPI_File_close(&fh);
				
				sprintf(fn,"lbf-lin.%.4ld",(TM+(ss/(int)utime)));
				MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_CREATE | MPI_MODE_WRONLY,info,&fh);
				MPI_File_set_view(fh, 0, MPI_DOUBLE, intplsarray, "native", info);
				MPI_File_write_all(fh, V.lin, 1, intpllarray, &status); 
				MPI_Get_count(&status,MPI_DOUBLE,&i);
				printf("V.lin %d\n",i);
				MPI_File_set_view(fh, (xdc-1)*(ydc-1)*23*YDIM*XDIM*sizeof(double), MPI_DOUBLE, intplsarray, "native", info);
				MPI_File_write_all(fh, V.linp, 1, intpllarray, &status); 
				MPI_Get_count(&status,MPI_DOUBLE,&i);
				printf("V.linp %d\n",i);
				MPI_File_set_view(fh, 2*(xdc-1)*(ydc-1)*23*YDIM*XDIM*sizeof(double), MPI_DOUBLE, intplsarray, "native", info);
				MPI_File_write_all(fh, V.linpp, 1, intpllarray, &status); 
				MPI_Get_count(&status,MPI_DOUBLE,&i);
				printf("V.linpp %d\n",i);
				MPI_File_close(&fh);
			}
			V=vel_print( xdf, ydf, zdf, xdc, ydc, zdc, V, (TM+(ss/(int)utime)));
			V=statistics(xdf, ydf, zdf, xdc, ydc, zdc, V, pd, (TM+(ss/(int)utime)), GR);
		}
		V=xcomputations(xdc,ydc,zdc,GR,V,pd,Gx,tauc,ubulk);
		iter = 1;
		for (iter =1; iter<GR; iter++)
		{
 			V=xfcomputations(xdf, ydf, zdf, (1./((double)GR)), V, pd, tauf, ubulk);
//			V=tinterpolate_ln(xdc,ydc,zdc,V,(iter),GR);
//			V=tinterpolate_qd(xdc,ydc,zdc,V,(iter),GR);
			V=tinterpolate_hrmt(xdc,ydc,zdc,V,(iter),GR);
			derivative(xdc,ydc,V.uint);
			derivative(xdc,ydc,V.lint);
			V=ctfdtransfer_opt(xdc,ydc, zdc, xdf, ydf, zdf, (iter), tauc, tauf, GR, V, pd);
			V=fexchange(V);
//			V=fine_force(xdf,ydf,zdf,zdc,GR,V,pd,ubulk);
		}
		V=xfcomputations(xdf, ydf, zdf, (1./((double)GR)), V, pd, tauf, ubulk);
//		V=tinterpolate_ln(xdc,ydc,zdc,V,GR,GR);
//		V=tinterpolate_qd(xdc,ydc,zdc,V,GR,GR);
		V=tinterpolate_hrmt(xdc,ydc,zdc,V,GR,GR);
		derivative(xdc,ydc,V.uint);
		derivative(xdc,ydc,V.lint);
		V=ctfdtransfer_opt(xdc,ydc, zdc, xdf, ydf, zdf, GR, tauc, tauf, GR, V, pd);
#ifdef NOFLTR
		V=ftcdtransfer(xdc,ydc, zdc, xdf, ydf, zdf, GR, tauc, tauf, GR, V, pd);
#endif
#ifdef FLTR
//		V=ftcdtransfer_flt(xdf, ydf, zdf, xdc, ydc, zdc, GR, tauc, tauf, GR, V, pd);
		V=ftcdtransfer_lesfltr(xdf, ydf, zdf, xdc, ydc, zdc, GR, tauc, tauf, GR, V, pd);
#endif
		V=exchange(V);
		V=fexchange(V);
	}
	
	t2=MPI_Wtime();
	V=output_twoD(xdf,ydf,zdf,xdc,ydc,zdc,gratio,V,pd);
	V=output_twoD2(xdf,ydf,zdf,xdc,ydc,zdc,gratio,V,pd);
	if (!pd.myrank)
	{
		tm=fopen("report.code","w");
		fprintf(tm,"xdc=%d, ydc=%d, zdc+2=%d, Reb=%f, tauc=%f, G=%.10f, Ub=%f, CFLnumber=%f\n",xdc,ydc,zdc,Reb,tauc,Fx,ubulk,CFL);
		fprintf(tm,"xdf=%d, ydf=%d, zdf+2=%d, Reb=%f, tauf=%f, G=%.10f, Ub=%f, CFLnumber=%f\n",xdf,ydf,zdf,Reb,tauf,Fx,ubulk,CFL);
		fprintf(tm,"NZ fine +2 = %d   GR=%d \n",NZs,GR);
		fprintf(tm,"number of processors: %d\n",pd.numproc);
		fprintf(tm,"Executed File's Name: %s\n",__FILE__);
		fprintf(tm,"Compilation Date: %s\n",__DATE__);
		fprintf(tm,"Compilation Time: %s\n",__TIME__);
		fprintf(tm,"Compiled wit standard C Compiler? (1 means yes): %d\n",__STDC__);
		fprintf(tm,"Compilation Option: xcut\n");
#ifdef RSLIP
		fprintf(tm,"slip option: RSLIP\n");
#endif
#ifdef PSLIP
		fprintf(tm,"slip option: PSLIP\n");
#endif
#ifdef SSLIP
		fprintf(tm,"slip option: SSLIP\n");
#endif
		fprintf(tm,"Total computational time is: %f s\n",t2-t4);
		fprintf(tm,"Unit non-dimensional time is: %f iterations\n",utime);
		fprintf(tm,"Starting time: %d\n",TSTR);
		fprintf(tm,"Current Time is: %d\n",TEND);
		fprintf(tm,"size of lbf-c: %lu lbf-f: %lu\n",(xdc-1)*XDIM*(ydc-1)*YDIM*zdc*19*2*sizeof(double),(xdf-1)*XDIM*(ydf-1)*YDIM*zdf*19*2*sizeof(double));
		fprintf(tm,"size of vel-c: %lu vel-f: %lu\n",(xdc-1)*XDIM*(ydc-1)*YDIM*zdc*3*sizeof(double),(xdf-1)*XDIM*(ydf-1)*YDIM*zdf*3*sizeof(double));
		fprintf(tm,"size of den-c: %lu den-f: %lu\n",(xdc-1)*XDIM*(ydc-1)*YDIM*zdc*sizeof(double),(xdf-1)*XDIM*(ydf-1)*YDIM*zdf*sizeof(double));
		fprintf(tm,"size of double: %lu\n",sizeof(double));
		fprintf(tm,"size of uin %lu\n",23*3*(xdc-1)*XDIM*(ydc-1)*YDIM*sizeof(double));
		fclose(tm);
		if (pd.numproc==1)
			V=output(xdc, ydc, zdc, xdf, ydf, zdf, gratio, V,Fx);
	}
	V=statistics(xdf, ydf, zdf, xdc, ydc, zdc, V, pd, (TM+(ss/(int)utime)), GR);
	
	V.rho[0]=0.;rho=0.;rhozero=0.;
	for (i=0;i<(xdf-1);i++)
	{
		for (j=0;j<(ydf-1);j++)
		{
			for (k=0;k<zdf;k++)
			{
				rhozero=0.;
				for (a=0;a<19;a++)
					rhozero += Fb(V.slftemp,i,j,k,a,xdf,ydf,zdf,19);	
				if (k==(zdf-1))
					coef = 0.5;
				else
					coef = 1.;
				rho+= (coef*rhozero);
			}
		}
	}
	for (i=0;i<(xdc-1);i++)
	{
		for (j=0;j<(ydc-1);j++)
		{
			for (k=1;k<(zdc-1);k++)
			{
				rhozero=0.;
				for (a=0;a<19;a++)
					rhozero += Fb(V.ftemp,i,j,k,a,xdc,ydc,zdc,19);	
				if ((k==1)||(k==(zdc-2)))
					coef=0.5*GR*GR*GR;
				else
					coef=1.*GR*GR*GR;
				rho+=(coef*rhozero);
			}
		}
	}
	for (i=0;i<(xdf-1);i++)
	{
		for (j=0;j<(ydf-1);j++)
		{
			for (k=0;k<zdf;k++)
			{
				rhozero=0.;
				for (a=0;a<19;a++)
					rhozero += Fb(V.suftemp,i,j,k,a,xdf,ydf,zdf,19);
				if (k==0)
					coef = 0.5;
				else
					coef = 1.;
				rho+= (coef*rhozero);
			}
		}
	}
	V.rho[0] = rho;
	printf("V.rho = %f xdf=%d ydf=%d zdf=%d NZs=%d\n",V.rho[0],xdf,ydf,zdf,NZs);
	rho /= ((double)((xdf-1.)*(ydf-1.)*XDIM*YDIM*(NZs-2.)));
	MPI_Allreduce(&rho,&rhozero,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    	V.rho[0]=rhozero;
	printf("Average Density in the channel at end: %.12f  myrank=%d\n",V.rho[0],pd.myrank);
	
	free(V.f);
	free(V.ftemp);
	free(V.leftbufr);
	free(V.leftbufs);
	free(V.rightbufr);
	free(V.rightbufs);
	free(V.dwbufr);
	free(V.dwbufs);
	free(V.upbufr);
	free(V.upbufs);
	free(V.s);

 	free(V.suf);
	free(V.suftemp);
	free(V.suleftbufr);
	free(V.suleftbufs);
	free(V.surightbufr);
	free(V.surightbufs);
	free(V.sudwbufr);
	free(V.sudwbufs);
	free(V.suupbufr);
	free(V.suupbufs);
	free(V.sus);

 	free(V.slf);
	free(V.slftemp);
	free(V.slleftbufr);
	free(V.slleftbufs);
	free(V.slrightbufr);
	free(V.slrightbufs);
	free(V.sldwbufr);
	free(V.sldwbufs);
	free(V.slupbufr);
	free(V.slupbufs);
	free(V.sls);

	free(V.lin);
	free(V.linp);
	free(V.linpp);
	free(V.uin);
	free(V.uinp);
	free(V.uinpp);
	
	free(V.lint);
	free(V.uint);

	for (i=0;i<12;i++)
	  MPI_Request_free(&V.req[i]);

	for (i=0;i<12;i++)
	  MPI_Request_free(&V.sureq[i]);

	for (i=0;i<12;i++)
	  MPI_Request_free(&V.slreq[i]);
	
	MPI_Info_free(&info);
	
	MPI_Finalize();
	
	return(0);
}

