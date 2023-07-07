#ifndef MYHEADER_H
#define MYHEADER_H

POINTER init(int xdc, int ydc, int zdc, int xdf, int ydf, int zdf, POINTER V);
POINTER solid_init(int xdc, int ydc, int zdc, int xdf, int ydf, int zdf, POINTER V);
POINTER slip_ridges_init(int xdf, int ydf, int zdf, POINTER V, PDATA pd);
POINTER slip_spanwise_init(int xdf, int ydf, int zdf, POINTER V, PDATA pd);
POINTER slip_posts_init(int xdf, int ydf, int zdf, POINTER V, PDATA pd);
POINTER encalc(int xdc, int ydc , int zdc, int xdf, int ydf , int zdf, int GR, POINTER V, PDATA pd, int TM, long cntr,DP ubulk);
POINTER encalc2(int xdc, int ydc , int zdc, int xdf, int ydf , int zdf, int GR, POINTER V, PDATA pd, int TM, long cntr,DP ubulk);
POINTER fine_force(int xdf, int ydf , int zdf, int zdc, int GR, POINTER V, PDATA pd, DP ubulk);
POINTER tinterpolate_qd(int xdc, int ydc, int zdc, POINTER V, int iter, int GR);
POINTER tinterpolate_hrmt(int xdc, int ydc, int zdc, POINTER V, int iter, int GR);
POINTER tinterpolate_ln(int xdc, int ydc, int zdc, POINTER V, int iter, int GR);
void derivative(int xd, int yd, double *addr);
POINTER xcomputations(int xdc, int ydc, int zdc, int GR, POINTER V,PDATA pd, DP Gx, DP tauc, DP ubulk);
POINTER xfcomputations(int xdf, int ydf, int zdf, double GR, POINTER V, PDATA pd, DP tauf, DP ubulk);
void bcucof(int xdc, int ydc, int ic, int jc, int a, int GR, double *addr, double *c);
void bcucof_opt(int xdc, int ydc, int ic, int jc, int GR, double *addr, double *c);
double bcuint(int xdc, int ydc, int GR, int i, int j, double *addr, int a);
POINTER ctfdtransfer_opt(int xdc, int ydc, int zdc, int xdf, int ydf, int zdf, int itr, double tauc, double tauf, int GR, POINTER V, PDATA pd);
POINTER ctfdtransfer(int xdc, int ydc, int zdc, int xdf, int ydf, int zdf, int itr, double tauc, double tauf, int GR, POINTER V, PDATA pd);
POINTER ftcdtransfer(int xdc, int ydc, int zdc, int xdf, int ydf, int zdf, int itr, double tauc, double tauf, int GR, POINTER V, PDATA pd);
POINTER ftcdtransfer_flt(int xdf, int ydf, int zdf, int xdc, int ydc, int zdc, int itr, double tauc, double tauf, int GR, POINTER V, PDATA pd);
POINTER ftcdtransfer_lesfltr(int xdf, int ydf, int zdf, int xdc, int ydc, int zdc, int itr, double tauc, double tauf, int GR, POINTER V, PDATA pd);
POINTER exchange(POINTER V);
POINTER fexchange(POINTER V);
POINTER output_twoD(int xdf, int ydf, int zdf, int xdc, int ydc, int zdc, int GR, POINTER V, PDATA pd);
POINTER output_twoD2(int xdf, int ydf, int zdf, int xdc, int ydc, int zdc, int GR, POINTER V, PDATA pd);
POINTER output(int xdc, int ydc, int zdc, int xdf, int ydf, int zdf, int GR, POINTER V, DP Fx);
POINTER wshearstr(int xdf, int ydf, int zdf, int zdc, int GR, POINTER V, PDATA pd, DP tauf, int TM, long cntr);
POINTER statistics(int xdf, int ydf, int zdf, int xdc, int ydc, int zdc, POINTER V, PDATA pd, long ts, int GR);
POINTER vel_print(int xdf, int ydf, int zdf, int xdc, int ydc, int zdc, POINTER V, long cntr);
void sp_interpolate_qd(int *ind, int *jnd, int *knd, int xd, int yd, int zd, double dx, double dy, double dz, double *adr, double *ans);
POINTER initializer(int xdf, int ydf, int zdf, double tauf, int xdc, int ydc, int zdc, double tauc, int xd, int yd, int zd, double tau, int GR, POINTER V, PDATA pd, int cntr);
POINTER inter_init(int xdc, int ydc, int zdc, POINTER V, PDATA pd);

#endif /* MYHEADER_H */