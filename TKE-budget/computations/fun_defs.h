/************************************************************************
 * By Amirreza Rastegari                                                *
 *   arstgri@gmail.com                                                  *
 *                                                                      *
 * Assumes a Multigrid data input                                       *
 * calculates rms, vorticity and mean velocities                        *
 * Calculates the energy budget2D                                       *
 * Skewness and Kurtosis                                                *
 * Assumes a 2D mean velocity i.e. U(y,z)                               *
 ************************************************************************/

POINTER reading2(int xdc, int ydc, int zdc, int xdf, int ydf, int zdf, POINTER V, PDATA pd, int ts);
POINTER reading(int xdc, int ydc, int zdc, int xdf, int ydf, int zdf, POINTER V, PDATA pd, long ts, DP ubulk);
POINTER uv_fluc(int xdf, int ydf, int zdf, POINTER V, PDATA pd, long ts, int position);
POINTER budget2D(int xdc, int ydc, int zdc, int xdf, int ydf, int zdf, int GR, POINTER V, PDATA pd, long ts);
