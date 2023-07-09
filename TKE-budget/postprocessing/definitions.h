// Modified postpc.c
// This should be used along with the postproc-budget
//This code gives all of the statistics and budgets

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define Reb 7200.0000
#define ubulk 0.0866034
#define TSTR 1130
#define TEND 1134
#define DT 1
#define TEXT "GR4"
typedef double DP;

#define xdv 512
#define ydv 256
#define zdv 197
/**************************************************************************
 * fine grid's parameters                                                 *
 * gratio : dx^c/dx^f ratio of the coarse to fine grid's resolution       *
 * ***********************************************************************/
#define gratio 4

// Number of z grid point on fine grid 
// a +1 is added inside the code to z for the wall points
#define zds 56