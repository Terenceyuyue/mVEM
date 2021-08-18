/*
  Matlab interface to fourier elimination algorithm
  
  Author: Colin Jones
  Email : cnj22@cam.ac.uk
  Date  : April 17, 2004
*/

#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "fourier.h"

int DIM;
CSet *vset;
CSet *hset;

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  mxArray *P;
  int      lnP;
  double  *H;
  int      n,m;
  int     *ax;
  double  *pdtmp;
  int      i,j;
  double   tol, qtol;

  /*  Check for proper number of arguments. */
  if(nrhs < 2) mexErrMsgTxt("Two arguments required");

  if(nrhs > 2) tol = mxGetScalar(prhs[2]);
  else         tol = 1e-6;
	
  // Convert quasi-tolerance to radians
  if(nrhs > 3) qtol = mxGetScalar(prhs[3])*0.0174532925;
  else         qtol = 0.01*0.0174532925;

  /* Get the arguments */
  H = mxGetPr(prhs[0]);
  n = mxGetN(prhs[0]);
  m = mxGetM(prhs[0]);
  DIM = n-1;

  vset = new CSet(DIM);
  hset = new CSet(m);

  ax  = vset->newset(); // axes to project onto

  pdtmp = mxGetPr(prhs[1]);
  for(i=0;i<mxGetN(prhs[1]);i++)
	ax[int(pdtmp[i])-1] = 1;
  /*
  printf("Rows: %i, Columns: %i\n",m,n);
  
  printf("H = \n");
  for(i=0;i<m;i++)
  {
    printf("[");
  	for(j=0;j<n;j++)
	  //	  printf("%6.2f ",H[j + m*i]);
	  printf("%6.2f ",H[i + m*j]);
  	printf("]\n");
  };
  printf("\n\n");

  printf("Projecting onto axes: "); vset->printset(ax); printf("\n");
  */
  P = fourier(H, m, ax, tol, qtol);

  if(nlhs == 1)
	plhs[0] = P;
  
  vset->freeset(ax);

  delete vset;
  delete hset;
}
