#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

/*
Fast set difference for integer arrays

%  2004/04/01
%     Colin Jones, Cambridge Control Laboratory, Cambridge UK
%     cnj22@cam.ac.uk
*/

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
  /* cnjsetdiff
  /* 
  /* input: A,B two integer arrays
  /* output: setdiff(A,B)
  /*
  /* Assumes inputs are sorted!
  */

	int    nA,nB,nC;
  double *A,*B,*C;
  int i,j,k,foundB,found;


  if(nrhs != 2 || nlhs != 1) mexErrMsgTxt("C = cnjsetdiff(A,B) Inputs must be pre-sorted!");

  /*if(mxIsDouble(prhs[0]) || mxIsDouble(prhs[1]))
	mexErrMsgTxt("Input must be integer");*/

  nA = mxGetN(prhs[0]);
  if(nA == 1) nA = mxGetM(prhs[0]);
  nB = mxGetN(prhs[1]);
  if(nB == 1) nB = mxGetM(prhs[1]);

  A = mxGetPr(prhs[0]);
  B = mxGetPr(prhs[1]);

  plhs[0] = mxCreateDoubleMatrix(nA,1,mxREAL);
  C = mxGetPr(plhs[0]);

  nC = 0;
  foundB = 0;
  for(i=0;i<nA;i++)
  {
    found = 0;
    j = foundB;
    while(j < nB)
    {
      if(A[i] == B[j])
      {
        found = 1;
        foundB = j+1;
        break;
      };
      j++;
    };
    if(found == 0)
      C[nC++] = A[i];
  };

  mxSetM(plhs[0],nC);
}
