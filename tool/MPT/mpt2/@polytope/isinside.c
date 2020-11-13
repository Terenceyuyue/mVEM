/*  isinside.c
  
  MEX version of isinside()
  
  Usage:
   [isin, inwhich, closest] = isinside(Pn, x0, Options)

  See "help polytope/isinside" for a complete description.

*/

/* Copyright (C) 2006 by Michal Kvasnica (kvasnica@control.ee.ethz.ch) */

/*  This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/ 

#include "mex.h"

typedef struct options_s {
    double abs_tol;
    int fastbreak;
} options_t;

typedef struct polytope_s {
    mxArray *H;
    mxArray *K;
} polytope_t;

int isinside(const polytope_t *P, const double *x0, const int dimx0, const double abs_tol)
{
    int ic, ir, nr, nc;
    double *H, *K, rowsum;
    
    H = mxGetPr(P->H);
    K = mxGetPr(P->K);
    nr = mxGetM(P->K);
    nc = mxGetN(P->H);

    for (ir=0; ir<nr; ir++)
    {
        rowsum = 0;
        /* compute H*x0 */
        for (ic=0; ic<nc; ic++)
            rowsum = rowsum + H[ir + ic*nr]*x0[ic];

        /* exit immediately if x0 violates a hyperplane */
        if ((rowsum - K[ir]) > abs_tol)
            return 0;
    }
    
    /* all hyperplanes satisfied, we are inside */
    return 1;
}


void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
    mxArray *tmp, *Array, *mptOptions;
    double *tmpr, *x0, *inwhichp, *outinwhichp, *rcheb;
    int dimx0, nc, isin = 0;
    long i, nPn, true_nPn, ninside = 0, cnt;
    options_t options;
    polytope_t **Pn;
    
    if ((nrhs < 2) || (nrhs > 3))
        mexErrMsgTxt("Wrong number of input arguments.");
    
    if (nlhs > 3)
        mexErrMsgTxt("Too many output arguments.");
    
    if ((nrhs == 3) && !mxIsStruct(prhs[2]))
        mexErrMsgTxt("Third input must be a structure.");

    if (strcmp(mxGetClassName(prhs[0]), "polytope"))
        mexErrMsgTxt("First input must be a polytope object.");
    
    /* check dimension of x0 */
    if (mxGetN(prhs[1]) > 1)
        mexErrMsgTxt("Wrong dimension of x0.");
    
    /* set default options, read abs_tol from mptOptions.abs_tol */
    if (!(mptOptions = mexGetVariablePtr("global", "mptOptions")))
        mexErrMsgTxt("You must call mpt_init first.");
    if (!(tmp = mxGetField(mptOptions, 0, "abs_tol")))
        mexErrMsgTxt("You must call mpt_init first.");
    options.abs_tol = *mxGetPr(tmp);
    options.fastbreak = 0;
    
    /* read options */
    if (nrhs == 3) 
    {
        if (tmp = mxGetField(prhs[2], 0, "abs_tol"))
        {
            tmpr = mxGetPr(tmp);
            options.abs_tol = tmpr[0];
        }
        
        if (tmp = mxGetField(prhs[2], 0, "fastbreak"))
        {
            if (mxIsLogical(tmp))
                tmpr = mxGetLogicals(tmp);
            else
                tmpr = mxGetPr(tmp);
            options.fastbreak = (int)tmpr[0];
        }
    }
    
    Array = mxGetField(prhs[0], 0, "Array");
    
    nPn = mxGetM(Array)*mxGetN(Array);
    if (nPn==0)
    {
        /* input is a single polytope */
        polytope_t *P;
        Pn = mxMalloc(sizeof(polytope_t*));
        P = mxMalloc(sizeof(polytope_t));
        P->H = mxGetField(prhs[0], 0, "H");
        P->K = mxGetField(prhs[0], 0, "K");
        Pn[0] = P;
        
    } else {
        /* input is a polytope array */
        Pn = mxMalloc(nPn*sizeof(polytope_t*));
        for (i=0; i<nPn; i++)
        {
            polytope_t *P;
            P = mxMalloc(sizeof(polytope_t));
            P->H = mxGetField(mxGetCell(Array, i), 0, "H");
            P->K = mxGetField(mxGetCell(Array, i), 0, "K");
            Pn[i] = P;
        }
    }
    
    dimx0 = mxGetM(prhs[1]);
    x0 = (double *)mxGetPr(prhs[1]);

    /* check dimension */
    nc = mxGetN(Pn[0]->H);
    if (nc != dimx0)
    {
        /* if P is an empty polytope, return {0, [], []} */
        rcheb = (double *)mxGetPr(mxGetField(prhs[0], 0, "RCheb"));
        if (*rcheb < -1) {
            plhs[0] = mxCreateLogicalScalar(0);
            plhs[1] = mxCreateDoubleMatrix(0, 0, mxREAL);
            plhs[2] = mxCreateDoubleScalar(0);
            mxFree(Pn);
            return;
        } else {
            mexErrMsgTxt("Wrong dimension of x0.");
        }
    }
    
    /* allocate output vector */
    true_nPn = nPn>0 ? nPn : 1;
    inwhichp = mxGetPr(mxCreateDoubleMatrix(true_nPn, 1, mxREAL));
    
    /* we can always do fastbreak if the user doesn't ask for "inwhich" */
    if (nlhs < 2)
        options.fastbreak = 1;
    
    for (i=0; i<true_nPn; i++)
    {
        if (isinside(Pn[i], x0, dimx0, options.abs_tol))
        {
            if (options.fastbreak)
            {
                /* at least one polytope found, exit quickly */
                plhs[0] = mxCreateLogicalScalar(1);
                plhs[1] = mxCreateDoubleScalar(i+1);
                plhs[2] = mxCreateDoubleMatrix(0, 0, mxREAL);
                mxFree(Pn);
                return;
            }
            inwhichp[i] = 1;
            isin = 1;
            ninside++;
        }
    }
    
    /* create the inwhich array which contains indicies of polytopes 
     which contain x0 */
    plhs[1] = mxCreateDoubleMatrix(ninside, 1, mxREAL);
    if (ninside > 0) 
    {
        outinwhichp = mxGetPr(plhs[1]);
        cnt = 0;
        for (i=0; i<true_nPn; i++)
            if (inwhichp[i])
                outinwhichp[cnt++] = i+1;
    }
    
    plhs[0] = mxCreateLogicalScalar(isin);
    if (!isin &  nlhs>2)
    {
        /* calcuate index of a region which is closest to x0. 
           But only do that if the user requests such information (i.e. nlhs=3) 
        */
        if (true_nPn == 1) {
            plhs[2] = mxCreateDoubleScalar(1);
        } else {
            plhs[2] = mxCreateDoubleScalar(0);
            mexCallMATLAB(1, &plhs[2], 2, prhs, "closestRegion");
        }
    } else {
        plhs[2] = mxCreateDoubleMatrix(0, 0, mxREAL);
    }
    mxFree(Pn);
    return;
}
