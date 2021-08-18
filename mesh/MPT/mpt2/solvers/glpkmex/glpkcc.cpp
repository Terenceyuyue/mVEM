/*

Copyright (C) 2001-2007 Nicolo' Giorgetti.

This file is part of GLPKMEX.

GLPKMEX is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2, or (at your option) any
later version.

This part of code is distributed with the FURTHER condition that it 
can be compiled and linked with the Matlab libraries and it can be 
used within the Matlab environment.

GLPKMEX is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with Octave; see the file COPYING.  If not, write to the Free
Software Foundation, 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

*/

#include <cfloat>
#include <csetjmp>
#include <ctime>
#include <cstring>

// From Matlab
#include "mex.h"

extern "C" {
  #include <glpk.h>
  #include <glplib.h>
} 

#define NIntP 18
#define NRealP 10

int glpIntParam[NIntP] = {
  0,
  1,
  0,
  1,
  0,
  -1,
  0,
  200,
  1,
  2,
  0,
  1,
  0,
  0,
  2,
  2,
  1,
  0
};

int IParam[NIntP] = {
  LPX_K_MSGLEV,
  LPX_K_SCALE,
  LPX_K_DUAL,
  LPX_K_PRICE,
  LPX_K_ROUND,
  LPX_K_ITLIM,
  LPX_K_ITCNT,
  LPX_K_OUTFRQ,
  LPX_K_MPSINFO,
  LPX_K_MPSOBJ,
  LPX_K_MPSORIG,
  LPX_K_MPSWIDE,
  LPX_K_MPSFREE,
  LPX_K_MPSSKIP,
  LPX_K_BRANCH,
  LPX_K_BTRACK,
  LPX_K_PRESOL,
  LPX_K_USECUTS
};


double glpRealParam[NRealP] = {
  0.07,
  1e-7,
  1e-7,
  1e-9,
  -DBL_MAX,
  DBL_MAX,
  -1,
  0.0,
  1e-6,
  1e-7
};

int RParam[NRealP] = {
  LPX_K_RELAX,
  LPX_K_TOLBND,
  LPX_K_TOLDJ,
  LPX_K_TOLPIV,
  LPX_K_OBJLL,
  LPX_K_OBJUL,
  LPX_K_TMLIM,
  LPX_K_OUTDLY,
  LPX_K_TOLINT,
  LPX_K_TOLOBJ
};

static jmp_buf mark;  /*-- Address for long jump */

/* -- Obsolete ---
int glpk_fault_hook (void *info, char *msg)
{
  char errmsg[1024];
  sprintf(errmsg,"CRITICAL ERROR in GLPK!\n%s\n",msg);
  mexErrMsgTxt(errmsg);
  longjmp( mark, -1 );
} */

static int glpk_print_hook (void *info, const char *msg)
{
  mexPrintf("%s",msg);
  return 1;
}

int glpk (int sense, int n, int m, double *c, int nz, int *rn, int *cn,
      	 double *a, double *b, char *ctype, int *freeLB, double *lb,
      	 int *freeUB, double *ub, int *vartype, int isMIP, int lpsolver,
      	 int save_pb, char *save_filename, char *filetype, 
         double *xmin, double *fmin, double *status,
      	 double *lambda, double *redcosts, double *time, double *mem)
{
  int typx = 0;
  int method;

  clock_t t_start = clock();

  // Obsolete
  //lib_set_fault_hook (NULL, glpk_fault_hook);

  //Redirect standard output
  if (glpIntParam[0] > 1) glp_term_hook (glpk_print_hook, NULL);
  else glp_term_hook (NULL, NULL);

  //-- Create an empty LP/MILP object
  glp_prob *lp = glp_create_prob ();

  //-- Set the sense of optimization
  if (sense == 1)
    glp_set_obj_dir (lp, GLP_MIN);
  else
    glp_set_obj_dir (lp, GLP_MAX);

  //-- Define the number of unknowns and their domains.
  glp_add_cols (lp, n);
  for (int i = 0; i < n; i++)
  {
    //-- Define type of the structural variables
    if (! freeLB[i] && ! freeUB[i])
      glp_set_col_bnds (lp, i+1, GLP_DB, lb[i], ub[i]);
    else
	  {
      if (! freeLB[i] && freeUB[i])
        glp_set_col_bnds (lp, i+1, GLP_LO, lb[i], ub[i]);
      else
      {
        if (freeLB[i] && ! freeUB[i])
		      glp_set_col_bnds (lp, i+1, GLP_UP, lb[i], ub[i]);
	      else
		      glp_set_col_bnds (lp, i+1, GLP_FR, lb[i], ub[i]);
	    }
	  }
  
  // -- Set the objective coefficient of the corresponding
  // -- structural variable. No constant term is assumed.
  glp_set_obj_coef(lp,i+1,c[i]);

  if (isMIP)
    glp_set_col_kind (lp, i+1, vartype[i]);
  }

  glp_add_rows (lp, m);

  for (int i = 0; i < m; i++)
  {
    /*  If the i-th row has no lower bound (types F,U), the
        corrispondent parameter will be ignored.
        If the i-th row has no upper bound (types F,L), the corrispondent
        parameter will be ignored.
        If the i-th row is of S type, the i-th LB is used, but
        the i-th UB is ignored.
    */

    switch (ctype[i])
    {
      case 'F': typx = GLP_FR; break;
      // upper bound
	    case 'U': typx = GLP_UP; break;
      // lower bound
	    case 'L': typx = GLP_LO; break;
      // fixed constraint
	    case 'S': typx = GLP_FX; break;
      // double-bounded variable
      case 'D': typx = GLP_DB; break;
	  }
      
    glp_set_row_bnds (lp, i+1, typx, b[i], b[i]);

  }
  // Load constraint matrix A
  glp_load_matrix (lp, nz, rn, cn, a);

  // Save problem
  if (save_pb) {
    if (!strcmp(filetype,"cplex")){
      if (lpx_write_cpxlp (lp, save_filename) != 0) {
	        mexErrMsgTxt("glpkcc: unable to write the problem");
	        longjmp (mark, -1);
      }
    }else{
      if (!strcmp(filetype,"fixedmps")){
        if (lpx_write_mps (lp, save_filename) != 0) {
          mexErrMsgTxt("glpkcc: unable to write the problem");
	        longjmp (mark, -1);  
        }
      }else{
        if (!strcmp(filetype,"freemps")){
          if (lpx_write_freemps (lp, save_filename) != 0) {
            mexErrMsgTxt("glpkcc: unable to write the problem");
	          longjmp (mark, -1);
          }
        }else{// plain text
          if (lpx_print_prob (lp, save_filename) != 0) {
            mexErrMsgTxt("glpkcc: unable to write the problem");
	          longjmp (mark, -1);
          } 
        } 
      }    
    } 
//    mexErrMsgTxt("problem saved");
  }
  //-- scale the problem data (if required)
  if (glpIntParam[1] && (! glpIntParam[16] || lpsolver != 1))
    lpx_scale_prob (lp);

  //-- build advanced initial basis (if required)
  if (lpsolver == 1 && ! glpIntParam[16])
    lpx_adv_basis (lp);

  glp_smcp sParam;
  glp_init_smcp(&sParam);
  
  //-- set control parameters
  if (lpsolver==1){
    //remap of control parameters for simplex method
    sParam.msg_lev=glpIntParam[0];	// message level
    // simplex method: primal/dual
    if (glpIntParam[2]==0) sParam.meth=GLP_PRIMAL;		
    else sParam.meth=GLP_DUALP;
    // pricing technique
    if (glpIntParam[3]==0) sParam.pricing=GLP_PT_STD;
    else sParam.pricing=GLP_PT_PSE;
    //sParam.r_test not available
    sParam.tol_bnd=glpRealParam[1];	// primal feasible tollerance
    sParam.tol_dj=glpRealParam[2];	// dual feasible tollerance
    sParam.tol_piv=glpRealParam[3];	// pivot tollerance
    sParam.obj_ll=glpRealParam[4];	// lower limit
    sParam.obj_ul=glpRealParam[5];	// upper limit
    // iteration limit
    if (glpIntParam[5]==-1) sParam.it_lim=INT_MAX;
    else sParam.it_lim=glpIntParam[5];   
    // time limit
    if (glpRealParam[6]==-1) sParam.tm_lim=INT_MAX;
    else sParam.tm_lim=(int) glpRealParam[6];	
    sParam.out_frq=glpIntParam[7];	// output frequency
    sParam.out_dly=(int) glpRealParam[7];	// output delay
    // presolver
    if (glpIntParam[16]) sParam.presolve=GLP_ON;
    else sParam.presolve=GLP_OFF;
  }else{
	for(int i = 0; i < NIntP; i++)
		lpx_set_int_parm (lp, IParam[i], glpIntParam[i]);
		
	for (int i = 0; i < NRealP; i++)
		lpx_set_real_parm (lp, RParam[i], glpRealParam[i]);
  }
  

  // Choose simplex method ('S') or interior point method ('T') to solve the problem
  if (lpsolver == 1)
    method = 'S';
  else
    method = 'T';
	
  int errnum;

  switch (method){
    case 'S': {
      if (isMIP){
	    method = 'I';
	    errnum = lpx_intopt (lp);
      }
      else{
		errnum = glp_simplex(lp, &sParam);
		errnum += 100; //this is to avoid ambiguity in the return codes.
	  }
    }
    break;

    case 'T': errnum = lpx_interior(lp); break;

    default:  xassert (method != method);
  }

  /*  errnum assumes the following results:
      errnum = 0 <=> No errors
      errnum = 1 <=> Iteration limit exceeded.
      errnum = 2 <=> Numerical problems with basis matrix.
  */
  if (errnum == LPX_E_OK || errnum==100){
    // Get status and object value
    if (isMIP)
    {
      *status = glp_mip_status (lp);
      *fmin = glp_mip_obj_val (lp);
    }
    else
    {
      if (lpsolver == 1)
      {
        *status = glp_get_status (lp);
        *fmin = glp_get_obj_val (lp);
	    }
      else
      {
        *status = glp_ipt_status (lp);
        *fmin = glp_ipt_obj_val (lp);
	    }
    }
    // Get optimal solution (if exists)
    if (isMIP)
    {
      for (int i = 0; i < n; i++)
        xmin[i] = glp_mip_col_val (lp, i+1);
    }
    else
    {
      /* Primal values */
      for (int i = 0; i < n; i++)
      {
        if (lpsolver == 1)
          xmin[i] = glp_get_col_prim (lp, i+1);
        else
		      xmin[i] = glp_ipt_col_prim (lp, i+1);
      }
      /* Dual values */
      for (int i = 0; i < m; i++)
      {
        if (lpsolver == 1) lambda[i] = glp_get_row_dual (lp, i+1);
	     else lambda[i] = glp_ipt_row_dual (lp, i+1);
      }
      /* Reduced costs */
      for (int i = 0; i < glp_get_num_cols (lp); i++)
      {
        if (lpsolver == 1) redcosts[i] = glp_get_col_dual (lp, i+1);
        else redcosts[i] = glp_ipt_col_dual (lp, i+1);
      }
    }

    *time = (clock () - t_start) / CLOCKS_PER_SEC;
    
/*
   	glp_long tpeak;
    lib_mem_usage(NULL, NULL, NULL, &tpeak);
    *mem=(double)(4294967296.0 * tpeak.hi + tpeak.lo) / (1024);
       */
    *time = (clock () - t_start) / CLOCKS_PER_SEC;

	  glp_delete_prob (lp);
    return 0;
  }

  glp_delete_prob (lp);

  *status = errnum;

  return errnum;
}

#define GLPK_GET_REAL_PARAM(PAR, NAME, IDX) \
  do \
    { \
    	mxArray *mxtmp=mxGetField(PAR,0,NAME); \
      if ( mxtmp != NULL) \
      { \
			double *rdtmp=mxGetPr(mxtmp); \
			glpRealParam[IDX] = *rdtmp; \
      } \
	 } \
    while(0)

#define GLPK_GET_INT_PARAM(PAR, NAME, VAL) \
  do \
    { \
    	mxArray *mxtmp=mxGetField(PAR,0,NAME); \
      if ( mxtmp != NULL) \
	   { \
		  double *rdtmp=mxGetPr(mxtmp); \
        \
        VAL =(int) *rdtmp; \
     } \
	} \
   while (0)

    
//-- Input arguments
#define	C_IN	     prhs[0]
#define	A_IN	     prhs[1]
#define	B_IN	     prhs[2]
#define LB_IN	     prhs[3]
#define UB_IN      prhs[4]
#define CTYPE_IN   prhs[5]
#define VARTYPE_IN prhs[6]
#define	SENSE_IN   prhs[7]
#define PARAM      prhs[8]

//-- Output Arguments
#define	 XMIN_OUT     plhs[0]
#define	 FMIN_OUT     plhs[1]
#define	 STATUS_OUT   plhs[2]
#define  EXTRA_OUT    plhs[3]


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  
  if(strcmp(glp_version(),"4.18")<0){
	mexErrMsgTxt("This MEX interface is compatible only with GLPK version 4.18 or higher.");
	}
  
  if (nrhs != 9){
      mexPrintf("MEX interface to GLPK Version %s\n",glp_version());
      mexPrintf("Internal interface for the GNU GLPK library.\n");
      mexPrintf("You should use the 'glpk' function instead.\n\n");
      mexPrintf("SYNTAX: [xopt, fmin, status, extra] = glpkcc(c, a, b, lb, ub, ctype, vartype, sense, param)\n");
      return;
  }

  //-- 1nd Input. A column array containing the objective function
  //--            coefficients.
  int mrowsc = mxGetM(C_IN);

  double *c=mxGetPr(C_IN);
  if (c == NULL) mexErrMsgTxt("glpkcc: invalid value of C");
  
  
  //-- 2nd Input. A matrix containing the constraints coefficients.
  // If matrix A is NOT a sparse matrix
  double *A = mxGetPr(A_IN); // get the matrix
  if(A==NULL) mexErrMsgTxt("glpkcc: invalid value of A");
  
  int mrowsA = mxGetM(A_IN);
  
  int *rn;
  int *cn;
  double *a;
  int nz = 0;
  
  if(!mxIsSparse(A_IN)){
     rn=(int *)mxCalloc(mrowsA*mrowsc+1,sizeof(int));
     cn=(int *)mxCalloc(mrowsA*mrowsc+1,sizeof(int));
	   a=(double *)mxCalloc(mrowsA*mrowsc+1,sizeof(double));

     for (int i = 0; i < mrowsA; i++){
      for (int j = 0; j < mrowsc; j++){
	     if (A[i+j*mrowsA] != 0){
	      nz++;
	      rn[nz] = i + 1;
	      cn[nz] = j + 1;
	      a[nz] = A[i+j*mrowsA];
	    }
	   }
     }
  }else{
	    /* NOTE: nnz is the actual number of nonzeros and is stored as the
          last element of the jc array where the size of the jc array is the
          number of columns + 1 */
	    nz = *(mxGetJc(A_IN) + mrowsc);
	    int *jc = mxGetJc(A_IN);
	    int *ir = mxGetIr(A_IN);
	    double *pr = mxGetPr(A_IN);

      rn=(int *)mxCalloc(nz+1,sizeof(int));
	    cn=(int *)mxCalloc(nz+1,sizeof(int));
	    a=(double *)mxCalloc(nz+1,sizeof(double));

      int nelc,count,row;
      count=0; row=0;
	    for(int i=1;i<=mrowsc;i++){
	      nelc=jc[i]-jc[i-1];
	      for(int j=0;j<nelc;j++){
		      count++;
		      rn[count]=ir[row]+1;
		      cn[count]=i;
		      a[count]=pr[row];
		      row++;
	      }
	    }
  }

  //-- 3rd Input. A column array containing the right-hand side value
  //	           for each constraint in the constraint matrix.
  double *b = mxGetPr(B_IN);
  
  if (b==NULL) mexErrMsgTxt("glpkcc: invalid value of b");
 

  //-- 4th Input. An array of length mrowsc containing the lower
  //--            bound on each of the variables.
  double *lb = mxGetPr(LB_IN);
  
  if (lb==NULL) mexErrMsgTxt("glpkcc: invalid value of lb");
      

  //-- LB argument, default: Free
  int *freeLB=(int *)mxCalloc(mrowsc,sizeof(int));
  for (int i = 0; i < mrowsc; i++) {
    if (lb[i]==-mxGetInf()){
      freeLB[i] = 1;
	 	}else freeLB[i] = 0;
  }

  //-- 5th Input. An array of at least length numcols containing the upper
  //--            bound on each of the variables.
  double *ub = mxGetPr(UB_IN);

  if (ub==NULL) mexErrMsgTxt("glpkcc: invalid value of ub");
      
  int *freeUB=(int *)mxCalloc(mrowsc,sizeof(int));
  for (int i = 0; i < mrowsc; i++)
  {
    if (ub[i]==mxGetInf())
		{
	     freeUB[i] = 1;
	  }else freeUB[i] = 0;
  }

  //-- 6th Input. A column array containing the sense of each constraint
  //--            in the constraint matrix.
  int size = mxGetNumberOfElements(CTYPE_IN) + 1;
  if (size==0) mexErrMsgTxt("glpkcc: invalid value of ctype");
  
  /* Allocate enough memory to hold the converted string. */
  char *ctype =(char *)mxCalloc(size, sizeof (char));

  /* Copy the string data from string_array_ptr and place it into buf. */
  if (mxGetString(CTYPE_IN, ctype, size) != 0)  mexErrMsgTxt("Could not convert string data.");
	  
  
  //-- 7th Input. A column array containing the types of the variables.
  size = mxGetNumberOfElements(VARTYPE_IN)+1;
  
  char *vtype = (char *)mxCalloc(size, sizeof (char));
  int *vartype = (int *)mxCalloc(size, sizeof (int));
  
  if (size==0) mexErrMsgTxt("glpkcc: invalid value of vartype");
    
  // Copy the string data from string_array_ptr and place it into buf.
  if (mxGetString(VARTYPE_IN, vtype, size) != 0)
	  mexErrMsgTxt("Could not convert string data.");
  
  int isMIP = 0;
  for (int i = 0; i < mrowsc ; i++)
  {
    switch (vtype[i]){
      case 'I': vartype[i] = GLP_IV; isMIP = 1; break;
      case 'B': vartype[i] = GLP_BV; isMIP = 1; break;
      default: vartype[i] = GLP_CV;   
    }
  }

  //-- 8th Input. Sense of optimization.
  int sense;
  
  double *tmp = mxGetPr(SENSE_IN);
  
  if (*tmp >= 0) sense = 1;
  else sense = -1;

  //-- 9th Input. A structure containing the control parameters.
  
  //-- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //-- Integer parameters
  //-- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //-- Level of messages output by the solver
  GLPK_GET_INT_PARAM (PARAM, "msglev", glpIntParam[0]);
  if (glpIntParam[0] < 0 || glpIntParam[0] > 3)
    {
      mexErrMsgTxt("glpkcc: param.msglev must be 0 (no output [default]) or 1 (error messages only) or 2 (normal output) or 3 (full output)");
    }
  
  //-- scaling option
  GLPK_GET_INT_PARAM (PARAM, "scale", glpIntParam[1]);
  if (glpIntParam[1] < 0 || glpIntParam[1] > 2)
    {
      mexErrMsgTxt("glpkcc: param.scale must be 0 (no scaling) or 1 (equilibration scaling [default]) or 2 (geometric mean scaling)");
    }

  //-- Dual dimplex option
  GLPK_GET_INT_PARAM (PARAM, "dual", glpIntParam[2]);
  if (glpIntParam[2] < 0 || glpIntParam[2] > 1)
    {
      mexErrMsgTxt("glpkcc: param.dual must be 0 (do NOT use dual simplex [default]) or 1 (use dual simplex)");
    }

  //-- Pricing option
  GLPK_GET_INT_PARAM (PARAM, "price", glpIntParam[3]);
  if (glpIntParam[3] < 0 || glpIntParam[3] > 1)
    {
      mexErrMsgTxt("glpkcc: param.price must be 0 (textbook pricing) or 1 (steepest edge pricing [default])");
    }

  //-- Solution rounding option
  GLPK_GET_INT_PARAM (PARAM, "round", glpIntParam[4]);
  if (glpIntParam[4] < 0 || glpIntParam[4] > 1)
    {
      mexErrMsgTxt("glpkcc: param.round must be 0 (report all primal and dual values [default]) or 1 (replace tiny primal and dual values by exact zero)");
    }

  //-- Simplex iterations limit
  GLPK_GET_INT_PARAM (PARAM, "itlim", glpIntParam[5]);

  //-- Simplex iterations count
  GLPK_GET_INT_PARAM (PARAM, "itcnt", glpIntParam[6]);

  //-- Output frequency, in iterations
  GLPK_GET_INT_PARAM (PARAM, "outfrq", glpIntParam[7]);

  //-- Branching heuristic option
  GLPK_GET_INT_PARAM (PARAM, "branch", glpIntParam[14]);
  if (glpIntParam[14] < 0 || glpIntParam[14] > 2)
    {
      mexErrMsgTxt("glpkcc: param.branch must be (MIP only) 0 (branch on first variable) or 1 (branch on last variable) or 2 (branch using a heuristic by Driebeck and Tomlin [default]");
    }

  //-- Backtracking heuristic option
  GLPK_GET_INT_PARAM (PARAM, "btrack", glpIntParam[15]);
  if (glpIntParam[15] < 0 || glpIntParam[15] > 2)
    {
      mexErrMsgTxt("glpkcc: param.btrack must be (MIP only) 0 (depth first search) or 1 (breadth first search) or 2 (backtrack using the best projection heuristic [default]");
    }

  //-- Presolver option
  GLPK_GET_INT_PARAM (PARAM, "presol", glpIntParam[16]);
  if (glpIntParam[16] < 0 || glpIntParam[16] > 1)
    {
      mexErrMsgTxt("glpkcc: param.presol must be 0 (do NOT use LP presolver) or 1 (use LP presolver [default])");
    }
  
  //-- Generating cuts
  GLPK_GET_INT_PARAM (PARAM, "usecuts", glpIntParam[17]);
  if (glpIntParam[17] < 0 || glpIntParam[17] > 1)
    {
      mexErrMsgTxt("glpkcc: param.usecuts must be 0 (do NOT generate cuts) or 1 (generate Gomory's cuts [default])");
    }


  //-- LPsolver option
  int lpsolver = 1;
  GLPK_GET_INT_PARAM (PARAM, "lpsolver", lpsolver);
  if (lpsolver < 1 || lpsolver > 2)
    {
      mexErrMsgTxt("glpkcc: param.lpsolver must be 1 (simplex method) or 2 (interior point method)");
    }

  //-- Save option
  int save_pb = 0;
  char *save_filename = NULL;
  char *filetype = NULL;
  GLPK_GET_INT_PARAM (PARAM, "save", save_pb);
  save_pb = (save_pb != 0);
  if (save_pb){   
    // -- Look for the name --
    mxArray *mxtmp=mxGetField(PARAM,0,"savefilename");
    if ( mxtmp != NULL ){
      int nl=mxGetNumberOfElements(mxtmp)+1;
      nl=nl+4; // increase size to consider then extension .xxx 
      save_filename=(char *)mxCalloc(nl,sizeof(char));
      if (mxGetString(mxtmp, save_filename, nl) != 0)
        mexErrMsgTxt("glpkcc: Could not load file name to save.");
    }else{
      // Default file name
      save_filename= (char *)mxCalloc(9, sizeof(char));
      strcpy(save_filename,"outpb");
    }
    
    // -- Look for the type --
    char save_filetype[4];
    mxArray *txtmp=mxGetField(PARAM,0,"savefiletype");
    if ( txtmp != NULL ){
      int nl=mxGetNumberOfElements(txtmp)+1; 
      filetype=(char *)mxCalloc(nl,sizeof(char));
      if (mxGetString(txtmp, filetype, nl) != 0)
        mexErrMsgTxt("glpkcc: Could not load file type.");
      if (!strcmp(filetype,"fixedmps") || !strcmp(filetype,"freemps")){
        strcpy(save_filetype,".mps");
      } else {
        if (!strcmp(filetype,"cplex")) strcpy(save_filetype,".lp");
        else {
          if (!strcmp(filetype,"plain")) strcpy(save_filetype,".txt");
        } 
      }  
    }else{
      filetype= (char *)mxCalloc(5, sizeof(char));
      strcpy(filetype,"cplex");
      strcpy(save_filetype,".lp"); // Default file type
    }  
    strcat(save_filename,save_filetype); // name.extension   
  }
  
  // MPS parameters
  //-- mpsinfo 
  GLPK_GET_INT_PARAM (PARAM, "mpsinfo", glpIntParam[8]);
  //-- mpsobj
  GLPK_GET_INT_PARAM (PARAM, "mpsobj", glpIntParam[9]);
  if (glpIntParam[9] < 0 || glpIntParam[9] > 2)
  {
    mexErrMsgTxt("glpkcc: param.mpsobj must be 0 (never output objective function row) or 1 (always output objective function row ) or 2 [default](output objective function row if the problem has no free rows)");
  }
  //-- mpsorig 
  GLPK_GET_INT_PARAM (PARAM, "mpsorig", glpIntParam[10]);
  //-- mpswide 
  GLPK_GET_INT_PARAM (PARAM, "mpswide", glpIntParam[11]);
  //-- mpsfree 
  GLPK_GET_INT_PARAM (PARAM, "mpsfree", glpIntParam[12]);
  
  

  //-- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //-- Real parameters
  //-- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //-- Ratio test option
  GLPK_GET_REAL_PARAM (PARAM, "relax", 0);

  //-- Relative tolerance used to check if the current basic solution
  //-- is primal feasible
  GLPK_GET_REAL_PARAM (PARAM, "tolbnd", 1);

  //-- Absolute tolerance used to check if the current basic solution
  //-- is dual feasible
  GLPK_GET_REAL_PARAM (PARAM, "toldj", 2);

  //-- Relative tolerance used to choose eligible pivotal elements of
  //--	the simplex table in the ratio test
  GLPK_GET_REAL_PARAM (PARAM, "tolpiv", 3);

  GLPK_GET_REAL_PARAM (PARAM, "objll", 4);

  GLPK_GET_REAL_PARAM (PARAM, "objul", 5);

  GLPK_GET_REAL_PARAM (PARAM, "tmlim", 6);

  GLPK_GET_REAL_PARAM (PARAM, "outdly", 7);

  GLPK_GET_REAL_PARAM (PARAM, "tolint", 8);

  GLPK_GET_REAL_PARAM (PARAM, "tolobj", 9);
  
  //-- Assign pointers to the output parameters
  const char **extranames=(const char **)mxCalloc(4,sizeof(*extranames));
  extranames[0]="lambda";
  extranames[1]="redcosts";
  extranames[2]="time";
  extranames[3]="memory";
  
  XMIN_OUT   = mxCreateDoubleMatrix(mrowsc, 1, mxREAL);
  FMIN_OUT   = mxCreateDoubleMatrix(1, 1, mxREAL);
  STATUS_OUT = mxCreateDoubleMatrix(1, 1, mxREAL);
  
  double *xmin   = mxGetPr(XMIN_OUT);
  double *fmin   = mxGetPr(FMIN_OUT);
  double *status = mxGetPr(STATUS_OUT);
  
  EXTRA_OUT  = mxCreateStructMatrix(1, 1, 4, extranames);
  mxArray *mxlambda   = mxCreateDoubleMatrix(mrowsA, 1, mxREAL);
  mxArray *mxredcosts = mxCreateDoubleMatrix(mrowsc, 1, mxREAL);
  mxArray *mxtime     = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxArray *mxmem      = mxCreateDoubleMatrix(1, 1, mxREAL);
  
  double *lambda = mxGetPr(mxlambda);
  double *redcosts= mxGetPr(mxredcosts);
  double *time   = mxGetPr(mxtime);
  double *mem    = mxGetPr(mxmem);
  
  int jmpret = setjmp (mark);

  if (jmpret == 0)
    glpk (sense, mrowsc, mrowsA, c, nz, rn,
	       cn, a, b, ctype, freeLB, lb, freeUB,
	       ub, vartype, isMIP, lpsolver, save_pb, save_filename, filetype,
	       xmin, fmin, status, lambda,
	       redcosts, time, mem);

  if (! isMIP)
    {
      mxSetField(EXTRA_OUT,0,extranames[0],mxlambda);
      mxSetField(EXTRA_OUT,0,extranames[1],mxredcosts);
    }

  mxSetField(EXTRA_OUT,0,extranames[2],mxtime);
  mxSetField(EXTRA_OUT,0,extranames[3],mxmem);

  mxFree(rn);
  mxFree(cn);
  mxFree(a);
  mxFree(freeLB);
  mxFree(freeUB);		
  mxFree(ctype);
  mxFree(vartype);
  mxFree(vtype);
  mxFree(extranames);
  mxFree(save_filename);
  mxFree(filetype);
  
  return;
}
