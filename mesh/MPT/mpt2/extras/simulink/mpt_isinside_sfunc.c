/*  mpt_isinside_sfunc.c
  
  C-code S-function which returns true if a given point belongs
  to certain polytope or to a polytope array
  
*/

/* Copyright (C) 2005 by Michal Kvasnica (kvasnica@control.ee.ethz.ch) */

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

#define S_FUNCTION_NAME  mpt_isinside_sfunc  /*Name of the S-function file*/
#define S_FUNCTION_LEVEL 2		/*Level 2 S-functions allow multi-port status*/

/* #define MPT_DEBUG */

#define USE_MX_ALLOC /* this tells to use mxCalloc instead of calloc */

#include "simstruc.h"			/*Header where different routines are defined */

#ifdef USE_MX_ALLOC
    #ifdef MATLAB_MEX_FILE
        #include "matrix.h"
    #endif
    #include <stdlib.h>
#else
    #include "malloc.h"
#endif

#ifdef MPT_DEBUG
#include <stdio.h>
#endif

/* number of input parameters */
#define NPARAMS 6

static double sub_mpt_isinside(SimStruct *S, double *X)
{
    int ix, ic, nc, isinside;
    unsigned int ireg;
    unsigned long abspos;
    double hx, region;
    double * MPT_H, * MPT_K, * MPT_NC, MPT_ABSTOL;
    int_T MPT_NR, MPT_NX;

    /* extract data from work vectors */
    MPT_NR = (unsigned int)ssGetIWorkValue(S, 0);
    MPT_NX = (int)ssGetIWorkValue(S, 1);
    MPT_ABSTOL = ssGetRWorkValue(S, 0);
    MPT_H = (double *)ssGetPWorkValue(S,0);
    MPT_K = (double *)ssGetPWorkValue(S,1);
    MPT_NC = (double *)ssGetPWorkValue(S,2);

    abspos = 0;
    region = 0;
    
    #ifdef MPT_DEBUG
    printf("isinside: X = [%f %f]\n", X[0], X[1]);
    #endif
    
    for (ireg=0; ireg<MPT_NR; ireg++) {
        isinside = 1;
        nc = (int)MPT_NC[ireg];
        for (ic=0; ic<nc; ic++) {
            hx = 0;
            for (ix=0; ix<MPT_NX; ix++) {
                hx = hx + MPT_H[abspos*MPT_NX+ic*MPT_NX+ix]*X[ix];
            }
            if ((hx - MPT_K[abspos+ic]) > MPT_ABSTOL) {
                /* constraint is violated, continue with next region */
                isinside = 0;
                break;
            } 
        }
        if (isinside==1) {
            region = ireg + 1;
            /* state belongs to this region, extract control law and exit */
            return region;
        }
        abspos = abspos + (int)MPT_NC[ireg];
    }

    return region;
}

#define IS_PARAM_DOUBLE(pVal) (mxIsNumeric(pVal) && !mxIsLogical(pVal) && !mxIsEmpty(pVal) && !mxIsSparse(pVal) && !mxIsComplex(pVal) && mxIsDouble(pVal))

#define MDL_CHECK_PARAMETERS
#if defined(MDL_CHECK_PARAMETERS) && defined(MATLAB_MEX_FILE)
   /* Function: mdlCheckParameters =============================================
    * Abstract:
    *    Validate our parameters to verify they are okay.
    */
static void mdlCheckParameters(SimStruct *S)
{
    
        return;
}
#endif /* MDL_CHECK_PARAMETERS */

static void mdlInitializeSizes(SimStruct *S)
{
    int MPT_NX;
    
    ssSetNumSFcnParams(S, NPARAMS);  /* Number of expected parameters */
    #ifdef MATLAB_MEX_FILE
    if (ssGetNumSFcnParams(S) == ssGetSFcnParamsCount(S)) {
        mdlCheckParameters(S);
        if (ssGetErrorStatus(S) != NULL) {
            return;
        }
    } else {
        return; /* Parameter mismatch will be reported by Simulink */
    }
    #endif

    MPT_NX = (int)mxGetPr(ssGetSFcnParam(S, 1))[0];
        
    /* no states, all computation will be done in output section */
    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, 0);

    /* real work vector will have 1 element - MPT_ABSTOL */
    ssSetNumRWork(S, 1);

    /* pointer work vector will have 3 elements - H, K, NC vectors */
    ssSetNumPWork(S, 3);
    
    /* integer work vector will have 2 elements:
       NR, NX
     */
    ssSetNumIWork(S, 2);
    
    /* one input port: state x(k) + vector of references */
    ssSetNumInputPorts(S, 1); 
    
    /* one output port: control action */
    ssSetNumOutputPorts(S, 1);

    /* width of input vector */
    ssSetInputPortWidth(S, 0, MPT_NX);
    
    ssSetInputPortDirectFeedThrough(S, 0, 1);
    
    /* width of output - index of a polytope */
    ssSetOutputPortWidth(S, 0, 1);
    
    /* ssSetNumSampleTimes(S, 1); */

    /* Take care when specifying exception free code - see sfuntmpl.doc */
    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);
    
}

#define MDL_START  /* Change to #undef to remove function */
#if defined(MDL_START) 
static void mdlStart(SimStruct *S)
{
    double *buffer;
    int i;
    size_t numel;

    /* set elements of the integer work vector */
    
    ssSetIWorkValue(S, 0, (int)mxGetPr(ssGetSFcnParam(S, 0))[0]); /* NR */
    ssSetIWorkValue(S, 1, (int)mxGetPr(ssGetSFcnParam(S, 1))[0]); /* NX */
    
    ssSetRWorkValue(S, 0, (double)mxGetPr(ssGetSFcnParam(S, 2))[0]); /* MPT_ABSTOL */
    
    /* allocate memory for H, K, NC vectors */
    for (i=0; i<3; i++)
    {
        numel = mxGetNumberOfElements(ssGetSFcnParam(S, i+3))+1;
        #ifdef MPT_DEBUG
        printf("allocating %f elements for array #%d\n", (double)numel, i+1);
        #endif
        
        #ifdef USE_MX_ALLOC
        buffer=(double *)mxCalloc(numel, sizeof(double));
        #else
        buffer=(double *)calloc(numel, sizeof(double));
        #endif
        
        if (buffer==NULL) {
            ssSetErrorStatus(S,"mpt_isinside_sfunc: could not allocate memory!\n");
            return;
        }
        buffer = mxGetPr(ssGetSFcnParam(S, i+3));
        ssSetPWorkValue(S,i,(void *)buffer);
    }

}
#endif


static void mdlInitializeSampleTimes(SimStruct *S)
{
    /* set sampling time */
    ssSetSampleTime(S, 0, INHERITED_SAMPLE_TIME);
    /* ssSetOffsetTime(S, 0, FIXED_IN_MINOR_STEP_OFFSET); */
    
}


#undef MDL_INITIALIZE_CONDITIONS

static void mdlInitializeConditions(SimStruct *S)
{
   
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
    
    real_T            	*R   = ssGetOutputPortRealSignal(S,0);
	InputRealPtrsType    Xin = ssGetInputPortRealSignalPtrs(S,0);
    double region;
    
    #ifdef MPT_DEBUG
    printf("X = [%f %f]\n", (double)*Xin[0], (double)*Xin[1]);
    #endif
    
    region = sub_mpt_isinside(S, *Xin);
    
    #ifdef MPT_DEBUG
    printf("mdlOutputs: %f\n", region);
    #endif
    *R = (real_T)region;
}


/* Function: mdlTerminate =====================================================
 * Abstract:
 *    No termination needed, but we are required to have this routine.
 */
static void mdlTerminate(SimStruct *S)
{
    #ifndef USE_MX_ALLOC
    double *buffer;
    int i;
    
    /* free allocated memory */
    for (i=0; i<3; i++)
    {
        buffer=(double *)ssGetPWorkValue(S,i);
        if (buffer!=NULL) {
            free(buffer);
        }
    }
    #endif
}

/*End of file necessary includes*/

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif

