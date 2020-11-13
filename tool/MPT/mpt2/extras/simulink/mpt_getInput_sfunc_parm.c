/*  mpt_getInput_sfunc_parm.c
  
  C-code S-function for simulation of explicit controllers.
  
  Parameters of the explicit controller are passed as input arguments
  (see 'help mpt_getSfuncParam')
  
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

#define S_FUNCTION_NAME  mpt_getInput_sfunc_parm   /*Name of the S-function file*/
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
#define NPARAMS 16


static double mpt_getInput(SimStruct *S, double *X, double *U);


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
    int MPT_NR, MPT_NU, MPT_NXT, MPT_NREF, MPT_TRACKING;
    
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

    MPT_NR = (int)mxGetPr(ssGetSFcnParam(S, 5))[0];
    MPT_NU = (int)mxGetPr(ssGetSFcnParam(S, 7))[0];
    MPT_NXT = (int)mxGetPr(ssGetSFcnParam(S, 9))[0];
    MPT_NREF = (int)mxGetPr(ssGetSFcnParam(S, 10))[0];
    MPT_TRACKING = (int)mxGetPr(ssGetSFcnParam(S, 13))[0];
        
    /* no states, all computation will be done in output section */
    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, 0);
    
    /* set dimension of the work vector - ABSTOL, TS, previous input */
    ssSetNumRWork(S, MPT_NU + 2);
    
    /* pointer work vector will have 5 elements - H, K, F, K, NC vectors */
    ssSetNumPWork(S, 5);
    
    /* integer work vector will have 8 elements:
       NR, NX, NU, NY, NXT, NREF, DUMODE, TRACKING, BREAK_WHEN_INFEASIBLE
     */
    ssSetNumIWork(S, 9);
    
    /* one input port: state x(k) + vector of references */
    ssSetNumInputPorts(S, 1); 
    
    /* one output port: control action */
    ssSetNumOutputPorts(S, 1);

    /* width of input vector - number of states of the original problem. 
     * note that tracking includes additional states, here we do not consider them
     */
    if (MPT_TRACKING>0) {
        ssSetInputPortWidth(S, 0, MPT_NXT + MPT_NREF);
    } else {
        /* dimension extended by one to allow empty reference to be passed */
        ssSetInputPortWidth(S, 0, MPT_NXT + 1);
    }
    
    ssSetInputPortDirectFeedThrough(S, 0, 1);
    
    /* width of output - number of control actions + active region*/
    ssSetOutputPortWidth(S, 0, MPT_NU + 1);
    
    ssSetNumSampleTimes(S, 1);

    /* Take care when specifying exception free code - see sfuntmpl.doc */
    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);
    
    #ifdef MPT_DEBUG
    printf("mdlInitializeSizes nr: %d\n", (int)mxGetPr(ssGetSFcnParam(S, 5))[0]);
    #endif
    
}

#define MDL_START  /* Change to #undef to remove function */
#if defined(MDL_START) 
static void mdlStart(SimStruct *S)
{
    double *buffer;
    int i;
    size_t numel;
    #ifdef MPT_DEBUG
    int MPT_NU;
    #endif

    /* set elements of the real work vector */
    ssSetRWorkValue(S, 0, mxGetPr(ssGetSFcnParam(S, 14))[0]); /* ABSTOL */
    ssSetRWorkValue(S, 1, mxGetPr(ssGetSFcnParam(S, 11))[0]); /* TS */
    
    /* set elements of the integer work vector */
    ssSetIWorkValue(S, 0, (int)mxGetPr(ssGetSFcnParam(S, 5))[0]); /* NR */
    ssSetIWorkValue(S, 1, (int)mxGetPr(ssGetSFcnParam(S, 6))[0]); /* NX */
    ssSetIWorkValue(S, 2, (int)mxGetPr(ssGetSFcnParam(S, 7))[0]); /* NU */
    ssSetIWorkValue(S, 3, (int)mxGetPr(ssGetSFcnParam(S, 8))[0]); /* NY */
    ssSetIWorkValue(S, 4, (int)mxGetPr(ssGetSFcnParam(S, 9))[0]); /* NXT */
    ssSetIWorkValue(S, 5, (int)mxGetPr(ssGetSFcnParam(S, 10))[0]); /* NREF */
    ssSetIWorkValue(S, 6, (int)mxGetPr(ssGetSFcnParam(S, 12))[0]); /* DUMODE */
    ssSetIWorkValue(S, 7, (int)mxGetPr(ssGetSFcnParam(S, 13))[0]); /* TRACKING */
    ssSetIWorkValue(S, 8, (int)mxGetPr(ssGetSFcnParam(S, 15))[0]); /* BREAK_WHEN_INFEASIBLE */
    
    #ifdef MPT_DEBUG
    printf("mdlStart nr: %d\n", (int)mxGetPr(ssGetSFcnParam(S, 5))[0]);
    printf("mdlStart dumode: %f\n", (double)mxGetPr(ssGetSFcnParam(S, 12))[0]);
    #endif
    
    /* allocate memory for H, K, F, G, NC vectors */
    for (i=0; i<5; i++)
    {
        numel = mxGetNumberOfElements(ssGetSFcnParam(S, i))+1;
        #ifdef MPT_DEBUG
        printf("allocating %f elements for array #%d\n", (double)numel, i+1);
        #endif
        
        #ifdef USE_MX_ALLOC
        buffer=(double *)mxCalloc(numel, sizeof(double));
        #else
        buffer=(double *)calloc(numel, sizeof(double));
        #endif
        
        if (buffer==NULL) {
            ssSetErrorStatus(S,"mpt_getInput_sfunc_parm: could not allocate memory!\n");
            return;
        }
        buffer = mxGetPr(ssGetSFcnParam(S, i));
        ssSetPWorkValue(S,i,(void *)buffer);
    }

}
#endif


static void mdlInitializeSampleTimes(SimStruct *S)
{
    real_T MPT_TS;
    
    MPT_TS = mxGetPr(ssGetSFcnParam(S, 11))[0];
    
    /* set sampling time */
    ssSetSampleTime(S, 0, MPT_TS);
    ssSetOffsetTime(S, 0, 0.0);
    
    #ifdef MPT_DEBUG
    printf("mdlInitializeSampleTimes Ts: %e nr: %d\n", mxGetPr(ssGetSFcnParam(S, 11))[0], (int)mxGetPr(ssGetSFcnParam(S, 5))[0]);
    #endif
}


#define MDL_INITIALIZE_CONDITIONS

static void mdlInitializeConditions(SimStruct *S)
{
    int_T MPT_NU;
    int i;

    #ifdef MPT_DEBUG
    printf("mdlInitializeConditions nr: %d\n", (int)mxGetPr(ssGetSFcnParam(S, 5))[0]);
    #endif
    
    MPT_NU = (int)mxGetPr(ssGetSFcnParam(S, 7))[0];

    /* set previous input to zero. note that first two elements
       of the RWork vector are reserved for the absolute tolerance and for TS
     */
    for (i=0; i<MPT_NU; i++)
        ssSetRWorkValue(S, i+2, 0);
    
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
    
    real_T            	*u   = ssGetOutputPortRealSignal(S,0);
	InputRealPtrsType    Xin = ssGetInputPortRealSignalPtrs(S,0);
    double region, U[100], X[100], REF[100], Uprev[100];
    int_T MPT_NU, MPT_NXT, MPT_NREF, MPT_DUMODE, MPT_TRACKING, MPT_NX;
    int_T MPT_BREAK_WHEN_INFEASIBLE;
    int i;
    #ifdef MPT_DEBUG
    int_T MPT_NR;
    #endif

    #ifdef MPT_DEBUG
    MPT_NR = ssGetIWorkValue(S, 0);
    #endif
    MPT_NX = ssGetIWorkValue(S, 1);
    MPT_NU = ssGetIWorkValue(S, 2);
    MPT_NXT = ssGetIWorkValue(S, 4);
    MPT_NREF = ssGetIWorkValue(S, 5);
    MPT_DUMODE = ssGetIWorkValue(S, 6);
    MPT_TRACKING = ssGetIWorkValue(S, 7);
    MPT_BREAK_WHEN_INFEASIBLE = ssGetIWorkValue(S, 8);
    
    for (i=0; i<MPT_NU; i++) {
        /* note that RWork[0]=ABSTOL, RWork[1]=TS, RWork[2..2+NU]=Uprev */
        Uprev[i] = ssGetRWorkValue(S, i+2); 
    }
    
    if (MPT_TRACKING>0) {
        /* extract references from input vector */
        for (i=0; i<MPT_NREF; i++)
        {
            REF[i] = (double)*Xin[MPT_NXT + i];
        }
    }
    
    #ifdef MPT_DEBUG
    printf("\n-----------------------\n");
    printf("Tracking: %d dumode: %d nr: %d\n", (int)MPT_TRACKING, (int)MPT_DUMODE, (int)MPT_NR);
    for (i=0; i<MPT_NU; i++)
        printf("Uprev(%d/%d) = %f\n", i+1, MPT_NU, Uprev[i]);
    
    for (i=0; i<MPT_NREF; i++)
        printf("Ref(%d/%d) = %f\n", i+1, MPT_NREF, REF[i]);
    #endif
    
    /* augment state vector to deal with tracking and deltaU formulation */
    if (MPT_TRACKING==1) {
        /* X = [Xin; U(k-1); REF] */
        for (i=0; i<MPT_NXT; i++)
            X[i] = (double)*Xin[i];
        for (i=0; i<MPT_NU; i++)
            X[MPT_NXT + i] = Uprev[i];
        for (i=0; i<MPT_NREF; i++)
            X[MPT_NXT + MPT_NU + i] = REF[i];
    } else if (MPT_TRACKING==2) {
        /* X = [Xin; REF] */
        for (i=0; i<MPT_NXT; i++)
            X[i] = (double)*Xin[i];
        for (i=0; i<MPT_NREF; i++)
            X[MPT_NXT + i] = REF[i];
    } else if (MPT_DUMODE==1) {
        /* X = [Xin; U(k-1)] */
        for (i=0; i<MPT_NXT; i++)
            X[i] = (double)*Xin[i];
        for (i=0; i<MPT_NU; i++)
            X[MPT_NXT + i] = Uprev[i];
    } else {
        for (i=0; i<MPT_NX; i++) {
            X[i] = (double)*Xin[i];
        }
    }
    
    #ifdef MPT_DEBUG
    for (i=0; i<MPT_NX; i++)
        printf("X(%d/%d) = %f\n", i+1, MPT_NX, X[i]);
    #endif
    
    /* get control law */
    region = mpt_getInput(S, X, U);
        
    #ifdef MPT_DEBUG
    printf("Region: %d\n", (int)region);
    #endif
        
    /* check if control law was found, if not, stop the simulation */
    if (region<1) {
        #ifdef MATLAB_MEX_FILE
        if (MPT_BREAK_WHEN_INFEASIBLE) {
            ssSetErrorStatus(S, "No feasible control law found!");
        }
        #else
        /* printf("No feasible control law found at time %f!\n", (double)tid); */
        #endif
        /* return; */
    }
    
    /* output control action */
    for (i=0; i<MPT_NU; i++) {
        if ((MPT_TRACKING==1) || (MPT_DUMODE==1)) {
            /* tracking controller generates deltaU! therefore we add U(k-1) */
            if (region>0) {
                u[i] = (real_T)U[i] + Uprev[i];
            } else {
                u[i] = (real_T)Uprev[i];
            }
        } else {
            if (region>0) {
                u[i] = (real_T)U[i];
            } else {
                u[i] = (real_T)Uprev[i];
            }
        }
        /* note that RWork[0]=ABSTOL, RWork[1]=TS, RWork[2..2+NU]=Uprev */
        ssSetRWorkValue(S, i+2, (double)u[i]);
    }
    
    u[MPT_NU] = (int_T)region;
    
    #ifdef MPT_DEBUG
    for (i=0; i<MPT_NU; i++)
        printf("U(%d/%d) = %f\n", i+1, MPT_NU, (double)u[i]);
    #endif
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
    for (i=0; i<5; i++)
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

static double mpt_getInput(SimStruct *S, double *X, double *U)
{
    int ix, iu, ic, nc, isinside;
    unsigned int ireg;
    unsigned long abspos;
    double hx, region;
    #ifdef MPT_DEBUG_2
    double diff;
    #endif
    double * MPT_H, * MPT_K, * MPT_F, * MPT_G, * MPT_NC, MPT_ABSTOL;
    int_T MPT_NR, MPT_NX, MPT_NU;
    
    /* extract data from work vectors */
    MPT_NR = (unsigned int)ssGetIWorkValue(S, 0);
    MPT_NX = (int)ssGetIWorkValue(S, 1);
    MPT_NU = (int)ssGetIWorkValue(S, 2);
    MPT_ABSTOL = ssGetRWorkValue(S, 0);
    MPT_H = (double *)ssGetPWorkValue(S,0);
    MPT_K = (double *)ssGetPWorkValue(S,1);
    MPT_F = (double *)ssGetPWorkValue(S,2);
    MPT_G = (double *)ssGetPWorkValue(S,3);
    MPT_NC = (double *)ssGetPWorkValue(S,4);

        
    abspos = 0;
    region = 0;
    
    /* initialize U to zero*/
    for (iu=0; iu<MPT_NU; iu++) {
        U[iu] = 0;
    }
    
    for (ireg=0; ireg<MPT_NR; ireg++) {
        #ifdef MPT_DEBUG_2
        printf("\tChecking region: %d\n", ireg+1);
        #endif
        
        isinside = 1;
        nc = (int)MPT_NC[ireg];
        for (ic=0; ic<nc; ic++) {
            hx = 0;
            for (ix=0; ix<MPT_NX; ix++) {
                hx = hx + MPT_H[abspos*MPT_NX+ic*MPT_NX+ix]*X[ix];
            }
            if ((hx - MPT_K[abspos+ic]) > MPT_ABSTOL) {
                /* constraint is violated, continue with next region */
                #ifdef MPT_DEBUG_2
                diff = hx - MPT_K[abspos+ic];
                printf("\t\tConstraint %d not satisfied with difference %e, abs_tol %e\n", ic + 1, diff, MPT_ABSTOL);
                #endif
                isinside = 0;
                break;
            } 
        }
        if (isinside==1) {
            /* state belongs to this region, extract control law and exit */
            region = ireg + 1;
            for (iu=0; iu<MPT_NU; iu++) {
                for (ix=0; ix<MPT_NX; ix++) {
                    U[iu] = U[iu] + MPT_F[ireg*MPT_NX*MPT_NU + iu*MPT_NX + ix]*X[ix];
                }
                U[iu] = U[iu] + MPT_G[ireg*MPT_NU + iu];
            }
            return region;
        }
        abspos = abspos + (int)MPT_NC[ireg];
    }
    #ifdef MPT_DEBUG_2
    printf("\tNo region found!!!\n");
    #endif
    return region;
}
