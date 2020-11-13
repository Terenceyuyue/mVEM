/*  mpt_getInput_sfunc.c
  
  Pure C-code S-function for simulation of explicit controllers.
  
  Requires the explicit controller to be described in a header file called
  mpt_getInput.h (see 'help mpt_exportc' for more details).
  
  Usage:
    mpt_exportc(ctrl)
    mex mpt_getInput_sfunc.c
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

#define S_FUNCTION_NAME  mpt_getInput_sfunc   /*Name of the S-function file*/
#define S_FUNCTION_LEVEL 2		/*Level 2 S-functions allow multi-port status*/

#include "simstruc.h"			/*Header where different routines are defined */

#include "mpt_getInput.c"

/* previous input */
static double Uprev[MPT_NU];

static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumSFcnParams(S, 0);
						
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        return; /* Parameter mismatch will be reported by Simulink */
    }

    /* no states, all computation will be done in output section */
    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, 0);
    
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
    
    /* width of output - number of control actions */
    ssSetOutputPortWidth(S, 0, MPT_NU);
    
    ssSetNumSampleTimes(S, 1);

    /* Take care when specifying exception free code - see sfuntmpl.doc */
    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);
    
}


static void mdlInitializeSampleTimes(SimStruct *S)
{
    /* set sampling time */
    ssSetSampleTime(S, 0, MPT_TS);
    ssSetOffsetTime(S, 0, 0.0);
}


#define MDL_INITIALIZE_CONDITIONS

static void mdlInitializeConditions(SimStruct *S)
{
    /* reset previous input to zero */
    int i;
    for (i=0; i<MPT_NU; i++) {
        Uprev[i] = 0;
    }
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
    real_T            	*u   = ssGetOutputPortRealSignal(S,0);
	InputRealPtrsType    Xin   = ssGetInputPortRealSignalPtrs(S,0);
    static double region, U[MPT_NU], X[MPT_NX], REF[MPT_NREF+1];
    int i;

    if (MPT_TRACKING>0) {
        /* extract references from input vector */
        for (i=0; i<MPT_NREF; i++)
        {
            REF[i] = *Xin[MPT_NXT + i];
        }
    }
    
    /* augment state vector to deal with tracking and deltaU formulation */
    if (MPT_TRACKING==1) {
        /* X = [Xin; U(k-1); REF] */
        for (i=0; i<MPT_NXT; i++)
            X[i] = *Xin[i];
        for (i=0; i<MPT_NU; i++)
            X[MPT_NXT + i] = Uprev[i];
        for (i=0; i<MPT_NREF; i++)
            X[MPT_NXT + MPT_NU + i] = REF[i];
    } else if (MPT_TRACKING==2) {
        /* X = [Xin; REF] */
        for (i=0; i<MPT_NXT; i++)
            X[i] = *Xin[i];
        for (i=0; i<MPT_NREF; i++)
            X[MPT_NXT + i] = REF[i];
    } else if (MPT_DUMODE==1) {
        /* X = [Xin; U(k-1)] */
        for (i=0; i<MPT_NXT; i++)
            X[i] = *Xin[i];
        for (i=0; i<MPT_NU; i++)
            X[MPT_NXT + i] = Uprev[i];
    } else {
        for (i=0; i<MPT_NX; i++) {
            X[i] = *Xin[i];
        }
    }
    
    /* get control law */
    region = mpt_getInput(X, U);
        
    /* check if control law was found, if not, stop the simulation */
    if (region<1) {
        ssSetErrorStatus(S, "No feasible control law found!");
    }
    
    /* output control action */
    for (i=0; i<MPT_NU; i++) {
        if ((MPT_TRACKING==1) || (MPT_DUMODE==1)) {
            /* tracking controller generates deltaU! therefore we add U(k-1) */
            u[i] = (real_T)U[i] + Uprev[i];
        } else {
            u[i] = (real_T)U[i];
        }
        Uprev[i] = (double)u[i];
    }
}


/* Function: mdlTerminate =====================================================
 * Abstract:
 *    No termination needed, but we are required to have this routine.
 */
static void mdlTerminate(SimStruct *S)
{
}

/*End of file necessary includes*/

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif