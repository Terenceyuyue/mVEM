/* this example illustrates how to embed explicit controllers directly
 * into your C/C++ applications
 *
 * Steps to compile:
 *   1. generate an explicit controller using 'mpt_control'
 *       >> Double_Integrator
 *       >> controller = mpt_control(sysStruct, probStruct);
 *
 *   2. export the explicit controller to C-code using 'mpt_exportc'
 *       >> mpt_exportc(controller);
 *
 *   3. compile this example
 *       >> !gcc mpt_example.c -o mpt_example
 *
 *
 * NOTE: compilation requires "mpt_getInput.c" which you can find in
 *       mpt/extras/simulink directory of your MPT installation
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

#include <stdio.h>

/* this will include also mpt_getInput.h */
#include "mpt_getInput.c"

#define STOP_TIME 20

int main()
{
    int i, iN;
    double X[MPT_NX], X0[MPT_NXT], Xn[MPT_NXT], U[MPT_NU], Uprev[MPT_NU];
    double reference[MPT_NREF+1], region;
    double A[2][2] = { {1.0, 1.0}, {0.0, 1.0}};
    double B[2] = {1.0, 0.5};
    
    
    /* check dimensions, this example is for 2 state variables and 1 input */
    if (MPT_NXT != 2)
    {
        printf("Wrong dimension of state vector!\n");
        return 0;
    }
    if (MPT_NU != 1)
    {
        printf("Wrong dimension of input vector!\n");
        return 0;
    }
        
    /* initialization */
    for (i=0; i<MPT_NXT; i++)
    {
        /* set initial conditions */
        X0[i] = 1;
    }

    for (i=0; i<MPT_NU; i++)
    {
        /* initialize previous input to zero */
        Uprev[i] = 0;
    }
    
    if (MPT_TRACKING>0)
    {
        /* set references if the controller is for tracking */
        reference[0] = 1;
        if (MPT_NREF>1)
        {
            for (i=1; i<MPT_NREF; i++)
            {
                reference[i] = 0;
            }
        }
    }
    
    
    for (iN=0; iN<STOP_TIME; iN++)
    {
        
        printf("time: %d X = [%f %f]\n", iN, X0[0], X0[1]);
        
        /* tracking controllers require the state vector to be augmented.
         * if the controller is not for tracking, the function will just
         * copy X0 to X
         */
        mpt_augmentState(X, X0, Uprev, reference);

        /* obtain control action for a given state. */
        region = mpt_getInput(X, U);
        
        if (region < 1)
        {
            printf("No feasible control law found!\n");
            return 0;
        }
        
        /* tracking controllers return U(k)-U(k-1) instead of U(k), therefore
         * in order to obtain the "true" control action, one needs to add
         * U(k-1) to the result obtained with mpt_getInput. This augmentation
         * can be automatically performed in the following function. we recommend
         * to always use this function, since it only performes the augmentation
         * when necessary.
         */
        mpt_augmentInput(U, Uprev);
                
        /* keep the control input, it will be used in subsequent iteration */
        *Uprev = *U;
        
        /* compute state update */
        Xn[0] = (A[0][0]) * (X0[0]) + (A[0][1]) * (X0[1]) + (B[0]) * (U[0]);
        Xn[1] = (A[1][0]) * (X0[0]) + (A[1][1]) * (X0[1]) + (B[1]) * (U[0]);
        X0[0] = Xn[0];
        X0[1] = Xn[1];
    }
    
}
