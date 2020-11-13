/*  mpt_getInput.c
  
  Identifies a control law associated to a given state X.

  Requires the explicit controller to be described in a header file called
  mpt_getInput.h (see 'help mpt_exportc' for more details).
  
  Usage:
    region = mpt_getInput(*X, *U)
    
    if "region" is smaller 1 (region < 1), there is no control law associated to
    a given state.

   Please note that all code in this file is provided under the terms of the
   GNU General Public License, which implies that if you include it directly
   into your commercial application, you will need to comply with the license.
   If you feel this is not a good solution for you or your company, feel free 
   to contact me at kvasnica@control.ee.ethz.ch, I can re-license this specific 
   piece of code to you free of charge.
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

#ifndef mpt_getInput_h
#include "mpt_getInput.h"
#endif

static float mpt_getInput(float *X, float *U)
{
    int ix, iu, ic, nc, isinside;
    unsigned long ireg, abspos;
    float hx, region;
    
    abspos = 0;
    region = 0;
    
    /* initialize U to zero*/
    for (iu=0; iu<MPT_NU; iu++) {
        U[iu] = 0;
    }
    
    for (ireg=0; ireg<MPT_NR; ireg++) {
        
        isinside = 1;
        nc = MPT_NC[ireg];
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
        abspos = abspos + MPT_NC[ireg];
    }
    return region;
}

static void mpt_augmentState(float *Xaug, float *X, float *Uprev, float *reference)
{
    /* augments states vector for tracking / deltaU formulation */
    int i;
    if (MPT_TRACKING==1) {
        /* Xaug = [X; U(k-1); reference] */
        for (i=0; i<MPT_NXT; i++)
            Xaug[i] = X[i];
        for (i=0; i<MPT_NU; i++)
            Xaug[MPT_NXT + i] = Uprev[i];
        for (i=0; i<MPT_NREF; i++)
            Xaug[MPT_NXT + MPT_NU + i] = reference[i];
    } else if (MPT_TRACKING==2) {
        /* X = [X; reference] */
        for (i=0; i<MPT_NXT; i++)
            Xaug[i] = X[i];
        for (i=0; i<MPT_NREF; i++)
            Xaug[MPT_NXT + i] = reference[i];
    } else if (MPT_DUMODE==1) {
        /* X = [Xin; U(k-1)] */
        for (i=0; i<MPT_NXT; i++)
            Xaug[i] = X[i];
        for (i=0; i<MPT_NU; i++)
            Xaug[MPT_NXT + i] = Uprev[i];
    } else {
        /* no augmentation necessary, just copy states */
        for (i=0; i<MPT_NX; i++) {
            Xaug[i] = X[i];
        }
    }
}


static void mpt_augmentInput(float *U, const float *Uprev)
{
    int i;
    if ((MPT_TRACKING==1) || (MPT_DUMODE>0))
    {
        /* controllers with MPT_TRACKING=1 or MPT_DUMODE=1 generate
         * deltaU instead of U, therefore to obtain the "true" control action,
         * we need to add previous input at this point
         */
        for (i=0; i<MPT_NU; i++)
        {
            U[i] = U[i] + Uprev[i];
        }
    } 
}
