/*  mpt_searchTree.c
  
  Identifies a control law associated to a given state X using a binary search tree.

  Usage:
    region = mpt_searchTree(*X, *U)
    
    if "region" is smaller 1 (region < 1), there is no control law associated to
    a given state.

   Please note that all code in this file is provided under the terms of the
   GNU General Public License, which implies that if you include it directly
   into your commercial application, you will need to comply with the license.
   If you feel this is not a good solution for you or your company, feel free 
   to contact me at kvasnica@control.ee.ethz.ch, I can re-license this specific 
   piece of code to you free of charge.
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

/* placeholder, do not edit or remove!!! */

#ifdef tmwtypes_h
  /* RTW is used, switch to real_T data types to avoid problems with TLC compilation */
  static long mpt_searchTree(const real_T *X, real_T *U)
#else
  static long mpt_searchTree(const float *X, float *U)
#endif
{
    int ix, iu;
    long node = 1, row;
    float hx, k;
    
    /* initialize U to zero*/
    for (iu=0; iu<MPT_NU; iu++) {
        U[iu] = 0;
    }
    
    /* find region which contains the state x0 */
    while (node > 0) {
        hx = 0;
        row = (node-1)*(MPT_NX+3);
        for (ix=0; ix<MPT_NX; ix++) {
            hx = hx + MPT_ST[row+ix]*X[ix];
        }
        k = MPT_ST[row+MPT_NX];
        
        if ((hx - k) < 0) {
            /* x0 on plus-side of the hyperplane */
            node = (long)MPT_ST[row+MPT_NX+2];
        } else {
            /* x0 on minus-side of the hyperplane */
            node = (long)MPT_ST[row+MPT_NX+1];
        }
    }
    
    node = -node;
    
    /* compute control action associated to state x0 */
    for (iu=0; iu<MPT_NU; iu++) {
        for (ix=0; ix<MPT_NX; ix++) {
            U[iu] = U[iu] + MPT_F[(node-1)*MPT_NX*MPT_NU + iu*MPT_NX + ix]*X[ix];
        }
        U[iu] = U[iu] + MPT_G[(node-1)*MPT_NU + iu];
    }
    
    return node;
}
