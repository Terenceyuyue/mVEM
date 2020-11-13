/*
  Fourier elimination algorithm
  
  Author: Colin Jones
  Email : cnj22@cam.ac.uk
  Date  : April 17, 2004

% ---------------------------------------------------------------------------
% Legal note:
%          This program is free software; you can redistribute it and/or
%          modify it under the terms of the GNU General Public
%          License as published by the Free Software Foundation; either
%          version 2.1 of the License, or (at your option) any later version.
%
%          This program is distributed in the hope that it will be useful,
%          but WITHOUT ANY WARRANTY; without even the implied warranty of
%          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%          General Public License for more details.
% 
%          You should have received a copy of the GNU General Public
%          License along with this library; if not, write to the 
%          Free Software Foundation, Inc., 
%          59 Temple Place, Suite 330, 
%          Boston, MA  02111-1307  USA
%
% ---------------------------------------------------------------------------
*/

#ifndef CTEST
#include "mex.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <memory.h>
#include <math.h>
#include "fourier.h"

//#define TOL 1e-6
//#define QUASI_TOL (1.0*0.0174532925) // how close the normalized normal needs to be to
                                     // another before they are equivalent (in radians)

#define printall 0

extern int DIM;
extern CSet *vset;
extern CSet *hset;

// C   - input constraints
// ax  - project onto
// P   - output constraints
// lnP - number of ontput constraints
mxArray* fourier(double *C, int lnC, int *ax, double TOL, double QUASI_TOL)
{
  int    *rem,lnrem;
  int     i,j,k,l;
  int    *O,*H,*E,*I;
  List    S = NULL;
  List    pos = NULL;
  List    neg = NULL;
  LList   test = NULL;
  LNode  *lntmp;
  double  dtmp, *bdiff;
  int     dredun,itmp;
  Node   *ntmp,*ntmp2;
  List    p,n;
  double  *a,b,scl;
  int    *tvset1,*tvset2;
  int     lnP;
  double *ptrP;
  mxArray *P;
  long    count,rej,quasi,crej,sclrej;

  int      rows,cols;
  mxArray *rnkResult;
  double  *rnkData;
  mxArray *rnkTest;
  int      rank;
  int     *cax = vset->newset();
  int     *indcax = vset->newset();
  vset->setcomp(ax,cax);
  vset->ind(cax,indcax);

  // Axes left to remove
  rem = vset->newset(); 
  vset->setcomp(ax,rem);
  lnrem = vset->card(rem);

  O = vset->newset(); // Officially eliminated variables
  H = hset->newset(); // Historical subsets
  E = vset->newset(); // Effectively eliminated variables
  I = vset->newset(); // Implicitly eliminated variables

  // temporary sets
  tvset1 = vset->newset();
  tvset2 = vset->newset();
  
  // Temp storage
  a = (double*)malloc(sizeof(double)*DIM);

  // Initialize constraints
  for(i=0;i<lnC;i++)
	{
	  ntmp = newnode();
	  scl = 0.0;
	  for(j=0;j<DIM;j++)
		{
		  ntmp->a[j] = C[i+lnC*j];
		  scl += ntmp->a[j] * ntmp->a[j];
		};
	  scl = sqrt(scl);
	  for(j=0;j<DIM;j++) ntmp->a[j] /= scl;
	  ntmp->b = C[i + lnC*DIM] / scl;
	  ntmp->H[i] = 1;
	  S = push(S,ntmp);
	}

  //  PrintConstraints(S);
  while(vset->card(rem))
	{
	  // Select the next variable to eliminate
	  //	  for(j=0;j<DIM;j++) if(rem[j]) {rem[j]=0; break;};
	  for(j=DIM-1;j>=0;j--) if(rem[j]) {rem[j]=0; break;};

	  // 2.1 Suppress from S all the inequations in which the coefficient of x_j is non-zero
	  ntmp = S;
	  while(ntmp)
		{
		  ntmp2 = ntmp->next;
		  if(ntmp->a[j] >  TOL)
			{
			  S = remove(S,ntmp);
			  pos = push(pos,ntmp);
			};
		  if(ntmp->a[j] < -TOL)
			{
			  S = remove(S,ntmp);
			  neg = push(neg,ntmp);
			}
		  ntmp = ntmp2;
		};

	  if(printall)
		{
		  printf("\n\npos\n");
		  PrintConstraints(pos);
		  
		  printf("\n\nneg\n");
		  PrintConstraints(neg);
		  
		  printf("\n\nS\n");
		  PrintConstraints(S);
		}

	  // 2.2
	  if(pos == NULL || neg==NULL) 
		{
		  freelist(pos); pos = NULL;
		  freelist(neg); neg = NULL;
		  continue;
		}
    
	  O[j] = 1;

	  count  = 1;
	  rej    = 0;
	  quasi  = 0;
	  crej   = 0;
	  sclrej = 0;
	  for(p=pos;p;p=p->next)
		for(n=neg;n;n=n->next)
		  {
			count = count + 1;
			if(0)//(count%10000)==0) 
			  {
				printf("\r");
				printf("%2i/%2i ",vset->card(cax)-vset->card(rem),vset->card(cax));
				printf("tot: %8i/%8i ",count,listlen(pos)*listlen(neg));
				printf("rej: %6i (%5.2f%%) ",rej,(double)rej/(double)count*100.0);
				printf("crej: %6i ",crej);
				printf("sclrej: %6i ",sclrej);
				printf("quasi: %6i (%5.2f) ",quasi,100.0*(double)quasi/(double)count);
				printf("keep: %6i",listlen(S));
			  };
				
			// 2.3.1
			H = hset->setunion(p->H,n->H,H);
			E = vset->setunion(p->E,n->E,E);
			E[j] = 1;
			I = vset->setunion(p->I,n->I,I);

			// 2.3.2
			vset->setintersect(I,O,tvset1);
			vset->setunion(E,tvset1,tvset2);
			if(vset->card(tvset2)+1 < hset->card(H))
			  {
				crej++;
				rej++;
				continue;
			  };
			/*
			// 2.3.3
			if(vset->card(E)+1 != hset->card(H))
			  {
				// 2.3.4 Testing H using matricial computation
				rows = hset->card(H);
				cols = vset->card(cax);
				rnkTest = mxCreateDoubleMatrix(rows,cols,mxREAL);
				rnkData = mxGetPr(rnkTest);
				k=0;
				for(i=0;i<lnC;i++)
				  {
					if(H[i])
					  {
						for(l=0;l<cols;l++)
						  rnkData[k+rows*l] = C[i+lnC*indcax[l]];
						k++;
					  };
				  };

				mexCallMATLAB(1,&rnkResult,1,&rnkTest,"rank");
				rank = int(mxGetScalar(rnkResult));
  				if(rank != hset->card(H)-1)
				  {
					rej = rej + 1;
					continue;
				  };
			  };
			*/

			// 2.3.5 Suppress some strongly redundant inequalities
			// test for quasi-redundancies(?)
			scl = 0.0;
			for(i=0;i<DIM;i++)
			  {
				a[i] = p->a[j] * n->a[i] - n->a[j] * p->a[i];
				scl += a[i]*a[i];
			  };
			scl = sqrt(scl);
			if(scl < TOL)
			  {
				sclrej++;
				rej++;
				continue;
			  };
			b = (p->a[j] * n->b - n->a[j] * p->b)/scl;
			for(i=0;i<DIM;i++) a[i] /= scl;

			if(printall)
			  {
				printf("p->a = ["); for(i=0;i<DIM;i++) printf("%5.2f ",p->a[i]); printf("]\n");
				printf("n->a = ["); for(i=0;i<DIM;i++) printf("%5.2f ",n->a[i]); printf("]\n");
				printf("SCL = %f\n",scl);
			  }

			// Eliminate quasi-redundancies...
			if(S != NULL)
			  {
				// Find all normals in S which are the same as d
				test = NULL;
				for(ntmp=S;ntmp;ntmp=ntmp->next)
				  {
					dtmp = 0;
					for(i=0;i<DIM;i++)
					  dtmp += (ntmp->a[i]*a[i]);
					if(acos(dtmp) < QUASI_TOL)
					  test = Lpush(test,ntmp);
				  };

				if(test != NULL)
				  {
					// Test if d is quasi-redundant
					itmp = Llistlen(test);
					bdiff  = (double*)malloc(sizeof(double)*itmp);
					dredun = 0;
					for(i=0,lntmp=test;lntmp;lntmp=lntmp->next,i++)
					  {
						bdiff[i] = 2.0*(lntmp->pNode->b - b)/(fabs(lntmp->pNode->b)+b);
						if(bdiff[i] < -0.01)
						  {
							dredun = 1;
							break;
						  };
					  };
					if(dredun == 1)
					  {
						free(bdiff);
					    quasi++;
						rej++;
						continue;
					  };
					
					// Find all redundancies in S wrt d
					for(i=0,lntmp=test;lntmp;lntmp=lntmp->next,i++)
					  {
						if(bdiff[i] > 0.01)
						  {
							// It's redundant - remove it
							S = remove(S,lntmp->pNode);
							freenode(lntmp->pNode);
						    quasi++;
							rej++;
						  };
					  };
					free(bdiff);
					Lfreelist(test);
				  };
			  };

			// 2.3.6 Put the equation from x and y in S
			ntmp = newnode();
			memcpy(ntmp->a,a,sizeof(double)*DIM);
			ntmp->b = b;
			hset->setcopy(ntmp->H,H);
			vset->setcopy(ntmp->E,E);
	  
			// you're implicitally eliminated if (e.g) you had x4 before and you
			// don't anymore by the elimination of x1
			for(i=0;i<DIM;i++)
			  if((fabs(p->a[i])>TOL || fabs(n->a[i])>TOL) && (fabs(a[i]) < TOL) && (rem[i] == 1))
				{
				  //printf("Got an implicit elimination!\n");
				  I[i] = 1;
				};
			vset->setcopy(ntmp->I,I);

			S = push(S,ntmp);
		  };
	  if(0)
	    {
	      printf("\r");
	      printf("%2i/%2i ",vset->card(cax)-vset->card(rem),vset->card(cax));
	      printf("tot: %8i/%8i ",count,listlen(pos)*listlen(neg));
	      printf("rej: %6i (%5.2f%%) ",rej,(double)rej/(double)count*100.0);
	      printf("crej: %6i ",crej);
	      printf("sclrej: %6i ",sclrej);
	      printf("quasi: %6i (%5.2f) ",quasi,100.0*(double)quasi/(double)count);
	      printf("keep: %6i",listlen(S));
	      printf("\n");
	    };

	  freelist(pos); pos = NULL;
	  freelist(neg); neg = NULL;

	  if(printall)
		{
		  printf("\n\nCurr S\n");
		  PrintConstraints(S);
		}
	}
  if(printall)
	{
	  printf("\n\nFinal S\n");
	  PrintConstraints(S);
	  printf("\n\n\n");
	}

  lnP = listlen(S);
  int pdim = vset->card(ax);

  indcax = vset->ind(ax,indcax);

  // Return the projection
  // Note: this is stored column-wise to match matlab standards
  P = mxCreateDoubleMatrix(lnP,pdim+1,mxREAL);
  ptrP = mxGetPr(P);

  for(i=0,ntmp=S;ntmp;ntmp=ntmp->next,i++)
	{
	  for(j=0;j<pdim;j++)
		ptrP[i+lnP*j] = ntmp->a[indcax[j]];
	  ptrP[i+lnP*pdim] = ntmp->b;
	};

  vset->freeset(rem);
  vset->freeset(O);
  hset->freeset(H);
  vset->freeset(E);
  vset->freeset(I);
  vset->freeset(tvset1);
  vset->freeset(tvset2);
  vset->freeset(indcax);
  vset->freeset(cax);
  freelist(S);
  free(a);

  if(printall)
	printf("\n");

  return P;
}


void PrintConstraints(List S)
{
  int j;
  Node *n;
  n = S;
  while(n != NULL)
	{
	  for(j=0;j<DIM;j++) printf("%5.2f ",n->a[j]);
	  printf("|%5.2f ",n->b);
	  //	  printf("H"); 
	  hset->printset(n->H);
	  //	  printf(" E"); vset->printset(n->E);
	  //	  printf(" I"); vset->printset(n->I);
	  printf("\n");
	  n = n->next;
	};
}

