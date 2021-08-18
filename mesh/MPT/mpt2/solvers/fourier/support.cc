#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include "fourier.h"

extern int DIM;
extern CSet *vset;
extern CSet *hset;

Node* newnode()
{
  int i;

  Node *tmp = (Node*)malloc(sizeof(Node));
  
  tmp->a = (double*)calloc(DIM,sizeof(double));
  for(i=0;i<DIM;i++) tmp->a[i] = 0.0;

  tmp->b = 0.0;
  tmp->H = hset->newset();
  tmp->E = vset->newset();
  tmp->I = vset->newset();
  
  tmp->prev = NULL;
  tmp->next = NULL;
  
  return tmp;
}

void freenode(Node *N)
{
  hset->freeset(N->H);
  vset->freeset(N->E);
  vset->freeset(N->I);
  free(N->a);
  free(N);
}

void freelist(List L)
{
  Node *next;
  while(L != NULL)
	{
	  next = L->next;
	  freenode(L);
	  L = next;
	};
}

List push(List L, Node *n)
{
  n->next = L;
  n->prev = NULL;
  if(L == NULL) return n;
  L->prev = n;
  return n;
}

// THIS MIGHT BE FUCKED...
Node *pop(List L)
{
  if(L == NULL) return NULL;
  Node *n = L;
  L = L->next;
  if(L != NULL) L->prev = NULL;
  n->next = NULL;
  return n;
}

// Remove node n from list L
List remove(List L, Node *n)
{
  if(L == NULL) return NULL;
  if(n->next != NULL) n->next->prev = n->prev;
  if(n->prev != NULL) n->prev->next = n->next;

  if(n == L) return L->next;
  return L;
}

int listlen(List L)
{
  Node *tmp;
  int i=0;
  for(tmp=L;tmp;tmp=tmp->next) i++;
  return i;
}

/* Set operations */
bool CSet::setisempty(int *A)
{
  int i;
  for(i=0;i<LN;i++)
	if(A[i] != 0)
	  return false;
  return true;
}

int* CSet::newset() 
{
  int i;
  int *tmp=(int*)calloc(LN,sizeof(int));
  for(i=0;i<LN;i++) tmp[i] = 0;
  return tmp;
}

void CSet::freeset(int *A) 
{
  free(A);
}

int* CSet::setcopy(int *C, int *A)
{
//  for(int i=0;i<LN;i++) C[i]=A[i];
  memcpy(C,A,sizeof(int)*LN);
  return C;
}

/* Returns values in A that are not in B */
int* CSet::setdiff(int *A, int *B, int *C)
{
  for(int i=0;i<LN;i++) if(A[i]==1 && B[i]==0) C[i]=1; else C[i]=0;
  return C;
}

int* CSet::setunion(int *A, int *B, int *C)
{
  for(int i=0;i<LN;i++) if(A[i] || B[i]) C[i]=1; else C[i]=0;
  return C;
}

int* CSet::setintersect(int *A, int *B, int *C)
{
  for(int i=0;i<LN;i++) if(A[i]==1 && B[i]==1) C[i]=1; else C[i]=0;
  return C;
}

int* CSet::setxor(int *A, int *B, int *C)
{
  for(int i=0;i<LN;i++) if(A[i]==1 && B[i]==0 || A[i]==0 && B[i]==1) C[i]=1; else C[i]=0;
  return C;
}

int* CSet::setcomp(int *A, int *C)
{
  for(int i=0;i<LN;i++) if(A[i]) C[i]=0; else C[i]=1;
  return C;
}

/* Cardinality */
int CSet::card(int *A)
{
  int tmp = 0;
  for(int i=0;i<LN;i++) tmp += A[i];
  return tmp;
}

void CSet::printset(int *A)
{
  printf("[");
  for(int i=0;i<LN;i++) if(A[i]) printf("%i ",i);
  printf("]");
}

/* Return the indicies of A in C */
int* CSet::ind(int *A, int *C)
{
  int i,j;
  for(i=0,j=0;i<LN;i++) if(A[i]) C[j++] = i;
  for(;j<LN;j++) C[j]=0;
  return C;
}

/*
void testset()
{
  int *A = newset();
  int *B = newset();
  int *C = newset();
  printf("Testing set operations:\n");
  A[2] = 1;
  A[5] = 1;
  B[1] = 1;
  B[5] = 1;
  printf("A: "); printset(A); printf("\n");
  printf("B: "); printset(B); printf("\n");
  
  printf("union: "); printset(setunion(A,B,C)); printf("\n");
  printf("setdiff: "); printset(setdiff(A,B,C)); printf("\n");
  printf("setintersect: "); printset(setintersect(A,B,C)); printf("\n");
  printf("setxor: "); printset(setxor(A,B,C)); printf("\n");
  printf("card(A): %i\n", card(A));
  
  freeset(A);
  freeset(B);
  freeset(C);
}
*/

LList  Lpush(LList L, Node *n)
{
  LNode *tmp = (LNode*)malloc(sizeof(LNode));
  tmp->next  = L;
  tmp->pNode = n;
  return tmp;
}

Node* Lpop(LList *L)
{
  if((*L) == NULL) return NULL;
  Node *n = (*L)->pNode;
  LNode *tmp = (*L);
  (*L) = (*L)->next;
  free(tmp);
  return n;
}

void Lfreelist(LList L)
{
  LNode *next;
  while(L != NULL)
	{
	  next = L->next;
	  free(L);
	  L = next;
	};
}

int Llistlen(LList L)
{
  int i=0;
  LNode *n;
  for(n=L;n;n=n->next)
	i++;
  return i;
}
