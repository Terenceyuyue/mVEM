#include "mex.h"

#ifndef _FOURIER__H
#define _FOURIER__H

// Set functions
class CSet {
public:
  CSet(int SetSize) {LN = SetSize;};

  bool setisempty(int *A);
  int* newset();
  void freeset(int *A);
  int* setcopy(int *C, int *A);
  int* setdiff(int *A, int *B, int *C);
  int* setunion(int *A, int *B, int *C);
  int* setintersect(int *A, int *B, int *C);
  int* setxor(int *A, int *B, int *C);
  int* setcomp(int *A, int *C);
  int  card(int *A);
  void printset(int *A);
  int* ind(int *A, int *C);

private:
  int LN;
};

// List structures
struct Node
{
  double *a;
  double  b;
  int    *H, *E, *I;
  Node   *prev,*next;
};

typedef Node* List;

// List functions
Node* newnode();
void  freenode(Node *N);
void  freelist(List L);
List  push(List L, Node *n);
Node* pop(List L);
List  remove(List L, Node *n);
int   listlen(List L);

// Container list for Node's
struct LNode
{
  Node *pNode;
  LNode *next;
};
typedef LNode *LList;
LList  Lpush(LList L, Node *n);
Node* Lpop(LList *L);
void   Lfreelist(LList L);
int    Llistlen(LList L);

// Fourier elim function
mxArray* fourier(double *H, int lnH, int *ax, double TOL, double QUASI_TOL);

void PrintConstraints(List S);
List elim_strong_redundancies(List S);

#endif
