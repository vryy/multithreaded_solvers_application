#ifndef __PRE_H__
#define __PRE_H__

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "vecdefs.h"
#include VECTOR_H

#define qsplit  qsplit_

void errexit( char *f_str, ... );

void *Malloc( int nbytes, char *msg );

typedef struct SpaFmt {
/*--------------------------------------------- 
| C-style CSR format - used internally
| for all matrices in CSR format 
|---------------------------------------------*/
  int n;
  int *nzcount;  /* length of each row */
  int **ja;      /* pointer-to-pointer to store column indices  */
  double **ma;   /* pointer-to-pointer to store nonzero entries */
} SparMat, *csptr;

typedef struct ILUfac {
    int n;
    csptr L;      /* L part elements                            */
    double *D;    /* diagonal elements                          */
    csptr U;      /* U part elements                            */
    int *work;    /* working buffer */
} ILUSpar, LDUmat, *iluptr;

class Preconditioner_double {
 public:
  Preconditioner_double();
  ~Preconditioner_double();
  virtual VECTOR_double solve (const VECTOR_double &x) const;
  virtual VECTOR_double trans_solve (const VECTOR_double &x) const;
};

int setupCS(csptr amat, int len, int job);
int cleanCS(csptr amat);
int setupILU(iluptr lu, int n);
int cleanILU(iluptr lu);

/* FORTRAN routines */
extern "C" void qsplit(double *a, int *ind, int *n, int *ncut);

#endif // __PRE_H__

