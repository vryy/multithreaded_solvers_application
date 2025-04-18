/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*             ********   ***                                 SparseLib++    */
/*          *******  **  ***       ***      ***                              */
/*           *****      ***     ******** ********                            */
/*            *****    ***     ******** ********              R. Pozo        */
/*       **  *******  ***   **   ***      ***                 K. Remington   */
/*        ********   ********                                 A. Lumsdaine   */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                                                                           */
/*                                                                           */
/*                     SparseLib++ : Sparse Matrix Library                   */
/*                                                                           */
/*               National Institute of Standards and Technology              */
/*                        University of Notre Dame                           */
/*              Authors: R. Pozo, K. Remington, A. Lumsdaine                 */
/*                                                                           */
/*                                 NOTICE                                    */
/*                                                                           */
/* Permission to use, copy, modify, and distribute this software and         */
/* its documentation for any purpose and without fee is hereby granted       */
/* provided that the above notice appear in all copies and supporting        */
/* documentation.                                                            */
/*                                                                           */
/* Neither the Institutions (National Institute of Standards and Technology, */
/* University of Notre Dame) nor the Authors make any representations about  */
/* the suitability of this software for any purpose.  This software is       */
/* provided ``as is'' without expressed or implied warranty.                 */
/*                                                                           */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "boost/progress.hpp"
#include "ilutpre_double.h"
#include "pre.h"

CompRow_ILUtPreconditioner_double::
CompRow_ILUtPreconditioner_double(const CompRow_Mat_double &A, int lfil, double tol)
{
  // Copy
  dim_[0] = A.dim(0);
  dim_[1] = A.dim(1);

/*----------------------------------------------------------------------------
 * ILUT preconditioner
 * incomplete LU factorization with dual truncation mechanism
 * NOTE : no pivoting implemented as yet in GE for diagonal elements
 *----------------------------------------------------------------------------
 * Parameters
 *----------------------------------------------------------------------------
 * on entry:
 * =========
 * lfil     = integer. The fill-in parameter. Each column of L and
 *            each column of U will have a maximum of lfil elements.
 *            WARNING: THE MEANING OF LFIL HAS CHANGED WITH RESPECT TO
 *            EARLIER VERSIONS. 
 *            lfil must be .ge. 0.
 * tol      = real*8. Sets the threshold for dropping small terms in the
 *            factorization. See below for details on dropping strategy.
 * fp       = file pointer for error log ( might be stdout )
 *
 * on return:
 * ==========
 * ierr     = return value.
 *            ierr  = 0   --> successful return.
 *            ierr  = -1  --> Illegal value for lfil
 *            ierr  = -2  --> zero diagonal or zero col encountered
 * lu->n    = dimension of the matrix
 *   ->L    = L part -- stored in SpaFmt format
 *   ->D    = Diagonals
 *   ->U    = U part -- stored in SpaFmt format
 *----------------------------------------------------------------------------
 * Notes:
 * ======
 * All the diagonals of the input matrix must not be zero
 *----------------------------------------------------------------------------
 * Dual drop-off strategy works as follows. 
 *
 * 1) Theresholding in L and U as set by tol. Any element whose size
 *    is less than some tolerance (relative to the norm of current
 *    row in u) is dropped.
 *
 * 2) Keeping only the largest lfil elements in the i-th column of L
 *    and the largest lfil elements in the i-th column of U.
 *
 * Flexibility: one can use tol=0 to get a strategy based on keeping the
 * largest elements in each column of L and U. Taking tol .ne. 0 but lfil=n
 * will give the usual threshold strategy (however, fill-in is then
 * impredictible).
 *--------------------------------------------------------------------------*/
  FILE* fp = stdout;
  int n = dim_[0]; 
  int len, lenu, lenl;
  int nzcount, *ja, *jbuf, *iw, i, j, k;
  int col, jpos, jrow, upos;
  double t, tnorm, tolnorm, fact, lxu, *wn, *ma, *w;
  csptr L, U;
  double *D;

  if( lfil < 0 ) {
    fprintf( fp, "ilut: Illegal value for lfil.\n" );
    exit(0);
  }

  lu = (iluptr) Malloc( sizeof(ILUSpar), "main" );
  setupILU( lu, n );
  L = lu->L;
  U = lu->U;
  D = lu->D;

  iw   = (int*) Malloc( n * sizeof(int), "ilut" );
  jbuf = (int*) Malloc( n * sizeof(int), "ilut" );
  wn   = (double*) Malloc( n * sizeof(double), "ilut" );
  w    = (double*) Malloc( n * sizeof(double), "ilut" );

  /* set indicator array jw to -1 */
  for( i = 0; i < n; ++i ) iw[i] = -1;

  /* beginning of main loop */
  boost::progress_display show_progress( n );
  for( i = 0; i < n; ++i ) {
    // create row-wise data structure from matrix A
    nzcount = A.row_ptr(i+1) - A.row_ptr(i);    //number of nonzeros of row i
    int* col_ind = (int*) malloc(nzcount * sizeof(int));
    for(j = 0; j < nzcount; ++j) col_ind[j] = A.col_ind(A.row_ptr(i) + j);
    double* val = (double*) malloc(nzcount * sizeof(double));
    for(j = 0; j < nzcount; ++j) val[j] = A.val(A.row_ptr(i) + j);
    
    // start ilut adapted from ilut.c
    ja = col_ind;          //column indices of row i
    ma = val;              //values of row i
    tnorm = 0;
    for( j = 0; j < nzcount; ++j ) {
      tnorm += fabs( ma[j] );       //absolute norm of row i
    }
    if( tnorm == 0.0 ) {
      fprintf( fp, "ilut: zero row encountered.\n" );
      exit(1);
    }
    tnorm /= (double)nzcount;       //absolute norm = sum(j)(abs(a_ij)) / nz
    tolnorm = tol * tnorm;

    /* unpack L-part and U-part of column of A in arrays w */
    lenu = 0;
    lenl = 0;
    jbuf[i] = i;
    w[i] = 0;
    iw[i] = i;
    for( j = 0; j < nzcount; ++j ) {
      col = ja[j];
      t = ma[j];
      if( col < i ) {
        iw[col] = lenl;
        jbuf[lenl] = col;
        w[lenl] = t;
        ++lenl;
      } else if( col == i ) {
        w[i] = t;
      } else {
        ++lenu;
        jpos = i + lenu;
        iw[col] = jpos;
        jbuf[jpos] = col;
        w[jpos] = t;
      }
    }

    j = -1;
    len = 0;
    /* eliminate previous rows */
    while( ++j < lenl ) {
/*----------------------------------------------------------------------------
 *  in order to do the elimination in the correct order we must select the
 *  smallest column index among jbuf[k], k = j+1, ..., lenl
 *--------------------------------------------------------------------------*/
      jrow = jbuf[j];
      jpos = j;
      /* determine smallest column index */
      for( k = j + 1; k < lenl; ++k ) {
        if( jbuf[k] < jrow ) {
          jrow = jbuf[k];
          jpos = k;
        }
      }
      if( jpos != j ) {
        col = jbuf[j];
        jbuf[j] = jbuf[jpos];
        jbuf[jpos] = col;
        iw[jrow] = j;
        iw[col]  = jpos;
        t = w[j];
        w[j] = w[jpos];
        w[jpos] = t;
      }

      /* get the multiplier */
      fact = w[j] * D[jrow];
      w[j] = fact;
      /* zero out element in row by resetting iw(n+jrow) to -1 */
      iw[jrow] = -1;

      /* combine current row and row jrow */
      nzcount = U->nzcount[jrow];
      ja = U->ja[jrow];
      ma = U->ma[jrow];
      for( k = 0; k < nzcount; ++k ) {
        col = ja[k];
        jpos = iw[col];
        lxu = - fact * ma[k];
        /* if fill-in element is small then disregard */
        if( fabs( lxu ) < tolnorm && jpos == -1 ) continue;

        if( col < i ) {
          /* dealing with lower part */
          if( jpos == -1 ) {
            /* this is a fill-in element */
            jbuf[lenl] = col;
            iw[col] = lenl;
            w[lenl] = lxu;
            ++lenl;
          } else {
            w[jpos] += lxu;
          }
        } else {
          /* dealing with upper part */
//          if( jpos == -1 ) {
          if( jpos == -1 && fabs(lxu) > tolnorm) {
            /* this is a fill-in element */
            ++lenu;
            upos = i + lenu;
            jbuf[upos] = col;
            iw[col] = upos;
            w[upos] = lxu;
          } else {
            w[jpos] += lxu;
          }
        }
      }
    }

    /* restore iw */
    iw[i] = -1;
    for( j = 0; j < lenu; ++j ) {
      iw[jbuf[i+j+1]] = -1;
    }

/*---------- case when diagonal is zero */
    if( w[i] == 0.0 ) {
      fprintf( fp, "zero diagonal encountered.\n" );
      for( j = i; j < n; ++j ) {
        L->ja[j] = NULL; 
        L->ma[j] = NULL;
        U->ja[j] = NULL; 
        U->ma[j] = NULL;
      }
      exit(1);
    }
/*-----------Update diagonal */    
    D[i] = 1 / w[i];

    /* update L-matrix */
//    len = min( lenl, lfil );
    len = lenl < lfil ? lenl : lfil;
    for( j = 0; j < lenl; ++j ) {
      wn[j] = fabs( w[j] );
      iw[j] = j;
    }
    qsplit( wn, iw, &lenl, &len );
    L->nzcount[i] = len;
    if( len > 0 ) {
      ja = L->ja[i] = (int*) Malloc( len * sizeof(int), "ilut" );
      ma = L->ma[i] = (double*) Malloc( len * sizeof(double), "ilut" );
    }
    for( j = 0; j < len; ++j ) {
      jpos = iw[j];
      ja[j] = jbuf[jpos];
      ma[j] = w[jpos];
    }
    for( j = 0; j < lenl; ++j ) iw[j] = -1;

    /* update U-matrix */
//    len = min( lenu, lfil );
    len = lenu < lfil ? lenu : lfil;
    for( j = 0; j < lenu; ++j ) {
      wn[j] = fabs( w[i+j+1] );
      iw[j] = i+j+1;
    }
    qsplit( wn, iw, &lenu, &len );
    U->nzcount[i] = len;
    if( len > 0 ) {
      ja = U->ja[i] = (int*) Malloc( len * sizeof(int), "ilut" );
      ma = U->ma[i] = (double*) Malloc( len * sizeof(double), "ilut" );
    }
    for( j = 0; j < len; ++j ) {
      jpos = iw[j];
      ja[j] = jbuf[jpos];
      ma[j] = w[jpos];
    }
    for( j = 0; j < lenu; ++j ) {
      iw[j] = -1;
    }
    
    // release data to avoid memory leaking
    delete [] col_ind;
    delete [] val;
    
    ++show_progress;
  }

  free( iw );
  free( jbuf );
  free( wn );
}


VECTOR_double
CompRow_ILUtPreconditioner_double::solve(const VECTOR_double &y) const 
{
  int M = y.size();
  VECTOR_double x(M);
  
/*----------------------------------------------------------------------
 *    performs a forward followed by a backward solve
 *    for LU matrix as produced by ilut
 *    y  = right-hand-side
 *    x  = solution on return
 *    lu = LU matrix as produced by ilut.
 *--------------------------------------------------------------------*/
    int n = lu->n, i, j, nzcount, *ja;
    double *D, *ma;
    csptr L, U;

    L = lu->L;
    U = lu->U;
    D = lu->D;

    /* Block L solve */
    for( i = 0; i < n; ++i ) {
        x(i) = y(i);
        nzcount = L->nzcount[i];
        ja = L->ja[i];
        ma = L->ma[i];
        for( j = 0; j < nzcount; ++j ) {
            x(i) -= x(ja[j]) * ma[j];
        }
    }
    /* Block -- U solve */
    for( i = n-1; i >= 0; --i ) {
        nzcount = U->nzcount[i];
        ja = U->ja[i];
        ma = U->ma[i];
        for( j = 0; j < nzcount; ++j ) {
            x(i) -= x(ja[j]) * ma[j];
        }
        x(i) *= D[i];
    }

  return x;
} 


VECTOR_double
CompRow_ILUtPreconditioner_double::trans_solve(const VECTOR_double &x) const 
{
  printf("trans_solve is not yet supported\n");
  exit(0);
}

