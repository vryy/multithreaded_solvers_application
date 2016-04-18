#include "pre.h"

void errexit( char *f_str, ... )
{
  va_list argp;
  char out1[256], out2[256];

  va_start(argp, f_str);
  vsprintf(out1, f_str, argp);
  va_end(argp);

  sprintf(out2, "Error! %s\n", out1);

  fprintf(stdout, out2);
  fflush(stdout);

  exit( -1 );
}

void *Malloc( int nbytes, char *msg )
{
  void *ptr;

  if (nbytes == 0)
    return NULL;

  ptr = (void *)malloc(nbytes);
  if (ptr == NULL)
    errexit( "Not enough mem for %s. Requested size: %d bytes", msg, nbytes );

  return ptr;
}

int setupCS(csptr amat, int len, int job)
{
/*----------------------------------------------------------------------
| Initialize SpaFmt structs.
|----------------------------------------------------------------------
| on entry:
|==========
| ( amat )  =  Pointer to a SpaFmt struct.
|     len   =  size of matrix
|     job   =  0: pattern only
|              1: data and pattern
|
| On return:
|===========
|
|  amat->n
|      ->*nzcount
|      ->**ja
|      ->**ma
|
| integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|--------------------------------------------------------------------*/
   amat->n = len;
   amat->nzcount = (int *)Malloc( len*sizeof(int), "setupCS" );
   amat->ja = (int **) Malloc( len*sizeof(int *), "setupCS" );
   if( job == 1 ) 
       amat->ma = (double **) Malloc( len*sizeof(double *), "setupCS" );
   else
       amat->ma = NULL;
   return 0;
}
/*---------------------------------------------------------------------
|     end of setupCS
|--------------------------------------------------------------------*/

int cleanCS(csptr amat)
{
/*----------------------------------------------------------------------
| Free up memory allocated for SpaFmt structs.
|----------------------------------------------------------------------
| on entry:
|==========
| ( amat )  =  Pointer to a SpaFmt struct.
|--------------------------------------------------------------------*/
   /*   */
  int i;
  if (amat == NULL) return 0;
  if (amat->n < 1) return 0;
  for (i=0; i<amat->n; i++) {
    if (amat->nzcount[i] > 0) {
      if( amat->ma ) free(amat->ma[i]);
      free(amat->ja[i]);
    }
  }    
  if (amat->ma) free(amat->ma);
  free(amat->ja);
  free(amat->nzcount);
  free(amat);
  return 0;
}
/*---------------------------------------------------------------------
|     end of cleanCS
|--------------------------------------------------------------------*/

int setupILU( iluptr lu, int n )
{
/*----------------------------------------------------------------------
| Initialize ILUSpar structs.
|----------------------------------------------------------------------
| on entry:
|==========
|   ( lu )  =  Pointer to a ILUSpar struct.
|       n   =  size of matrix
|
| On return:
|===========
|
|    lu->n
|      ->L     L matrix, SpaFmt format
|      ->D     Diagonals
|      ->U     U matrix, SpaFmt format
|      ->work  working buffer of length n
|      ->bf    buffer
|
| integer value returned:
|             0   --> successful return.
|            -1   --> memory allocation error.
|--------------------------------------------------------------------*/
    lu->n  = n;
    lu->D = (double *)Malloc( sizeof(double) * n, "setupILU" );
    lu->L = (csptr)Malloc( sizeof(SparMat), "setupILU" );
    setupCS( lu->L, n, 1 );
    lu->U = (csptr)Malloc( sizeof(SparMat), "setupILU" );
    setupCS( lu->U, n, 1 );
    lu->work = (int *)Malloc( sizeof(int) * n, "setupILU" );
    return 0;
}
/*---------------------------------------------------------------------
|     end of setupILU
|--------------------------------------------------------------------*/

int cleanILU( iluptr lu )
{
/*----------------------------------------------------------------------
| Free up memory allocated for ILUSpar structs.
|----------------------------------------------------------------------
| on entry:
|==========
|   ( lu )  =  Pointer to a ILUSpar struct.
|--------------------------------------------------------------------*/
  if( NULL == lu ) return 0;
  if( lu->D ) {
    free( lu->D );
  }
  cleanCS( lu->L );
  cleanCS( lu->U );
  if( lu->work ) free( lu->work );
  free( lu );
  return 0;
}
/*---------------------------------------------------------------------
|     end of cleanILU
|--------------------------------------------------------------------*/

Preconditioner_double::Preconditioner_double() {}
Preconditioner_double::~Preconditioner_double() {}
VECTOR_double Preconditioner_double::solve (const VECTOR_double &x) const {return NULL;}
VECTOR_double Preconditioner_double::trans_solve (const VECTOR_double &x) const {return NULL;}

