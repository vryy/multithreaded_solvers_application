
/*
 * -- SuperLU MT routine (version 2.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
 * and Xerox Palo Alto Research Center.
 * September 10, 2007
 *
 */
#include <math.h>
#include "pcsp_defs.h"

void
cCreate_CompCol_Matrix(SuperMatrix *A, int m, int n, int nnz, complex *nzval,
		      int *rowind, int *colptr,
		      Stype_t stype, Dtype_t dtype, Mtype_t mtype)
{
    NCformat *Astore;

    A->Stype = stype;
    A->Dtype = dtype;
    A->Mtype = mtype;
    A->nrow = m;
    A->ncol = n;
    A->Store = (void *) SUPERLU_MALLOC( sizeof(NCformat) );
    Astore = (NCformat *) A->Store;
    Astore->nnz = nnz;
    Astore->nzval = nzval;
    Astore->rowind = rowind;
    Astore->colptr = colptr;
}

void
cCreate_CompCol_Permuted(SuperMatrix *A, int m, int n, int nnz, complex *nzval,
			 int *rowind, int *colbeg, int *colend,
			 Stype_t stype, Dtype_t dtype, Mtype_t mtype)
{
    NCPformat *Astore;

    A->Stype = stype;
    A->Dtype = dtype;
    A->Mtype = mtype;
    A->nrow = m;
    A->ncol = n;
    A->Store = (void *) SUPERLU_MALLOC( sizeof(NCPformat) );
    Astore = (NCPformat *) A->Store;
    Astore->nnz = nnz;
    Astore->nzval = nzval;
    Astore->rowind = rowind;
    Astore->colbeg = colbeg;
    Astore->colend = colend;
}
/*
 * Convert a row compressed storage into a column compressed storage.
 */
void
cCompRow_to_CompCol(int m, int n, int nnz, 
		    complex *a, int *colind, int *rowptr,
		    complex **at, int **rowind, int **colptr)
{
    register int i, j, col, relpos;
    int *marker;

    /* Allocate storage for another copy of the matrix. */
    *at = (complex *) complexMalloc(nnz);
    *rowind = (int *) intMalloc(nnz);
    *colptr = (int *) intMalloc(n+1);
    marker = (int *) intCalloc(n);
    
    /* Get counts of each column of A, and set up column pointers */
    for (i = 0; i < m; ++i)
	for (j = rowptr[i]; j < rowptr[i+1]; ++j) ++marker[colind[j]];
    (*colptr)[0] = 0;
    for (j = 0; j < n; ++j) {
	(*colptr)[j+1] = (*colptr)[j] + marker[j];
	marker[j] = (*colptr)[j];
    }

    /* Transfer the matrix into the compressed column storage. */
    for (i = 0; i < m; ++i) {
	for (j = rowptr[i]; j < rowptr[i+1]; ++j) {
	    col = colind[j];
	    relpos = marker[col];
	    (*rowind)[relpos] = i;
	    (*at)[relpos] = a[j];
	    ++marker[col];
	}
    }

    SUPERLU_FREE(marker);
}


/* Copy matrix A into matrix B. */
void
cCopy_CompCol_Matrix(SuperMatrix *A, SuperMatrix *B)
{
    NCformat *Astore, *Bstore;
    int      ncol, nnz, i;

    B->Stype = A->Stype;
    B->Dtype = A->Dtype;
    B->Mtype = A->Mtype;
    B->nrow  = A->nrow;;
    B->ncol  = ncol = A->ncol;
    Astore   = (NCformat *) A->Store;
    Bstore   = (NCformat *) B->Store;
    Bstore->nnz = nnz = Astore->nnz;
    for (i = 0; i < nnz; ++i)
	((complex *)Bstore->nzval)[i] = ((complex *)Astore->nzval)[i];
    for (i = 0; i < nnz; ++i) Bstore->rowind[i] = Astore->rowind[i];
    for (i = 0; i <= ncol; ++i) Bstore->colptr[i] = Astore->colptr[i];
}


int cPrint_CompCol_Matrix(SuperMatrix *A)
{
    NCformat     *Astore;
    register int i;
    float       *dp;
    
    printf("\nCompCol matrix: ");
    printf("Stype %d, Dtype %d, Mtype %d\n", A->Stype,A->Dtype,A->Mtype);
    Astore = (NCformat *) A->Store;
    dp = (float *) Astore->nzval;
    printf("nrow %d, ncol %d, nnz %d\n", A->nrow,A->ncol,Astore->nnz);
    printf("\nnzval: ");
    for (i = 0; i < 2*Astore->nnz; ++i) printf("%f  ", dp[i]);
    printf("\nrowind: ");
    for (i = 0; i < Astore->nnz; ++i) printf("%d  ", Astore->rowind[i]);
    printf("\ncolptr: ");
    for (i = 0; i <= A->ncol; ++i) printf("%d  ", Astore->colptr[i]);
    printf("\nend CompCol matrix.\n");

    return 0;
}

int cPrint_Dense_Matrix(SuperMatrix *A)
{
    DNformat     *Astore;
    register int i;
    float       *dp;
    
    printf("\nDense matrix: ");
    printf("Stype %d, Dtype %d, Mtype %d\n", A->Stype,A->Dtype,A->Mtype);
    Astore = (DNformat *) A->Store;
    dp = (float *) Astore->nzval;
    printf("nrow %d, ncol %d, lda %d\n", A->nrow,A->ncol,Astore->lda);
    printf("\nnzval: ");
    for (i = 0; i < 2*A->nrow; ++i) printf("%f  ", dp[i]);
    printf("\nend Dense matrix.\n");

    return 0;
}

void
cCreate_Dense_Matrix(SuperMatrix *X, int m, int n, complex *x, int ldx,
		    Stype_t stype, Dtype_t dtype, Mtype_t mtype)
{
    DNformat    *Xstore;
    
    X->Stype = stype;
    X->Dtype = dtype;
    X->Mtype = mtype;
    X->nrow = m;
    X->ncol = n;
    X->Store = (void *) SUPERLU_MALLOC( sizeof(DNformat) );
    Xstore = (DNformat *) X->Store;
    Xstore->lda = ldx;
    Xstore->nzval = (complex *) x;
}

void
cCopy_Dense_Matrix(int M, int N, complex *X, int ldx, complex *Y, int ldy)
{
/*
 *
 *  Purpose
 *  =======
 *
 *  Copies a two-dimensional matrix X to another matrix Y.
 */
    int    i, j;
    
    for (j = 0; j < N; ++j)
        for (i = 0; i < M; ++i)
            Y[i + j*ldy] = X[i + j*ldx];
}

void
cCreate_SuperNode_Matrix(SuperMatrix *L, int m, int n, int nnz, complex *nzval,
			int *nzval_colptr, int *rowind, int *rowind_colptr,
			int *col_to_sup, int *sup_to_col,
			Stype_t stype, Dtype_t dtype, Mtype_t mtype)
{
    SCformat *Lstore;

    L->Stype = stype;
    L->Dtype = dtype;
    L->Mtype = mtype;
    L->nrow = m;
    L->ncol = n;
    L->Store = (void *) SUPERLU_MALLOC( sizeof(SCformat) );
    Lstore = L->Store;
    Lstore->nnz = nnz;
    Lstore->nsuper = col_to_sup[n];
    Lstore->nzval = nzval;
    Lstore->nzval_colptr = nzval_colptr;
    Lstore->rowind = rowind;
    Lstore->rowind_colptr = rowind_colptr;
    Lstore->col_to_sup = col_to_sup;
    Lstore->sup_to_col = sup_to_col;

}

void
cCreate_SuperNode_Permuted(SuperMatrix *L, int m, int n, int nnz,
			   complex *nzval, 
			   int *nzval_colbeg, int *nzval_colend,
			   int *rowind, int *rowind_colbeg, int *rowind_colend,
			   int *col_to_sup, 
			   int *sup_to_colbeg, int *sup_to_colend,
			   Stype_t stype, Dtype_t dtype, Mtype_t mtype)
{
    SCPformat *Lstore;

    L->Stype = stype;
    L->Dtype = dtype;
    L->Mtype = mtype;
    L->nrow = m;
    L->ncol = n;
    L->Store = (void *) SUPERLU_MALLOC( sizeof(SCPformat) );
    Lstore = L->Store;
    Lstore->nnz = nnz;
    Lstore->nsuper = col_to_sup[n];
    Lstore->nzval = nzval;
    Lstore->nzval_colbeg = nzval_colbeg;
    Lstore->nzval_colend = nzval_colend;
    Lstore->rowind = rowind;
    Lstore->rowind_colbeg = rowind_colbeg;
    Lstore->rowind_colend = rowind_colend;
    Lstore->col_to_sup = col_to_sup;
    Lstore->sup_to_colbeg = sup_to_colbeg;
    Lstore->sup_to_colend = sup_to_colend;

}


/*
 * Diagnostic print of column "jcol" in the U/L factor.
 */
void
cprint_lu_col(int pnum, char *msg, int pcol, int jcol, int w, int pivrow,
	      int *xprune, GlobalLU_t *Glu)
{
    int     i, k, fsupc;
    int     *xsup, *supno;
    int     *xlsub, *xlsub_end, *lsub;
    complex  *lusup;
    int     *xlusup, *xlusup_end;

    xsup    = Glu->xsup;
    supno   = Glu->supno;
    lsub    = Glu->lsub;
    xlsub   = Glu->xlsub;
    xlsub_end = Glu->xlsub_end;
    lusup   = Glu->lusup;
    xlusup  = Glu->xlusup;
    xlusup_end = Glu->xlusup_end;
    
    printf("(%d)%s fstcol %d,col %d,w %d: pivrow %d, supno %d, xprune %d\n", 
	   pnum, msg, pcol, jcol, w, pivrow, supno[jcol], xprune[jcol]);
    
    printf("(%d)\tU-col: xusub %d - %d\n",
	   pnum, Glu->xusub[jcol], Glu->xusub_end[jcol]);
    for (i = Glu->xusub[jcol]; i < Glu->xusub_end[jcol]; i++)
	printf("(%d)\t%d\t%8e\n", pnum, Glu->usub[i], Glu->ucol[i]);
    fsupc = xsup[supno[jcol]];
    k = xlusup[jcol];
    printf("(%d)\tL-col in s-node: xlsub %d - %d, xlusup %d - %d\n",
	   pnum, xlsub[fsupc],xlsub_end[fsupc],xlusup[jcol],xlusup_end[jcol]);
    for (i = xlsub[fsupc]; i < xlsub_end[fsupc]; ++i)
	printf("(%d)\t%d\t%.8e\n", pnum, lsub[i], lusup[k++]);

    fflush(stdout);
}

/*Dan fix above printf 's*/

/*
 * Check whether vec[*] == 0. For the two vectors dense[*] and tempv[*],
 * this invariant should be mantained before and after calling some
 * numeric update routines, such as "panel_bmod" and "column_bmod". 
 */
void
ccheck_zero_vec(int pnum, char *msg, int n, complex *vec)
{
    register int i, nonzero;

    nonzero = FALSE;
    for (i = 0; i < n; ++i) {
        if ((vec[i].r != 0.0) || (vec[i].i != 0.0))
        {
            printf("(%d) vec[%d] = %.10e; should be zero!\n",
                   pnum, i, vec[i]);
            nonzero = TRUE;
        }
    }
    if ( nonzero ) {
	printf("(%d) %s\n", pnum, msg);
	SUPERLU_ABORT("Not a zero vector.");
    } else {
        printf(".. Normal exit ccheck_zero_vec() ..\n");
    }
}


void
cGenXtrue(int n, int nrhs, complex *x, int ldx)
{
    int  i, j;
    for (j = 0; j < nrhs; ++j) {
	for (i = 0; i < n; ++i) {
            x[i + j*ldx].r = 1.0;
            x[i + j*ldx].i = 0.0;
        }
    }
}

/*
 * Let rhs[i] = sum of i-th row of A, so the solution vector is all 1's
 */
void
cFillRHS(trans_t trans, int nrhs, complex *x, int ldx, SuperMatrix *A, SuperMatrix *B)
{
    NCformat *Astore;
    complex   *Aval;
    DNformat *Bstore;
    complex   *rhs;
    complex one = {1.0, 0.0};
    complex zero = {0.0, 0.0};
    int      ldc;
    char     trans_c[1];

    Astore = A->Store;
    Aval   = (complex *) Astore->nzval;
    Bstore = B->Store;
    rhs    = Bstore->nzval;
    ldc    = Bstore->lda;
    
    if ( trans == NOTRANS ) *trans_c = 'N';
    else *trans_c = 'T';
    
    sp_cgemm(trans_c, A->nrow, nrhs, A->ncol, one, A,
	     x, ldx, zero, rhs, ldc);
}

/* 
 * Fills a double precision array with a given value.
 */
void 
cfill(complex *a, int alen, complex dval)
{
    register int i;
    for (i = 0; i < alen; i++) a[i] = dval;
}



/* 
 * Check the inf-norm of the error vector 
 */
void cinf_norm_error(int nrhs, SuperMatrix *X, complex *xtrue)
{
    DNformat *Xstore;
    float err, xnorm;
    complex *Xmat, *soln_work;
    complex temp;
    int i, j;

    Xstore = X->Store;
    Xmat = Xstore->nzval;

    for (j = 0; j < nrhs; j++) {
      soln_work = &Xmat[j*Xstore->lda];
      err = xnorm = 0.0;
      for (i = 0; i < X->nrow; i++) {
        c_sub(&temp, &soln_work[i], &xtrue[i]);
        err = SUPERLU_MAX(err, c_abs(&temp));
        xnorm = SUPERLU_MAX(xnorm, c_abs(&soln_work[i]));
      }
      err = err / xnorm;
      printf("||X - Xtrue||/||X|| = %e\n", err);
    }
}



/* Print performance of the code. */
void
cPrintPerf(SuperMatrix *L, SuperMatrix *U, superlu_memusage_t *superlu_memusage, 
	float rpg, float rcond, float *ferr,
	float *berr, char *equed, Gstat_t *Gstat)
{
    SCPformat *Lstore;
    NCPformat *Ustore;
    double   *utime;
    flops_t  *ops;
    
    utime = Gstat->utime;
    ops   = Gstat->ops;
    
    if ( utime[FACT] != 0. )
	printf("Factor flops = %e\tMflops = %8.2f\n", ops[FACT],
	       ops[FACT]*1e-6/utime[FACT]);
    printf("Identify relaxed snodes	= %8.2f\n", utime[RELAX]);
    if ( utime[SOLVE] != 0. )
	printf("Solve flops = %.0f, Mflops = %8.2f\n", ops[SOLVE],
	       ops[SOLVE]*1e-6/utime[SOLVE]);
    
    Lstore = (SCPformat *) L->Store;
    Ustore = (NCPformat *) U->Store;
    printf("\t#NZ in factor L = %d\n", Lstore->nnz);
    printf("\t#NZ in factor U = %d\n", Ustore->nnz);
    printf("\t#NZ in L+U = %d\n", Lstore->nnz + Ustore->nnz - L->ncol);
	
    printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n",
	   superlu_memusage->for_lu/1e6, superlu_memusage->total_needed/1e6,
	   superlu_memusage->expansions);
	
    printf("\tFactor\tMflops\tSolve\tMflops\tEtree\tEquil\tRcond\tRefine\n");
    printf("PERF:%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f\n",
	   utime[FACT], ops[FACT]*1e-6/utime[FACT],
	   utime[SOLVE], ops[SOLVE]*1e-6/utime[SOLVE],
	   utime[ETREE], utime[EQUIL], utime[RCOND], utime[REFINE]);
    
    printf("\tRpg\t\tRcond\t\tFerr\t\tBerr\t\tEquil?\n");
    printf("NUM:\t%e\t%e\t%e\t%e\t%s\n",
	   rpg, rcond, ferr[0], berr[0], equed);
    
#if 0

    printf("\tTRSV (total%%)\tGEMV (total%%)\tfloat_time%%\tmax_n\tmax_m\tmin_n\tmin_m\tavg_n\tavg_m\n");
    printf("BLAS:\t%.0f  %.2f\t%.0f  %.2f\t%.2f\t\t%d\t%d\t%d\t%d\t%.0f\t%.0f\n",
	   ops[TRSV], ops[TRSV]/ops[FACT], ops[GEMV], ops[GEMV]/ops[FACT],
	   utime[FLOAT]/utime[FACT],
	   max_blas_n, max_gemv_m, min_blas_n, min_gemv_m,
	   (float)sum_blas_n/num_blas, (float)sum_gemv_m/num_blas);
    printf("\tRCOND\tREFINE\tFERR\n");
    printf("SOLVES:\t%d\t%d\t%d\n", no_solves[RCOND],
	   no_solves[REFINE], no_solves[FERR]);
    
    flops_dist_for_matlab();

#endif
    
}


int print_complex_vec(char *what, int n, int *ind, complex *vec)
{
    int i;
    printf("%s: n %d\n", what, n);
    for (i = 0; i < n; ++i) printf("%d\t%f%f\n", ind[i], vec[i].r, vec[i].i);
    return 0;
}

