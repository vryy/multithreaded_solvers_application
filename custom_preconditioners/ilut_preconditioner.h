/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 25 Aug 2014 $
//   Revision:            $Revision: 1.2 $
//
//

#if !defined(KRATOS_MULTITHREADED_SOLVERS_APPLICATION_ILUT_PRECONDITIONER_H_INCLUDED)
#define KRATOS_MULTITHREADED_SOLVERS_APPLICATION_ILUT_PRECONDITIONER_H_INCLUDED




// System includes
#include <cstdio>


// External includes
#include <boost/smart_ptr.hpp>
#include <boost/numeric/ublas/vector.hpp>

// Project includes
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "utilities/progress.h"
#include "linear_solvers/preconditioner.h"
#include "linear_solvers/linear_solver.h"

/* FORTRAN routines */
#define qsplit qsplit_
extern "C" void qsplit(double *a, int *ind, int *n, int *ncut);

namespace Kratos
{


///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

///@name  Preconditioners
///@{

/// ITSOLPreconditioner class.
/**   */
template<class TSparseSpaceType, class TDenseSpaceType, class TModelPartType>
class ILUtPreconditioner : public Preconditioner<TSparseSpaceType, TDenseSpaceType, TModelPartType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ILUtPreconditioner
    KRATOS_CLASS_POINTER_DEFINITION (ILUtPreconditioner);

    typedef Preconditioner<TSparseSpaceType, TDenseSpaceType, TModelPartType> BaseType;

    typedef typename BaseType::SparseMatrixType SparseMatrixType;

    typedef typename BaseType::VectorType VectorType;

    typedef typename BaseType::DenseMatrixType DenseMatrixType;

    typedef typename BaseType::SizeType SizeType;

    typedef typename BaseType::IndexType IndexType;

    typedef typename BaseType::DataType DataType;

    typedef typename BaseType::ValueType ValueType;

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TModelPartType> LinearSolverType;

    typedef typename LinearSolverType::Pointer LinearSolverPointerType;

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
        double *D;   /* diagonal elements                          */
        csptr U;      /* U part elements                            */
        int *work;    /* working buffer */
    } ILUSpar, LDUmat, *iluptr;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ILUtPreconditioner(ValueType lfil, ValueType droptol)
    {
        if(lfil < 0.0 || lfil > 1.0)
            KRATOS_ERROR << "level of fill must be in [0.0 1.0]";
        mlfil = lfil;
        mdroptol = droptol;
        mlu = NULL;
    }


    /// Copy constructor.
    ILUtPreconditioner(const ILUtPreconditioner& Other)
    {
        mlfil = Other.mlfil;
        mdroptol = Other.mdroptol;
    }

    /// Destructor.
    ~ILUtPreconditioner() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    ILUtPreconditioner& operator=(const ILUtPreconditioner& Other)
    {
        mlfil = Other.mlfil;
        mdroptol = Other.mdroptol;
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    void Initialize(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
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
        int n = TSparseSpaceType::Size1(rA);
        int lfil = (int)(mlfil * n);
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

        mlu = (iluptr) Malloc( sizeof(ILUSpar), "mlu" );
        setupILU( mlu, n );
        L = mlu->L;
        U = mlu->U;
        D = mlu->D;

        iw   = (int*) Malloc( n * sizeof(int), "iw" );
        jbuf = (int*) Malloc( n * sizeof(int), "jbuf" );
        wn   = (double*) Malloc( n * sizeof(double), "wn" );
        w    = (double*) Malloc( n * sizeof(double), "w" );

        /* set indicator array jw to -1 */
        for( i = 0; i < n; ++i ) iw[i] = -1;

        /* beginning of main loop */
        Kratos::progress_display show_progress( n );
        for( i = 0; i < n; ++i ) {
            // create row-wise data structure from matrix A
            nzcount = rA.index1_data()[i + 1] - rA.index1_data()[i];    //number of nonzeros of row i
            int* col_ind = (int*) Malloc(nzcount * sizeof(int), "col_ind");
            for(j = 0; j < nzcount; ++j) col_ind[j] = rA.index2_data()[rA.index1_data()[i] + j];
            double* val = (double*) Malloc(nzcount * sizeof(double), "val");
            for(j = 0; j < nzcount; ++j) val[j] = rA.value_data()[rA.index1_data()[i] + j];

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
            tnorm /= (double) nzcount;       //absolute norm = sum(j)(abs(a_ij)) / nz
            tolnorm = mdroptol * tnorm;

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
    //                  if( jpos == -1 ) {
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
            len = (lenl < lfil ? lenl : lfil);
            for( j = 0; j < lenl; ++j ) {
                wn[j] = fabs( w[j] );
                iw[j] = j;
            }
            qsplit( wn, iw, &lenl, &len );
            L->nzcount[i] = len;
            if( len > 0 ) {
                ja = L->ja[i] = (int*) Malloc( len * sizeof(int), "jaL" );
                ma = L->ma[i] = (double*) Malloc( len * sizeof(double), "maL" );
            }
            for( j = 0; j < len; ++j ) {
                jpos = iw[j];
                ja[j] = jbuf[jpos];
                ma[j] = w[jpos];
            }
            for( j = 0; j < lenl; ++j ) iw[j] = -1;

            /* update U-matrix */
            len = (lenu < lfil ? lenu : lfil);
            for( j = 0; j < lenu; ++j ) {
                wn[j] = fabs( w[i+j+1] );
                iw[j] = i+j+1;
            }
            qsplit( wn, iw, &lenu, &len );
            U->nzcount[i] = len;
            if( len > 0 ) {
                ja = U->ja[i] = (int*) Malloc( len * sizeof(int), "jaU" );
                ma = U->ma[i] = (double*) Malloc( len * sizeof(double), "maU" );
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

    bool AdditionalPhysicalDataIsNeeded() override
    {
        return false;
    }

    void Mult(SparseMatrixType& rA, VectorType& rX, VectorType& rY) override
    {
        TSparseSpaceType::Mult(rA, rX, rY);
        ApplyLeft(rY);
    }

    /*
     * calculate preconditioned_X = A^{-1} * X;
     @param rX Unknows of preconditioner system
     */
    VectorType& ApplyLeft(VectorType& rX) override
    {
        /*----------------------------------------------------------------------
         *    performs a forward followed by a backward solve
         *    for LU matrix as produced by ilut
         *    y  = right-hand-side
         *    x  = solution on return
         *    lu = LU matrix as produced by ilut.
         *--------------------------------------------------------------------*/
        int n = mlu->n, i, j, nzcount, *ja;
        double *D, *ma;
        csptr L, U;

        double* x = new double[n];

        L = mlu->L;
        U = mlu->U;
        D = mlu->D;

        /* Block L solve */
        for( i = 0; i < n; ++i ) {
            x[i] = rX(i);
            nzcount = L->nzcount[i];
            ja = L->ja[i];
            ma = L->ma[i];
            for( j = 0; j < nzcount; ++j ) {
                x[i] -= x[ja[j]] * ma[j];
            }
        }
        /* Block -- U solve */
        for( i = n-1; i >= 0; --i ) {
            nzcount = U->nzcount[i];
            ja = U->ja[i];
            ma = U->ma[i];
            for( j = 0; j < nzcount; ++j ) {
                x[i] -= x[ja[j]] * ma[j];
            }
            x[i] *= D[i];
        }

        std::copy(x, x + n, rX.begin());
        delete [] x;
        return rX;
    }

    VectorType& Finalize(VectorType& rX) override
    {
        cleanILU( mlu );
        return rX;
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Return information about this object.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "ILUtPreconditioner, ";
        PrintData(buffer);
        return buffer.str();
    }

    /// Print information about this object.
    void  PrintInfo(std::ostream& OStream) const override
    {
        OStream << Info();
    }

    void PrintData(std::ostream& OStream) const override
    {
        OStream << "Level of fill  = " << mlfil;
        OStream << ", Drop tolerance = " << mdroptol;
    }

    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{
    iluptr mlu;    /* ilu preconditioner structure */
    double mlfil;
    ValueType mdroptol;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

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

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Unaccessible methods
    ///@{


    ///@}

}; // Class ILUtPreconditioner

///@}


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_MULTITHREADED_SOLVERS_APPLICATION_ILUT_PRECONDITIONER_H_INCLUDED defined
