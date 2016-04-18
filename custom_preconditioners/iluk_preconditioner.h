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
//   Date:                $Date: 27 Aug 2014 $
//   Revision:            $Revision: 1.1 $
//
//

#if !defined(KRATOS_MULTITHREADED_SOLVERS_APPLICATION_ILUK_PRECONDITIONER_H_INCLUDED)
#define KRATOS_MULTITHREADED_SOLVERS_APPLICATION_ILUK_PRECONDITIONER_H_INCLUDED




// System includes
#include <cstdio>


// External includes
#include <boost/smart_ptr.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/progress.hpp>

// Project includes
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "linear_solvers/preconditioner.h"
#include "linear_solvers/linear_solver.h"

#ifndef min
#define min(a,b) (((a)>(b))?(b):(a))
#endif
#ifndef max
#define max(a,b) (((a)>(b))?(a):(b))
#endif

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
template<class TSparseSpaceType, class TDenseSpaceType>
class  ILUkPreconditioner : public Preconditioner<TSparseSpaceType, TDenseSpaceType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of  ILUkPreconditioner
    KRATOS_CLASS_POINTER_DEFINITION ( ILUkPreconditioner);

    typedef Preconditioner<TSparseSpaceType, TDenseSpaceType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
    
    typedef typename TSparseSpaceType::MatrixPointerType SparseMatrixPointerType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType> LinearSolverType;
    
    typedef typename LinearSolverType::Pointer LinearSolverPointerType;

    typedef std::size_t  SizeType;
    
    typedef std::size_t  IndexType;
    
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
     ILUkPreconditioner(int lfil)
    {
        mlfil = lfil;
        mlu = NULL;
    }


    /// Copy constructor.
     ILUkPreconditioner(const ILUkPreconditioner& Other)
    {
        mlfil = Other.mlfil;
    }


    /// Destructor.
    virtual ~ILUkPreconditioner()
    {
    }



    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
     ILUkPreconditioner& operator=(const  ILUkPreconditioner& Other)
    {
        mlfil = Other.mlfil;
        return *this;
    }

    

    ///@}
    ///@name Operations
    ///@{

    virtual void Initialize(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
        /*----------------------------------------------------------------------------
         * ILUK preconditioner
         * incomplete LU factorization with level of fill dropping
         *----------------------------------------------------------------------------
         * Parameters
         *----------------------------------------------------------------------------
         * on entry:
         * =========
         * mlfil    = level of fill: all entries with level of fill > lofM are
         *            dropped. Setting lofM = 0 gives BILU(0).
         * lu       = pointer to a ILUKSpar struct -- see heads.h for details
         *            on format
         * fp       = file pointer for error log ( might be stderr )
         *
         * on return:
         * ==========
         * ierr     = return value.
         *            ierr  = 0   --> successful return.
         *            ierr  = -1  --> error in lofC
         *            ierr  = -2  --> zero diagonal found
         * lu->n    = dimension of the matrix
         *   ->L    = L part -- stored in SpaFmt format
         *   ->D    = Diagonals
         *   ->U    = U part -- stored in SpaFmt format
         *----------------------------------------------------------------------------
         * Notes:
         * ======
         * All the diagonals of the input matrix must not be zero
         *--------------------------------------------------------------------------*/
        FILE* fp = stderr;
        int n = TSparseSpaceType::Size1(rA);
        int ierr;
        int *jw, i, j, k, col, jpos, jrow, nzcount;
        csptr L, U;
        double *D;

        mlu = (iluptr) Malloc( sizeof(ILUSpar), "mlu" );
        setupILU( mlu, n );
        L = mlu->L;
        U = mlu->U;
        D = mlu->D;

        boost::progress_display show_progress(2 * n);
        /* symbolic factorization to calculate level of fill index arrays */
        if( ( ierr = lofC( mlfil, rA, mlu, fp, show_progress ) ) != 0 ) {
            fprintf( fp, "iluk error: lofC\n" );
            exit(0);
        }

        jw = mlu->work;
        /* set indicator array jw to -1 */
        for( j = 0; j < n; ++j ) jw[j] = -1;

        /* beginning of main loop */
        for( i = 0; i < n; ++i ) {
            /* set up the i-th row accroding to the nonzero information from
               symbolic factorization */
            mallocRow( mlu, i );

            /* setup array jw[], and initial i-th row */
            for( j = 0; j < L->nzcount[i]; ++j ) {  /* initialize L part   */
                col = L->ja[i][j];
                jw[col] = j;
                L->ma[i][j] = 0;
            }
            jw[i] = i;
            D[i] = 0; /* initialize diagonal */
            for( j = 0; j < U->nzcount[i]; ++j ) {  /* initialize U part   */
                col = U->ja[i][j];
                jw[col] = j;
                U->ma[i][j] = 0;
            }

            /* copy row from rA into lu */
            nzcount = rA.index1_data()[i + 1] - rA.index1_data()[i];    //number of nonzeros of row i
            for( j = 0; j < nzcount; ++j ) {
                col = rA.index2_data()[rA.index1_data()[i] + j]; //csmat->ja[i][j];
                jpos = jw[col];
                if( col < i )
                    L->ma[i][jpos] = rA.value_data()[rA.index1_data()[i] + j];// csmat->ma[i][j];
                else if( col == i )
                    D[i] = rA.value_data()[rA.index1_data()[i] + j]; //csmat->ma[i][j];
                else
                    U->ma[i][jpos] = rA.value_data()[rA.index1_data()[i] + j]; //csmat->ma[i][j];
            }

            /* eliminate previous rows */
            for( j = 0; j < L->nzcount[i]; ++j ) {
                jrow = L->ja[i][j];
                /* get the multiplier for row to be eliminated (jrow) */
                L->ma[i][j] *= D[jrow];

                /* combine current row and row jrow */
                for( k = 0; k < U->nzcount[jrow]; ++k ) {
                    col = U->ja[jrow][k];
                    jpos = jw[col];
                    if( jpos == -1 ) continue;
                    if( col < i )
                        L->ma[i][jpos] -= L->ma[i][j] * U->ma[jrow][k];
                    else if( col == i )
                        D[i] -= L->ma[i][j] * U->ma[jrow][k];
                    else
                        U->ma[i][jpos] -= L->ma[i][j] * U->ma[jrow][k];
                }
            }

            /* reset double-pointer to -1 ( U-part) */
            for( j = 0; j < L->nzcount[i]; ++j )
            {
                col = L->ja[i][j];
                jw[col] = -1;
            }
            jw[i] = -1;
            for( j = 0; j < U->nzcount[i]; ++j )
            {
                col = U->ja[i][j];
                jw[col] = -1;
            }

            if( D[i] == 0 ) {
                for( j = i+1; j < n; ++j ) {
                    L->ma[j] = NULL;
                    U->ma[j] = NULL;
                }
                fprintf( fp, "iluk fatal error: Zero diagonal found...\n" );
                exit(1);
            }
            D[i] = 1.0 / D[i];
            ++show_progress;
        }
    }
    
    
    virtual bool AdditionalPhysicalDataIsNeeded()
    {
        return false;
    }


    virtual void Mult(SparseMatrixType& rA, VectorType& rX, VectorType& rY)
    {
        TSparseSpaceType::Mult(rA, rX, rY);
        ApplyLeft(rY);
    }


    /*
     * calculate preconditioned_X = A^{-1} * X;
     @param rX Unknows of preconditioner system
     */
    virtual VectorType& ApplyLeft(VectorType& rX)
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

    virtual VectorType& Finalize(VectorType& rX)
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
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << " ILUkPreconditioner, ";
        PrintData(buffer);
        return buffer.str();
    }


    /// Print information about this object.
    virtual void  PrintInfo(std::ostream& OStream) const
    {
        OStream << Info();
    }


    virtual void PrintData(std::ostream& OStream) const
    {
        OStream << "Level of fill  = " << mlfil;
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
    int mlfil;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    int lofC( int lofM, SparseMatrixType& rA, iluptr lu, FILE *fp, boost::progress_display& progress )
    {
        /*--------------------------------------------------------------------
         * symbolic ilu factorization to calculate structure of ilu matrix
         * for specified level of fill
         *--------------------------------------------------------------------
         * on entry:
         * =========
         * lofM     = level of fill, lofM >= 0
         * rA       = compressed sparse row matrix structure
         * lu       = pointer to a ILUSpar struct -- see heads.h for details
         *            on format
         * fp       = file pointer for error log ( might be stderr )
         *--------------------------------------------------------------------
         * on return:
         * ==========
         * ierr     = return value.
         *            ierr  = 0   --> successful return.
         *            ierr != 0   --> error
         * lu->n    = dimension of the block matrix
         *   ->L    = L part -- stored in SpaFmt format, patterns only in lofC
         *   ->U    = U part -- stored in SpaFmt format, patterns only in lofC
         *------------------------------------------------------------------*/
        int n = TSparseSpaceType::Size1(rA);
        int *levls = NULL, *jbuf = NULL, *iw = lu->work;
        int **ulvl;  /*  stores lev-fils for U part of ILU factorization*/
        csptr L = lu->L, U = lu->U;
        /*--------------------------------------------------------------------
         * n        = number of rows or columns in matrix
         * inc      = integer, count of nonzero(fillin) element of each row
         *            after symbolic factorization
         * ju       = entry of U part of each row
         * lvl      = buffer to store levels of each row
         * jbuf     = buffer to store column index of each row
         * iw       = work array
         *------------------------------------------------------------------*/
        int i, j, k, col, ip, it, jpiv;
        int incl, incu, jmin, kmin, nzcount;
      
        levls  = (int*) Malloc( n*sizeof(int), "lofC: levls" );
        jbuf = (int*) Malloc( n*sizeof(int), "lofC: jbuf" ); 
        ulvl = (int**) Malloc( n*sizeof(int *), "lofC: ulvl" );

        /* initilize iw */
        for( j = 0; j < n; ++j ) iw[j] = -1;
        for( i = 0; i < n; ++i ) {
            incl = 0;
            incu = i;
            nzcount = rA.index1_data()[i + 1] - rA.index1_data()[i];    //number of nonzeros of row i
        /*-------------------- assign lof = 0 for matrix elements */
            for( j = 0; j < nzcount; ++j ) {
                col = rA.index2_data()[rA.index1_data()[i] + j]; //csmat->ja[i][j];
                if( col < i ) {
        /*-------------------- L-part  */
	                jbuf[incl] = col;
	                levls[incl] = 0;
	                iw[col] = incl++;
                } 
                else if (col > i) { 
        /*-------------------- U-part  */
	                jbuf[incu] = col;
	                levls[incu] = 0;
	                iw[col] = incu++;
                } 
            }
        /*-------------------- symbolic k,i,j Gaussian elimination  */ 
            jpiv = -1; 
            while (++jpiv < incl) {
                k = jbuf[jpiv] ; 
        /*-------------------- select leftmost pivot */
                kmin = k;
                jmin = jpiv; 
                for( j = jpiv + 1; j< incl; ++j) {
	                if( jbuf[j] < kmin ) {
	                    kmin = jbuf[j];
	                    jmin = j;
	                }
                }
        /*-------------------- swap  */  
                if( jmin != jpiv ) {
	                jbuf[jpiv] = kmin; 
	                jbuf[jmin] = k; 
	                iw[kmin] = jpiv;
	                iw[k] = jmin; 
	                j = levls[jpiv] ;
	                levls[jpiv] = levls[jmin];
	                levls[jmin] = j;
	                k = kmin; 
                }
        /*-------------------- symbolic linear combinaiton of rows  */
                for( j = 0; j < U->nzcount[k]; ++j ) {
	                col = U->ja[k][j];
	                it = ulvl[k][j]+levls[jpiv]+1 ; 
	                if( it > lofM ) continue; 
	                ip = iw[col];
	                if( ip == -1 ) {
	                    if( col < i) {
	                        jbuf[incl] = col;
	                        levls[incl] = it;
	                        iw[col] = incl++;
                        }
	                    else if( col > i ) {
	                        jbuf[incu] = col;
	                        levls[incu] = it;
	                        iw[col] = incu++;
	                    } 
                    }
                    else
	                    levls[ip] = min(levls[ip], it); 
                }
            }   /* end - while loop */
        /*-------------------- reset iw */
            for( j = 0; j < incl; ++j ) iw[jbuf[j]] = -1;
            for( j = i; j < incu; ++j ) iw[jbuf[j]] = -1;
        /*-------------------- copy L-part */ 
            L->nzcount[i] = incl;
            if(incl > 0 ) {
                L->ja[i] = (int*) Malloc( incl*sizeof(int), "lofC: L->ja[i]" );
                memcpy( L->ja[i], jbuf, sizeof(int)*incl);
            }
        /*-------------------- copy U - part        */ 
            k = incu-i; 
            U->nzcount[i] = k; 
            if( k > 0 ) {
                U->ja[i] = (int*) Malloc( sizeof(int)*k, "lofC: U->ja[i]" );
                memcpy( U->ja[i], jbuf+i, sizeof(int)*k );
        /*-------------------- update matrix of levels */
                ulvl[i] = (int*) Malloc( k*sizeof(int), "lofC: ulvl[i]" ); 
                memcpy( ulvl[i], levls+i, k*sizeof(int) );
            }
            ++progress;
        }
      
        /*-------------------- free temp space and leave --*/
        free(levls);
        free(jbuf);
        for(i = 0; i < n-1; ++i ) {
            if (U->nzcount[i]) free(ulvl[i]);
        }
        free(ulvl);

        return 0;
    }

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

    int mallocRow( iluptr lu, int nrow )
    {
        /*----------------------------------------------------------------------
        | Prepare space of a row according to the result of level structure
        |----------------------------------------------------------------------
        | on entry:
        |==========
        |   ( lu )  =  Pointer to a ILUSpar struct.
        |     nrow  =  the current row to deal with
        |
        | On return:
        |===========
        |
        |    lu->L->ma[nrow][...]
        |      ->U->ma[nrow][...]
        |
        | integer value returned:
        |             0   --> successful return.
        |            -1   --> memory allocation error.
        |--------------------------------------------------------------------*/
        int nzcount = lu->L->nzcount[nrow];
        lu->L->ma[nrow] = (double *)Malloc( sizeof(double)*nzcount, "mallocRow" );
        nzcount = lu->U->nzcount[nrow];
        lu->U->ma[nrow] = (double *)Malloc( sizeof(double)*nzcount, "mallocRow" );
        return 0;
    }
    /*---------------------------------------------------------------------
    |     end of mallocRow
    |--------------------------------------------------------------------*/

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

}; // Class  ILUkPreconditioner

///@}


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType>
inline std::istream& operator >> (std::istream& IStream,  ILUkPreconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
{
    return IStream;
}


/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType>
inline std::ostream& operator << (std::ostream& OStream, const  ILUkPreconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
{
    rThis.PrintInfo(OStream);
    OStream << std::endl;
    rThis.PrintData(OStream);
    return OStream;
}
///@}


}  // namespace Kratos.

#undef min
#undef max

#endif // KRATOS_MULTITHREADED_SOLVERS_APPLICATION_ILUK_PRECONDITIONER_H_INCLUDED defined

