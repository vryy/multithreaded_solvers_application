/*
* =======================================================================*
* kkkk   kkkk  kkkkkkkkkk   kkkkk    kkkkkkkkkk kkkkkkkkkk kkkkkkkkkK    *
* kkkk  kkkk   kkkk   kkkk  kkkkkk   kkkkkkkkkk kkkkkkkkkk kkkkkkkkkK    *
* kkkkkkkkk    kkkk   kkkk  kkkkkkk     kkkk    kkk    kkk  kkkk         *
* kkkkkkkkk    kkkkkkkkkkk  kkkk kkk	kkkk    kkk    kkk    kkkk       *
* kkkk  kkkk   kkkk  kkkk   kkkk kkkk   kkkk    kkk    kkk      kkkk     *
* kkkk   kkkk  kkkk   kkkk  kkkk  kkkk  kkkk    kkkkkkkkkk  kkkkkkkkkk   *
* kkkk    kkkk kkkk    kkkk kkkk   kkkk kkkk    kkkkkkkkkk  kkkkkkkkkk 	 *
*                                                                        *
* krATos: a fREe opEN sOURce CoDE for mULti-pHysIC aDaptIVe SoLVErS,     *
* aN extEnsIBLe OBjeCt oRiEnTEd SOlutION fOR fInITe ELemEnt fORmULatIONs *
* Copyleft by 2003 ciMNe                                                 *
* Copyleft by 2003 originary authors Copyleft by 2003 your name          *
* This library is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License as         *
* published by the Free Software Foundation; either version 2.1 of       *
* the License, or any later version.                                     *
*                                                                        *
* This library is distributed in the hope that it will be useful, but    *
* WITHOUT ANY WARRANTY; without even the implied warranty of             *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                   *
* See the GNU Lesser General Public License for more details.            *
*                                                                        *
* You should have received a copy of the GNU Lesser General Public       *
* License along with this library; if not, write to International Centre *
* for Numerical Methods in Engineering (CIMNE),                          *
* Edifici C1 - Campus Nord UPC, Gran Capità s/n, 08034 Barcelona.        *
*                                                                        *
* You can also contact us to the following email address:                *
* kratos@cimne.upc.es                                                    *
* or fax number: +34 93 401 65 17                                        *
*                                                                        *
* Created at Institute for Structural Mechanics                          *
* Ruhr-University Bochum, Germany                                        *
* Last modified by:    $Author: hbui $  				                 *
* Date:                $Date: 21 Aug 2014 $			                     *
* Revision:            $Revision: 1.4 $ 				                 *
*========================================================================*
* International Center of Numerical Methods in Engineering - CIMNE	     *
* Barcelona - Spain 							                         *
*========================================================================*
*/

#if !defined(KRATOS_MULTITHREADED_SOLVERS_APPLICATION_PARDISO_SOLVER_H_INCLUDED )
#define  KRATOS_MULTITHREADED_SOLVERS_APPLICATION_PARDISO_SOLVER_H_INCLUDED

// #define BOOST_NUMERIC_BINDINGS_SUPERLU_PRINT

// External includes

#include <stdio.h>
#include <cstdlib>
#include <cmath>

/* PARDISO prototype. */
extern "C" void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
extern "C" void pardiso     (void   *, int    *,   int *, int *,    int *, int *,
                             double *, int    *,    int *, int *,   int *, int *,
                             int *, double *, double *, int *, double *);
extern "C" void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
extern "C" void pardiso_chkvec     (int *, int *, double *, int *);
extern "C" void pardiso_printstats (int *, int *, double *, int *, int *, int *, double *, int *);

#include <boost/timer.hpp>

#include "utilities/openmp_utils.h"

#include "boost/smart_ptr.hpp"
#include "includes/ublas_interface.h"

#include <boost/numeric/bindings/traits/sparse_traits.hpp>
#include <boost/numeric/bindings/traits/matrix_traits.hpp>
#include <boost/numeric/bindings/traits/vector_traits.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_sparse.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>

// Project includes
#include "includes/define.h"
#include "linear_solvers/direct_solver.h"

#define ENABLE_PROFILING

namespace ublas = boost::numeric::ublas;

namespace Kratos
{
template< class TSparseSpaceType, class TDenseSpaceType,
          class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class PardisoSolver : public DirectSolver< TSparseSpaceType, TDenseSpaceType, TReordererType>
{
public:
    /**
     * Counted pointer of SuperLUSolver
     */
    typedef boost::shared_ptr<PardisoSolver> Pointer;

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    /**
     * @param niter number of iterative refinements allowed
     */
    PardisoSolver(unsigned int niter)
    {
        mRefinements        = niter;
        mReusePerm          = false;
        mIsInitialized      = false;
    }

    PardisoSolver()
    {
        mRefinements        = 0;
        mReusePerm          = false;
        mIsInitialized      = false;
    }

    /**
     * Destructor
     */
    virtual ~PardisoSolver() {}

    void SetReusePermutation(bool value)
    {
        mReusePerm = value;
    }

    void ResetReusePermutation()
    {
        mReusePerm = true;
        mPermutationReady = false;
    }

    /**
     * Normal solve method.
     * Solves the linear system Ax=b and puts the result on SystemVector& rX.
     * rX is also th initial guess for iterative methods.
     * @param rA. System matrix
     * @param rX. Solution vector.
     * @param rB. Right hand side vector.
     */
    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
        double start_solver = OpenMPUtils::GetCurrentTime();
        typedef boost::numeric::bindings::traits::sparse_matrix_traits<SparseMatrixType> matraits;
        typedef boost::numeric::bindings::traits::vector_traits<VectorType> mbtraits;

        /* Size checking */
        int n = matraits::size1 (rA);
        assert (n == matraits::size2 (rA));
        assert (n == mbtraits::size1 (rB));
        assert (n == mbtraits::size1 (rX));

//        /* nonzeros in rA */
//        double* a = matraits::value_storage(rA);
//        /* RHS and solution vectors */
//        double *b = mbtraits::storage(rB);
//        double *x = mbtraits::storage(rX);

        /* nonzeros in rA */
        double *a = rA.value_data().begin();
        /* RHS and solution vectors */
        double *b = &rB[0];
        double *x = &rX[0];

        /* manual index vector generation */
        int *index1_vector = new (std::nothrow) int[rA.index1_data().size()];
        int *index2_vector = new (std::nothrow) int[rA.index2_data().size()];
        std::cout << "Size of the problem: " << n << std::endl;
        std::cout << "Size of index1_vector: " << rA.index1_data().size() << std::endl;
        std::cout << "Size of index2_vector: " << rA.index2_data().size() << std::endl;
        for(unsigned int i = 0; i < rA.index1_data().size(); ++i )
            index1_vector[i] = (int)(rA.index1_data()[i]) + 1;
        for(unsigned int i = 0; i < rA.index2_data().size(); ++i )
            index2_vector[i] = (int)(rA.index2_data()[i]) + 1;

        /**
         *  Matrix type flag:
         * 1    real and structurally symmetric
         * 2    real and symmetic positive definite
         * -2   real and symmetric indefinite
         * 3    complex and structurally symmetric
         * 4    complex and Hermitian positive definite
         * -4   complex and Hermitian indefinite
         * 6    complex and symmetic
         * 11   real and nonsymmetric
         * 13   complex and nonsymmetric
         */
        int mtype = 1;

        /* Number of right hand sides */
        int nrhs = 1;

        /* Internal solver memory pointer pt, */
        /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
        /* or void *pt[64] should be OK on both architectures */
        void *pt[64];

        /* Pardiso control parameters */
        int iparm[64];
        double dparm[64];
        int maxfct, mnum, phase, error = 0, msglvl;

        /* Auxiliary variables */
        int i;
        double ddum;    // Double dummy
        int* perm;  // Permutation array

        /* -------------------------------------------------------------------- */
        /* ..  Setup Pardiso control parameters and initialize the solvers      */
        /*     internal adress pointers. This is only necessary for the FIRST   */
        /*     call of the PARDISO solver.                                      */
        /* ---------------------------------------------------------------------*/
        int solver = 0; /* use sparse direct solver */
        pardisoinit (pt, &mtype, &solver, iparm, dparm, &error);
        if (error != 0)
        {
            if (error == -10 )
                printf("No license file found \n");
            if (error == -11 )
                printf("License is expired \n");
            if (error == -12 )
                printf("Wrong username or hostname \n");
            return 1;
        }
        else
            printf("[PARDISO]: License check was successful ... \n");

        /* -------------------------------------------------------------------- */
        /*  .. pardiso_chk_matrix(...)                                          */
        /*     Checks the consistency of the given matrix.                      */
        /*     Use this functionality only for debugging purposes               */
        /* -------------------------------------------------------------------- */
        pardiso_chkmatrix (&mtype, &n, a, index1_vector, index2_vector, &error);
        if (error != 0)
            KRATOS_THROW_ERROR(std::logic_error, "ERROR in consistency of matrix:", error);

        /* -------------------------------------------------------------------- */
        /* ..  pardiso_chkvec(...)                                              */
        /*     Checks the given vectors for infinite and NaN values             */
        /*     Input parameters (see PARDISO user manual for a description):    */
        /*     Use this functionality only for debugging purposes               */
        /* -------------------------------------------------------------------- */
        pardiso_chkvec (&n, &nrhs, b, &error);
        if (error != 0)
            KRATOS_THROW_ERROR(std::logic_error, "ERROR in right hand side:", error);

        /* -------------------------------------------------------------------- */
        /* .. pardiso_printstats(...)                                           */
        /*    prints information on the matrix to STDOUT.                       */
        /*    Use this functionality only for debugging purposes                */
        /* -------------------------------------------------------------------- */
        pardiso_printstats (&mtype, &n, a, index1_vector, index2_vector, &nrhs, b, &error);
        if (error != 0)
            KRATOS_THROW_ERROR(std::logic_error, "ERROR right hand side:", error);

        /* -------------------------------------------------------------------- */
        /* Setup Pardiso control parameters.                                    */
        /* -------------------------------------------------------------------- */
//        for (i = 0; i < 64; i++) iparm[i] = 0;
        iparm[0] = 1; /* No solver default, user supply entries of iparm */
        //iparm[1] = 0; /* Fill-in reordering by minimum degree algorithm*/
        iparm[1] = 2; /* Fill-in reordering from METIS */
//        iparm[1] = 3; /* Fill-in reordering from METIS with OpenMP support*/ //not available in original pardiso

        /* Numbers of processors, value of OMP_NUM_THREADS */
        iparm[2] = OpenMPUtils::GetNumThreads(); //omp_get_max_threads();
        std::cout << "Number of threads/procs (for MKL): " << iparm[2] << std::endl;

        if( mRefinements > 0 )
            iparm[3] = 1; /* iterative-direct algorithm */
        else
            iparm[3] = 0; /* no iterative-direct algorithm */

        if(mReusePerm == false)
        {
            iparm[4] = 0; /* No user fill-in reducing permutation */
            perm = new int[1];
        }
        else
        {
            perm = new int[n];
            if(mPermutationReady == false)
            {
                iparm[4] = 2; /* return permutation array to perm */
            }
            else
            {
                iparm[4] = 1; /* use the ready permutation array */
                std::copy(mPerm.begin(), mPerm.end(), perm);
            }
        }

        iparm[5] = 0;               /* Write solution into x */
        iparm[6] = 0;               /* Not in use */
        iparm[7] = mRefinements;    /* Max numbers of iterative refinement steps */
        iparm[8] = 0;               /* Not in use */
        iparm[9] = 13;              /* Perturb the pivot elements with 1E-13 */
        iparm[10] = 1;              /* Use nonsymmetric permutation and scaling MPS */
        iparm[11] = 0;              /* Solve Ax = b; set to 1 if solve A^T x = b' */
        iparm[12] = 1;              /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
        iparm[13] = 0;              /* Output: Number of perturbed pivots */
        iparm[14] = 0;              /* Not in use */
        iparm[15] = 0;              /* Not in use */
        iparm[16] = 0;              /* Not in use */
        iparm[17] = -1;             /* Output: Number of nonzeros in the factor LU */
        iparm[18] = -1;             /* Output: Mflops for LU factorization */
        iparm[19] = 0;              /* Output: Numbers of CG Iterations */
        iparm[23] = 1;              /* Parallel numerical factorization */
        iparm[24] = 1;              /* Parallel forward/backward solve */
        iparm[27] = 1;              /* Parallel reordering for Metis */
        iparm[31] = 0;              /* Use sparse direct solver; set to 1 for iterative solver */
        maxfct = 1;                 /* Maximum number of numerical factorizations. */
        mnum = 1;                   /* Which factorization to use. */
        msglvl = 0;                 /* Print statistical information to the screen */
        error = 0;                  /* Initialize error flag */

        /* -------------------------------------------------------------------- */
        /* .. Reordering and Symbolic Factorization. This step also allocates   */
        /* all memory that is necessary for the factorization.                  */
        /* -------------------------------------------------------------------- */
        #ifdef ENABLE_PROFILING
        double tmp_timing_solver = OpenMPUtils::GetCurrentTime();
        #endif
        phase = 11;
        pardiso (pt, &maxfct, &mnum, &mtype, &phase,
                 &n, a, index1_vector, index2_vector, perm, &nrhs,
                 iparm, &msglvl, &ddum, &ddum, &error, dparm);
        if (error != 0)
            KRATOS_THROW_ERROR(std::logic_error, "ERROR during symbolic factorization:", error);

        if(mReusePerm == true)
        {
            if(mPermutationReady == false)
            {
                mPerm.resize(n);
                std::copy(perm, perm + n, mPerm.begin());
                mPermutationReady = true;
            }
        }

        #ifdef ENABLE_PROFILING
        printf("Reordering completed ... %f s\n", OpenMPUtils::GetCurrentTime() - tmp_timing_solver);
        #else
        printf("Reordering completed ...\n");
        #endif
        printf("Peak memory symbolic factorization          = %d KBs\n", iparm[14]);
        printf("Permanent memory symbolic factorization     = %d KBs\n", iparm[15]);
        printf("Number of nonzeros in factors               = %d\n", iparm[17]);
        printf("Number of factorization MFLOPS              = %d\n", iparm[18]);

        /* -------------------------------------------------------------------- */
        /* .. Numerical factorization.                                          */
        /* -------------------------------------------------------------------- */
        #ifdef ENABLE_PROFILING
        tmp_timing_solver = OpenMPUtils::GetCurrentTime();
        #endif
        phase = 22;
        pardiso (pt, &maxfct, &mnum, &mtype, &phase,
                 &n, a, index1_vector, index2_vector, perm, &nrhs,
                 iparm, &msglvl, &ddum, &ddum, &error, dparm);
        if (error != 0)
            KRATOS_THROW_ERROR(std::logic_error, "ERROR during numerical factorization:", error);

        #ifdef ENABLE_PROFILING
        printf("Factorization completed ... %f s\n", OpenMPUtils::GetCurrentTime() - tmp_timing_solver);
        #else
        printf("Factorization completed ...\n");
        #endif
        printf("Number of perturbed pivots                  = %d\n", iparm[13]);
        printf("Memory numerical factorization and solution = %d KBs\n", iparm[16]);

        /* -------------------------------------------------------------------- */
        /* .. Back substitution and iterative refinement.                       */
        /* -------------------------------------------------------------------- */
        #ifdef ENABLE_PROFILING
        tmp_timing_solver = OpenMPUtils::GetCurrentTime();
        #endif
        phase = 33;
        iparm[7] = 2; /* Max numbers of iterative refinement steps. */
        pardiso (pt, &maxfct, &mnum, &mtype, &phase,
                 &n, a, index1_vector, index2_vector, perm, &nrhs,
                 iparm, &msglvl, b, x, &error, dparm);
        if (error != 0)
            KRATOS_THROW_ERROR(std::logic_error, "ERROR during solution:", error);

        #ifdef ENABLE_PROFILING
        printf("Solve completed ... %f s\n", OpenMPUtils::GetCurrentTime() - tmp_timing_solver);
        #else
        printdf("Solve completed ...\n");
        #endif
        printf("Number of performed iterative refinement steps = %d\n", iparm[6]);
        if(iparm[19] < 0)
            printf("Warning: Iterations executed, but CG/CGS failed, iparm[19] = %d\n", iparm[19]);
        else
            printf("Number of CGS iteration                        = %d\n", iparm[19]);
        printf("Number of positive eigenvalues                 = %d\n", iparm[21]); //remarks: only for symmetric matrix (mtype = -2), as of pardiso 5.0.0
        printf("Number of negative eigenvalues                 = %d\n", iparm[22]); //remarks: only for symmetric matrix (mtype = -2), as of pardiso 5.0.0

        /* -------------------------------------------------------------------- */
        /* .. Termination and release of memory. */
        /* -------------------------------------------------------------------- */
        phase = -1; /* Release internal memory. */
        pardiso (pt, &maxfct, &mnum, &mtype, &phase,
                 &n, &ddum, index1_vector, index2_vector, perm, &nrhs,
                 iparm, &msglvl, &ddum, &ddum, &error, dparm);
        delete [] index1_vector;
        delete [] index2_vector;
        delete [] perm;

        std::cout << "#### SOLVER TIME: " << OpenMPUtils::GetCurrentTime() - start_solver << " ####" << std::endl;
        return true;
    }

    /**
     * Multi solve method for solving a set of linear systems with same coefficient matrix.
     * Solves the linear system Ax=b and puts the result on SystemVector& rX.
     * rX is also th initial guess for iterative methods.
     * @param rA. System matrix
     * @param rX. Solution vector.
     * @param rB. Right hand side vector.
     */
    bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
    {
        //TODO
        return false;
    }

    /// Return information about this object.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "PARDISO solver";
        return buffer.str();
    }

    /**
     * Print information about this object.
     */
    void  PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "PARDISO solver finished.";
    }

    /**
     * Print object's data.
     */
    void  PrintData(std::ostream& rOStream) const
    {
    }

private:

    int mRefinements;
    bool mReusePerm; //reuse the permutation vector for subsequent steps
    bool mPermutationReady; //indicate that the permutation vector is ready or not
    boost::numeric::ublas::vector<int> mPerm; //array to store the permutation arrary in case the permutation is reused
    bool mIsInitialized;

    /**
     * Assignment operator.
     */
    PardisoSolver& operator=(const PardisoSolver& Other);

    /**
     * Copy constructor.
     */
//             PardisoSolver(const ParallelSuperLUSolver& Other);

}; // Class PardisoSolver


///**
// * input stream function
// */
//template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
//inline std::istream& operator >> (std::istream& rIStream, PardisoSolver< TSparseSpaceType,
//                                  TDenseSpaceType, TReordererType>& rThis)
//{
//    return rIStream;
//}

///**
// * output stream function
// */
//template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
//inline std::ostream& operator << (std::ostream& rOStream,
//                                  const PardisoSolver<TSparseSpaceType,
//                                  TDenseSpaceType, TReordererType>& rThis)
//{
//    rThis.PrintInfo(rOStream);
//    rOStream << std::endl;
//    rThis.PrintData(rOStream);

//    return rOStream;
//}


}  // namespace Kratos.

#endif // KRATOS_MULTITHREADED_SOLVERS_APPLICATION_PARDISO_SOLVER_H_INCLUDED  defined


