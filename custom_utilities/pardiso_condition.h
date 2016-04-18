#if !defined(KRATOS_MULTITHREADED_SOLVERS_APPLICATION_PARDISO_CONDITION_H_INCLUDED )
#define  KRATOS_MULTITHREADED_SOLVERS_APPLICATION_PARDISO_CONDITION_H_INCLUDED

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

#include "utilities/openmp_utils.h"

#define ENABLE_PROFILING

namespace ublas = boost::numeric::ublas;

namespace Kratos
{
template< class TSparseMatrixType, class TVectorType>
class PardisoCondition
{
public:
    /**
     * Counted pointer of SuperLUSolver
     */
    typedef boost::shared_ptr<PardisoCondition> Pointer;

    /**
     * @param niter number of iterative refinements allowed
     */
    PardisoCondition(unsigned int niter)
    {
        mRefinements = niter;
    }

    PardisoCondition()
    {
        mRefinements = 0;
    }

    /**
     * Destructor
     */
    virtual ~PardisoCondition() {}

    /**
     * Normal solve method.
     * Solves the linear system Ax=b and puts the result on SystemVector& rX.
     * rX is also th initial guess for iterative methods.
     * @param rA. System matrix
     * @param rX. Solution vector.
     * @param rB. Right hand side vector.
     */
    double hager_norm1_inv_A(TSparseMatrixType& rA)
    {
        double start_solver = OpenMPUtils::GetCurrentTime();

        /* Size checking */
        unsigned int i;
        int n = rA.size1();
        assert (n == rA.size2());

        /* nonzeros in rA */
        double* a = rA.value_data().begin();

        /* manual index vector generation */
        int* index1_vector = new (std::nothrow) int[rA.index1_data().size()];
        int* index2_vector = new (std::nothrow) int[rA.index2_data().size()];
        std::cout << "Size of the problem: " << n << std::endl;
        std::cout << "Size of index1_vector: " << rA.index1_data().size() << std::endl;
        std::cout << "Size of index2_vector: " << rA.index2_data().size() << std::endl;
        for(i = 0; i < rA.index1_data().size(); ++i )
            index1_vector[i] = (int)(rA.index1_data()[i]) + 1;
        for(i = 0; i < rA.index2_data().size(); ++i )
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
        int mtype = 11;
        
        /* Number of right hand sides */
        int nrhs = 1;
        
        /* Internal solver memory pointer pt, */
        /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
        /* or void *pt[64] should be OK on both architectures */
        void *pt[64];
        
        /* Pardiso control parameters */
        int iparm[64];
        double dparm[64];
        int maxfct, mnum, phase, error, msglvl;
        
        /* Auxiliary variables */
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
            
        iparm[4] = 0; /* No user fill-in reducing permutation */
        perm = new int[1];
        iparm[5] = 0;               /* Write solution into x */
        iparm[6] = 0;               /* Not in use */
        iparm[7] = mRefinements;    /* Max numbers of iterative refinement steps */
        iparm[8] = 0;               /* Not in use */
        iparm[9] = 13;              /* Perturb the pivot elements with 1E-13 */
        iparm[10] = 0;              /* Use nonsymmetric permutation and scaling MPS */
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
        /* .. W. Hager algorithm.                                               */
        /* -------------------------------------------------------------------- */
        #ifdef ENABLE_PROFILING
        tmp_timing_solver = OpenMPUtils::GetCurrentTime();
        #endif
        phase = 33;
        iparm[7] = 2; /* Max numbers of iterative refinement steps. */
        
        double* b = new double[n];
        double* x = new double[n];
        
        for(i = 0; i < n; ++i)
            b[i] = 1.0 / (double)(n);
        
        int i1 = -1, i2, cnt = 0;
        double c1 = 0.0, c2;
        while(true)
        {
            /* Solve Ax = b */
            iparm[11] = 0;
            pardiso (pt, &maxfct, &mnum, &mtype, &phase,
                 &n, a, index1_vector, index2_vector, perm, &nrhs,
                 iparm, &msglvl, b, x, &error, dparm);
            if (error != 0)
                KRATOS_THROW_ERROR(std::logic_error, "ERROR during solution:", error);
            
            /* Compute norm_l1 of x */
            c2 = 0.0;
            for(i = 0; i < n; ++i)
                c2 += fabs(x[i]);
            
            /* Compute rhs vector of sign */
            for(i = 0; i < n; ++i)
                if(x[i] < 0.0)
                    b[i] = -1.0;
                else
                    b[i] = 1.0;
            
            /* Solve A^Tx = b */
            iparm[11] = 1;
            pardiso (pt, &maxfct, &mnum, &mtype, &phase,
                 &n, a, index1_vector, index2_vector, perm, &nrhs,
                 iparm, &msglvl, b, x, &error, dparm);
            if (error != 0)
                KRATOS_THROW_ERROR(std::logic_error, "ERROR during solution:", error);

            /* Find the index of maximum entry in x */
            i2 = 0;
            for(i = 1; i < n; ++i)
                if(fabs(x[i2] ) < fabs(x[i]))
                    i2 = i;

            if(0 <= i1)
                if ( i1 == i2 || c2 <= c1 )
                    break;
            i1 = i2;
            c1 = c2;

            for(i = 0; i < n; ++i)
                b[i] = 0.0;
            b[i1] = 1.0;
            ++cnt;
        }
        std::cout << std::endl;
        std::cout << "Condition iteration converged after " << cnt << " iteration(s)" << std::endl;
//        KRATOS_WATCH(norm_1(rA))
        #ifdef ENABLE_PROFILING
        printf("Condition completed ... %f s\n", OpenMPUtils::GetCurrentTime() - tmp_timing_solver);
        #else
        printdf("Condition completed ...\n");
        #endif
        
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
        delete [] b;
        delete [] x;
        
        std::cout << "#### CONDITION SOLVER TIME: " << OpenMPUtils::GetCurrentTime() - start_solver << " ####" << std::endl;
        return c2;
    }

    /// Return information about this object.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "PARDISO condition number utility";
        return buffer.str();
    }
    
    /**
     * Print information about this object.
     */
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "PARDISO condition finished.";
    }

    /**
     * Print object's data.
     */
    void PrintData(std::ostream& rOStream) const
    {
    }

private:

    int mRefinements;
    
    /**
     * Assignment operator.
     */
    PardisoCondition& operator=(const PardisoCondition& Other);

    /**
     * Copy constructor.
     */
//    PardisoCondition(const ParallelSuperLUSolver& Other);

}; // Class PardisoCondition

///**
// * input stream function
// */
//template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
//inline std::istream& operator >> (std::istream& rIStream, PardisoCondition< TSparseSpaceType,
//                                  TDenseSpaceType, TReordererType>& rThis)
//{
//    return rIStream;
//}

///**
// * output stream function
// */
//template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
//inline std::ostream& operator << (std::ostream& rOStream,
//                                  const PardisoCondition<TSparseSpaceType,
//                                  TDenseSpaceType, TReordererType>& rThis)
//{
//    rThis.PrintInfo(rOStream);
//    rOStream << std::endl;
//    rThis.PrintData(rOStream);

//    return rOStream;
//}

}  // namespace Kratos.

#endif // KRATOS_MULTITHREADED_SOLVERS_APPLICATION_PARDISO_CONDITION_H_INCLUDED  defined 


