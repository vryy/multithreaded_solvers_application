/*
* =======================================================================*
* kkkk   kkkk  kkkkkkkkkk   kkkkk    kkkkkkkkkk kkkkkkkkkk kkkkkkkkkK    *
* kkkk  kkkk   kkkk   kkkk  kkkkkk   kkkkkkkkkk kkkkkkkkkk kkkkkkkkkK    *
* kkkkkkkkk    kkkk   kkkk  kkkkkkk     kkkk    kkk    kkk  kkkk         *
* kkkkkkkkk    kkkkkkkkkkk  kkkk kkk    kkkk    kkk    kkk    kkkk       *
* kkkk  kkkk   kkkk  kkkk   kkkk kkkk   kkkk    kkk    kkk      kkkk     *
* kkkk   kkkk  kkkk   kkkk  kkkk  kkkk  kkkk    kkkkkkkkkk  kkkkkkkkkk   *
* kkkk    kkkk kkkk    kkkk kkkk   kkkk kkkk    kkkkkkkkkk  kkkkkkkkkk      *
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
* Last modified by:    $Author: hbui $                                   *
* Date:                $Date: 19 Aug 2014 $                              *
* Revision:            $Revision: 1.5 $                                  *
*========================================================================*
* International Center of Numerical Methods in Engineering - CIMNE       *
* Barcelona - Spain                                                      *
*========================================================================*
*/

#if !defined(KRATOS_MULTITHREADED_SOLVERS_APPLICATION_SUPERLU_MT_SOLVER_H_INCLUDED )
#define  KRATOS_MULTITHREADED_SOLVERS_APPLICATION_SUPERLU_MT_SOLVER_H_INCLUDED

// External includes

#include "includes/ublas_interface.h"
#include "utilities/openmp_utils.h"

#include "SRC/pdsp_defs.h"

// Project includes
#include "includes/define.h"
#include "linear_solvers/direct_solver.h"

namespace ublas = boost::numeric::ublas;

namespace Kratos
{
template< class TSparseSpaceType, class TDenseSpaceType,
          class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class SuperLUMTSolver : public DirectSolver< TSparseSpaceType,
    TDenseSpaceType, TReordererType>
{
public:
    /**
     * Counted pointer of SuperLUMTSolver
     */
    KRATOS_CLASS_POINTER_DEFINITION(SuperLUMTSolver);

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    /**
     * Default constructor
     */
    SuperLUMTSolver() {}

    /**
     * Destructor
     */
    ~SuperLUMTSolver() override {}

    /**
     * Normal solve method.
     * Solves the linear system Ax=b and puts the result on SystemVector& rX.
     * rX is also th initial guess for iterative methods.
     * @param rA. System matrix
     * @param rX. Solution vector.
     * @param rB. Right hand side vector.
     */
    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        std::cout << "matrix size in solver:  " << rA.size1() << std::endl;
        std::cout << "RHS size in solver SLU_MT: " << rB.size() << std::endl;

        if(this->IsNotConsistent(rA, rX, rB))
            return false;

        SuperMatrix Aslu, L, U;
        SuperMatrix B, X;
        SCPformat   *Lstore;
        NCPformat   *Ustore;
        int         nprocs;
        fact_t      fact;
        trans_t     trans;
        yes_no_t    refact, usepr;
        equed_t     equed;
        int         *perm_c; /* column permutation vector */
        int         *perm_r; /* row permutations from partial pivoting */
        void        *work;
        superlumt_options_t superlumt_options;
        int         info, lwork, nrhs, panel_size, relax;
        int         permc_spec;
        int         firstfact;
        double      *R, *C;
        double      *ferr, *berr;
        double      u, drop_tol, rpg, rcond;
        superlu_memusage_t superlu_memusage;
        double      start, stop;

        /* Default parameters to control factorization. */
        nprocs      = OpenMPUtils::GetNumThreads();
        fact        = EQUILIBRATE;
        trans       = TRANS;
        equed       = NOEQUIL;
        refact      = NO;
        panel_size  = sp_ienv(1);
        relax       = sp_ienv(2);
        u           = 1.0;
        usepr       = NO;
        drop_tol    = 0.0;
        lwork       = 0;
        nrhs        = 1;

        printf("SuperLU_MT is called with %d threads\n", nprocs);
        if ( lwork > 0 )
        {
            work = SUPERLU_MALLOC(lwork);
            printf("Use work space of size LWORK = %d bytes\n", lwork);
            if ( !work )
                SUPERLU_ABORT("DLINSOLX: cannot allocate work[]");
        }

        //create a copy of the matrix
        start = OpenMPUtils::GetCurrentTime();
        std::vector<int> index1_vector(rA.index1_data().size());
        std::vector<int> index2_vector(rA.index2_data().size());
//         double *values_vector = new (std::nothrow) double[rA.value_data().size()];

        for( unsigned int i = 0; i < rA.index1_data().size(); i++ )
            index1_vector[i] = (int)rA.index1_data()[i];

        for( unsigned int i = 0; i < rA.index2_data().size(); i++ )
            index2_vector[i] = (int)rA.index2_data()[i];

        firstfact = (fact == FACTORED || refact == YES);

        dCreate_CompCol_Matrix (&Aslu, rA.size1(), rA.size2(),
                                rA.nnz(),
                                rA.value_data().begin(),
                                index2_vector.data(), //can not avoid a copy as ublas uses unsigned int internally
                                index1_vector.data(), //can not avoid a copy as ublas uses unsigned int internally
                                SLU_NC, SLU_D, SLU_GE
                               );

        dCreate_Dense_Matrix (&B, rB.size(), nrhs, &rB[0], rB.size(), SLU_DN, SLU_D, SLU_GE);
        dCreate_Dense_Matrix (&X, rX.size(), nrhs, &rX[0], rX.size(), SLU_DN, SLU_D, SLU_GE);

        //allocate memory for permutation arrays
        if (!(perm_r = intMalloc(rA.size1()))) SUPERLU_ABORT("Malloc fails for perm_r[].");
        if (!(perm_c = intMalloc(rA.size2()))) SUPERLU_ABORT("Malloc fails for perm_c[].");
        if (!(R = (double *) SUPERLU_MALLOC(Aslu.nrow * sizeof(double))))
            SUPERLU_ABORT("SUPERLU_MALLOC fails for R[].");
        if (!(C = (double *) SUPERLU_MALLOC(Aslu.ncol * sizeof(double))))
            SUPERLU_ABORT("SUPERLU_MALLOC fails for C[].");
        if (!(ferr = (double *) SUPERLU_MALLOC(nrhs * sizeof(double))))
            SUPERLU_ABORT("SUPERLU_MALLOC fails for ferr[].");
        if (!(berr = (double *) SUPERLU_MALLOC(nrhs * sizeof(double))))
            SUPERLU_ABORT("SUPERLU_MALLOC fails for berr[].");

        stop = OpenMPUtils::GetCurrentTime();
        printf("SuperLU_MT: allocation for system of equations comleted: %f s\n", stop - start);
        start = stop;

        /*
         * Get column permutation vector perm_c[], according to permc_spec:
         *   permc_spec = 0: natural ordering
         *   permc_spec = 1: minimum degree ordering on structure of A'*A
         *   permc_spec = 2: minimum degree ordering on structure of A'+A
         *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
         */
        permc_spec = 1;
        get_perm_c(permc_spec, &Aslu, perm_c);

        stop = OpenMPUtils::GetCurrentTime();
        printf("SuperLU_MT: reordering completed: %f s\n", stop - start);
        start = stop;

        superlumt_options.nprocs = nprocs;
        superlumt_options.fact = fact;
        superlumt_options.trans = trans;
        superlumt_options.refact = refact;
        superlumt_options.panel_size = panel_size;
        superlumt_options.relax = relax;
        superlumt_options.diag_pivot_thresh = u;
        superlumt_options.usepr = usepr;
        superlumt_options.drop_tol = drop_tol;
        superlumt_options.SymmetricMode = NO;
        superlumt_options.PrintStat = NO;
        superlumt_options.perm_c = perm_c;
        superlumt_options.perm_r = perm_r;
        superlumt_options.work = work;
        superlumt_options.lwork = lwork;

        /*
         * Solve the system and compute the condition number
         * and error bounds using pdgssvx.
         */
        pdgssvx(nprocs, &superlumt_options, &Aslu, perm_c, perm_r,
            &equed, R, C, &L, &U, &B, &X, &rpg, &rcond,
            ferr, berr, &superlu_memusage, &info);

        stop = OpenMPUtils::GetCurrentTime();
        printf("SuperLU_MT: solve completed: %f s\n", stop - start);
        start = stop;

        printf("pdgssvx(): info %d\n", info);

        //print matrix analysis result and deallocate memory
        if ( info == 0 || info == rX.size() + 1 )
        {
            printf("Recip. pivot growth = %e\n", rpg);
            printf("Recip. condition number = %e\n", rcond);
            printf("%8s%16s%16s\n", "rhs", "FERR", "BERR");
            for (unsigned int i = 0; i < nrhs; ++i)
                printf("%8d%16e%16e\n", i+1, ferr[i], berr[i]);

            Lstore = (SCPformat *) L.Store;
            Ustore = (NCPformat *) U.Store;
            printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
            printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
            printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz - rX.size());
            printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n",
               superlu_memusage.for_lu/1e6, superlu_memusage.total_needed/1e6,
               superlu_memusage.expansions);

            fflush(stdout);
        }
        else if ( info > 0 && lwork == -1 )
            printf("** Estimated memory: %d bytes\n", info - rX.size());

        SUPERLU_FREE (perm_r);
        SUPERLU_FREE (perm_c);
        SUPERLU_FREE (R);
        SUPERLU_FREE (C);
        SUPERLU_FREE (ferr);
        SUPERLU_FREE (berr);
//        Destroy_CompCol_Matrix(&A);
        Destroy_SuperMatrix_Store(&Aslu);
        Destroy_SuperMatrix_Store(&B);
        if ( lwork >= 0 )
        {
            Destroy_SuperNode_SCP(&L);
            Destroy_CompCol_NCP(&U);
        }

//         delete [] b_vector;

        //CHECK WITH VALGRIND IF THIS IS NEEDED ...or if it is done by the lines above
        //deallocate tempory storage used for the matrix
//                 if(b_vector!=NULL) delete [] index1_vector;
// //           if(b_vector!=NULL) delete [] index2_vector;
//           if(b_vector!=NULL) delete [] values_vector;
//         if(b_vector!=NULL) delete [] b_vector;

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
    bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB) override
    {
        //TODO
        bool is_solved = false;

        KRATOS_ERROR << "This solver can be used for single RHS only";

        return is_solved;
    }

    /// Return information about this object.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "SuperLU_MT solver";
        return buffer.str();
    }

    /**
     * Print information about this object.
     */
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SuperLU_MT solver finished.";
    }

    /**
     * Print object's data.
     */
    void PrintData(std::ostream& rOStream) const override
    {
    }

private:

    /**
     * Assignment operator.
     */
    SuperLUMTSolver& operator=(const SuperLUMTSolver& Other);

    /**
     * Copy constructor.
     */
    SuperLUMTSolver(const SuperLUMTSolver& Other);

}; // Class SkylineLUFactorizationSolver

}  // namespace Kratos.

#endif // KRATOS_MULTITHREADED_SOLVERS_APPLICATION_SUPERLU_SOLVER_H_INCLUDED  defined
