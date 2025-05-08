/*
 * see multithreaded_solvers_application/LICENSE.txt                      *
 *========================================================================*
 * Created at Institute for Structural Mechanics                          *
 * Ruhr-University Bochum, Germany                                        *
 * Last modified by:    $Author: hbui $                                   *
 * Date:                $Date: Apr 19,2012 $                              *
 * Revision:            $Revision: 1.0 $                                  *
 *========================================================================*
 */

#if !defined(KRATOS_MULTITHREADED_SOLVERS_APPLICATION_UMFPACK_SOLVER_H_INCLUDED )
#define  KRATOS_MULTITHREADED_SOLVERS_APPLICATION_UMFPACK_SOLVER_H_INCLUDED

// System includes
#include <cstdio>
#include <cstdlib>
#include <cmath>

// External includes
#include <boost/numeric/bindings/traits/sparse_traits.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_sparse.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <umfpack.h>

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "utilities/openmp_utils.h"

namespace Kratos
{

template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class UmfPackSolver: public DirectSolver<TSparseSpaceType, TDenseSpaceType, ModelPart, TReordererType>
{
public:
    /**
     * Counted pointer of SuperLUSolver
     */
    KRATOS_CLASS_POINTER_DEFINITION(UmfPackSolver);

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, ModelPart, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef typename BaseType::DataType DataType;

    /**
     * @param niter number of iterative refinements allowed
     */
    UmfPackSolver() {}

    /**
     * Destructor
     */
    ~UmfPackSolver() override {}

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
        double start_solver = OpenMPUtils::GetCurrentTime();

        int n = rA.size1();
        assert(n == rA.size2());
        assert(n == rB.size());
        assert(n == rX.size());

        /* nonzeros in rA */
        DataType* a = rA.value_data().begin();

        /* manual index vector generation */
        std::vector<int> index1_vector(rA.index1_data().size());
        std::vector<int> index2_vector(rA.index2_data().size());
        std::cout << "Size of the problem: " << n << std::endl;
        std::cout << "Size of index1_vector: " << rA.index1_data().size() << std::endl;
        std::cout << "Size of index2_vector: " << rA.index2_data().size() << std::endl;
        // index1_vector is rowptr
        for (unsigned int i = 0; i < rA.index1_data().size(); ++i)
            index1_vector[i] = (int) (rA.index1_data()[i]);
        // index2_vector is colind
        for (unsigned int i = 0; i < rA.index2_data().size(); ++i)
            index2_vector[i] = (int) (rA.index2_data()[i]);

        void *Symbolic, *Numeric;
        int status;
        status = umfpack_di_symbolic(n, n, index1_vector.data(), index2_vector.data(), a, &Symbolic, NULL, NULL);
        if (status != UMFPACK_OK)
            KRATOS_ERROR << "Error with symbolic factorization, error code = " << status;

        std::vector<double> Control(UMFPACK_CONTROL);
        std::vector<double> Info(UMFPACK_INFO);
        umfpack_di_defaults( Control.data() );
        status = umfpack_di_numeric(index1_vector.data(), index2_vector.data(), a, Symbolic, &Numeric, Control.data(), Info.data());
        if (status != UMFPACK_OK)
            KRATOS_ERROR << "Error with numeric factorization, error code = " << status;

        //since UMFPACK uses compressed sparse column by default, so the system of A'x = b will be solved by default. We use UMFPACK_At to solve the transpose of A', which is A itself.
        status = umfpack_di_solve(UMFPACK_At, index1_vector.data(), index2_vector.data(), a, &rX[0], &rB[0], Numeric, NULL, NULL);
        if (status != UMFPACK_OK)
            KRATOS_ERROR << "Error with back solve, error code = " << status;

        // clean up
        if (Symbolic) umfpack_di_free_symbolic(&Symbolic);
        if (Numeric) umfpack_di_free_numeric(&Numeric);

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
    bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB) override
    {
        KRATOS_ERROR << "This solver can be used for single RHS only";
        return false;
    }

    /// Return information about this object.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "UmfPack solver";
        return buffer.str();
    }

    /**
     * Print information about this object.
     */
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "UmfPack solver finished.";
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
    UmfPackSolver& operator=(const UmfPackSolver& Other);

    /**
     * Copy constructor.
     */

}; // Class UmfPackSolver

} // namespace Kratos.

#endif // KRATOS_MULTITHREADED_SOLVERS_APPLICATION_UMFPACK_SOLVER_H_INCLUDED  defined
