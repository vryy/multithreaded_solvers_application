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

namespace Kratos {
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class UmfPackSolver: public DirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
{
public:
    /**
     * Counted pointer of SuperLUSolver
     */
    KRATOS_CLASS_POINTER_DEFINITION(UmfPackSolver);

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    /**
     * @param niter number of iterative refinements allowed
     */
    UmfPackSolver() {}

    /**
     * Destructor
     */
    virtual ~UmfPackSolver() {}

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

        int n = rA.size1();
        assert(n == rA.size2());
        assert(n == rB.size());
        assert(n == rX.size());

        /* nonzeros in rA */
        double* a = rA.value_data().begin();

        /* manual index vector generation */
        int* index1_vector = new int[rA.index1_data().size()];
        int* index2_vector = new int[rA.index2_data().size()];
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
        umfpack_di_symbolic(n, n, index1_vector, index2_vector, a, &Symbolic, NULL, NULL);
        umfpack_di_numeric(index1_vector, index2_vector, a, Symbolic, &Numeric, NULL, NULL);
        umfpack_di_free_symbolic(&Symbolic);
        umfpack_di_solve(UMFPACK_At, index1_vector, index2_vector, a, &rX[0], &rB[0], Numeric, NULL, NULL); //since UMFPACK uses compressed sparse column by default, so the system of A'x = b should be solved
        umfpack_di_free_numeric(&Numeric);
        
        delete [] index1_vector;
        delete [] index2_vector;

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
        KRATOS_THROW_ERROR(std::logic_error, "ERROR: This solver can be used for single RHS only", "");
        return false;
    }

    /// Return information about this object.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "UmfPack solver";
        return buffer.str();
    }

    /**
     * Print information about this object.
     */
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "UmfPack solver finished.";
    }

    /**
     * Print object's data.
     */
    void PrintData(std::ostream& rOStream) const
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

};
// Class UmfPackSolver

/**
 * input stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::istream& operator >>(
        std::istream& rIStream,
        UmfPackSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis) {
    return rIStream;
}

/**
 * output stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator <<(
        std::ostream& rOStream,
        const UmfPackSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_MULTITHREADED_SOLVERS_APPLICATION_UMFPACK_SOLVER_H_INCLUDED  defined
