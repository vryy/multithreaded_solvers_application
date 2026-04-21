/*
 * see multithreaded_solvers_application/LICENSE.txt                      *
 *========================================================================*
 * Created at Institute for Structural Mechanics                          *
 * Ruhr-University Bochum, Germany                                        *
 * Last modified by:    $Author: hbui $                                   *
 * Date:                $Date: Feb 25,2017 $                              *
 * Revision:            $Revision: 1.0 $                                  *
 *========================================================================*
 */

#if !defined(KRATOS_MULTITHREADED_SOLVERS_APPLICATION_DIAGONAL_FIT_SOLVER_H_INCLUDED )
#define  KRATOS_MULTITHREADED_SOLVERS_APPLICATION_DIAGONAL_FIT_SOLVER_H_INCLUDED

// System includes
#include <cstdio>
#include <cstdlib>
#include <cmath>

// External includes
#include <boost/numeric/bindings/traits/sparse_traits.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_sparse.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"

#define CHECK_INCOMPATIBLE_ROWS

namespace Kratos {

/*
This solver tries to change the diagonal before solving. If a row is zero, then the diagonal is set by the average value.
If the row is zero but the corresponding RHS is nonzero (compared with a tolerance), error will be thrown.
It is essential to remove the singularity of the linear system. Especially when it's arised from methods like immersed boundary method.
Additionally, it also change the sign of the row if the negative diagonal is detected.
*/
template<class TSparseSpaceType, class TDenseSpaceType, class TModelPartType, class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class DiagonalFitSolver: public DirectSolver<TSparseSpaceType, TDenseSpaceType, TModelPartType, TReordererType>
{
public:
    /**
     * Counted pointer of SuperLUSolver
     */
    KRATOS_CLASS_POINTER_DEFINITION(DiagonalFitSolver);

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TModelPartType, TReordererType> BaseType;

    typedef typename BaseType::SparseMatrixType SparseMatrixType;

    typedef typename BaseType::VectorType VectorType;

    typedef typename BaseType::DenseMatrixType DenseMatrixType;

    typedef typename BaseType::SizeType SizeType;

    typedef typename BaseType::IndexType IndexType;

    typedef typename BaseType::DataType DataType;

    typedef typename BaseType::ValueType ValueType;

    /**
     * Default Constructor
     */
    DiagonalFitSolver(typename BaseType::Pointer pLinearSolver)
    : mpLinearSolver(pLinearSolver), mTol(1.0e-12)
    {}

    /**
     * Extended Constructor
     */
    DiagonalFitSolver(typename BaseType::Pointer pLinearSolver, const ValueType Tol)
    : mpLinearSolver(pLinearSolver), mTol(Tol)
    {}

    /**
     * Destructor
     */
    ~DiagonalFitSolver() override {}

    /** Some solvers may require a minimum degree of knowledge of the structure of the matrix. To make an example
     * when solving a mixed u-p problem, it is important to identify the row associated to v and p.
     * another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers
     * which require knowledge on the spatial position of the nodes associated to a given dof.
     * This function tells if the solver requires such data
     */
    bool AdditionalPhysicalDataIsNeeded() override
    {
        #ifdef CHECK_INCOMPATIBLE_ROWS
        return true;
        #else
        return mpLinearSolver->AdditionalPhysicalDataIsNeeded();
        #endif
    }

    /** Some solvers may require a minimum degree of knowledge of the structure of the matrix. To make an example
     * when solving a mixed u-p problem, it is important to identify the row associated to v and p.
     * another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers
     * which require knowledge on the spatial position of the nodes associated to a given dof.
     * This function is the place to eventually provide such data
     */
    void ProvideAdditionalData(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB,
        typename TModelPartType::DofsArrayType& rdof_set,
        TModelPartType& r_model_part
    ) override
    {
        mpLinearSolver->ProvideAdditionalData(rA, rX, rB, rdof_set, r_model_part);

        #ifdef CHECK_INCOMPATIBLE_ROWS
        std::set<std::size_t> zero_rows = GetZeroRows(rA, mTol);
        std::set<std::size_t> incompatible_rows;
        std::set<std::size_t> incompatible_nodes;

        for(auto it = zero_rows.begin(); it != zero_rows.end(); ++it)
            if(rB(*it) != 0.0)
                incompatible_rows.insert(*it);

        std::cout << "incompatible_rows (at tol = " << mTol << "):";
        for(auto it = incompatible_rows.begin(); it != incompatible_rows.end(); ++it)
            std::cout << " " << *it;
        std::cout << std::endl;

        if (incompatible_rows.size() > 0)
            KRATOS_ERROR << "Incompatiable rows (diagonal=0, rhs!=0) are detected in the linear system";

        for(auto it = incompatible_rows.begin(); it != incompatible_rows.end(); ++it)
        {
            for(auto it2 = rdof_set.begin(); it2 != rdof_set.end(); ++it2)
            {
                if(it2->EquationId() == *it)
                {
                    incompatible_nodes.insert(it2->Id());
                    break;
                }
            }
        }

        std::cout << "incompatible_nodes (at tol = " << mTol << "):";
        for(std::set<std::size_t>::iterator it = incompatible_nodes.begin(); it != incompatible_nodes.end(); ++it)
            std::cout << " " << *it;
        std::cout << std::endl;
        #endif
    }

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
        std::size_t n = rA.size1();
        double diagonal_ave = 0.0;
        std::size_t diagonal_count = 0;
        std::size_t num_negative_rows = 0;
        std::set<std::size_t> zero_rows = GetZeroRows(rA, mTol);

        // reverse sign of the negative rows if needed
        for(std::size_t i = 0; i < n; ++i)
        {
            if(rA(i, i) < 0.0)
            {
                noalias( row(rA, i) ) = -row(rA, i);
                rB(i) = -rB(i);
                ++num_negative_rows;
            }

            diagonal_ave += rA(i, i);
            ++diagonal_count;
        }
        diagonal_ave /= diagonal_count;

        std::cout << num_negative_rows << " negative rows detected; reverse sign completed" << std::endl;

        // adjust the zero rows
        std::cout << "zero_rows:";
        for(std::set<std::size_t>::iterator it = zero_rows.begin(); it != zero_rows.end(); ++it)
        {
            rA(*it, *it) = diagonal_ave;
            std::cout << " " << *it;
        }
        std::cout << std::endl;

        // print the zero rows
        std::cout << zero_rows.size() << " zero rows (at tol = " << mTol << ") detected; set diagonal to " << diagonal_ave
                  << " and respective rhs to zero" << std::endl;

        // solve the system
        mpLinearSolver->Solve(rA, rX, rB);

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
        KRATOS_THROW_ERROR(std::logic_error, "ERROR: This solver can be used for single RHS only", "");
        return false;
    }

    /// Return information about this object.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "DiagonalFitSolver";
        return buffer.str();
    }

    /**
     * Print information about this object.
     */
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "DiagonalFitSolver finished.";
    }

    /**
     * Print object's data.
     */
    void PrintData(std::ostream& rOStream) const override
    {
    }

private:

    typename BaseType::Pointer mpLinearSolver;
    ValueType mTol;

    template<typename TMatrixType>
    std::set<std::size_t> GetZeroRows(const TMatrixType& rA, const ValueType tol) const
    {
        std::set<IndexType> zero_rows;

        for(IndexType i = 0; i < rA.size1(); ++i)
        {
            if( std::abs(rA(i, i)) < tol )
            {
                if( norm_2( row(rA, i) ) < tol )
                    zero_rows.insert(i);
            }
        }

        return zero_rows;
    }

    /**
     * Assignment operator.
     */
    DiagonalFitSolver& operator=(const DiagonalFitSolver& Other);

    /**
     * Copy constructor.
     */

}; // Class DiagonalFitSolver

} // namespace Kratos.

#endif // KRATOS_MULTITHREADED_SOLVERS_APPLICATION_DIAGONAL_FIT_SOLVER_H_INCLUDED  defined
