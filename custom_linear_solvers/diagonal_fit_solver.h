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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// External includes 
#include <boost/timer.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/numeric/bindings/traits/sparse_traits.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_sparse.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"

#define CHECK_INCOMPATIBLE_ROWS

namespace ublas = boost::numeric::ublas;

namespace Kratos {

/*
This solver tries to change the diagonal before solving. If a row is zero, then the diagonal is set by the average value and the rhs is set to zero, regardless the previous value. It is essential to remove the singularity of the linear system. Especially when it's arised from methods like immersed boundary method.
Additionally, it also change the sign of the row if the negative diagonal is detected.
*/
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class DiagonalFitSolver: public DirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
{
public:
    /**
     * Counted pointer of SuperLUSolver
     */
    KRATOS_CLASS_POINTER_DEFINITION(DiagonalFitSolver);

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    /**
     * @param niter number of iterative refinements allowed
     */
    DiagonalFitSolver(typename BaseType::Pointer pLinearSolver) : mpLinearSolver(pLinearSolver) {}

    /**
     * Destructor
     */
    virtual ~DiagonalFitSolver() {}

    /** Some solvers may require a minimum degree of knowledge of the structure of the matrix. To make an example
     * when solving a mixed u-p problem, it is important to identify the row associated to v and p.
     * another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers
     * which require knowledge on the spatial position of the nodes associated to a given dof.
     * This function tells if the solver requires such data
     */
    virtual bool AdditionalPhysicalDataIsNeeded()
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
    virtual void ProvideAdditionalData(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB,
        typename ModelPart::DofsArrayType& rdof_set,
        ModelPart& r_model_part
    )
    {
        mpLinearSolver->ProvideAdditionalData(rA, rX, rB, rdof_set, r_model_part);

        #ifdef CHECK_INCOMPATIBLE_ROWS
        std::set<std::size_t> zero_rows = GetZeroRows(rA, 1.0e-12);
        std::set<std::size_t> incompatible_rows;
        std::set<std::size_t> incompatible_nodes;

        for(std::set<std::size_t>::iterator it = zero_rows.begin(); it != zero_rows.end(); ++it)
            if(rB(*it) != 0.0)
                incompatible_rows.insert(*it);

        std::cout << "incompatible_rows:";
        for(std::set<std::size_t>::iterator it = incompatible_rows.begin(); it != incompatible_rows.end(); ++it)
            std::cout << " " << *it;
        std::cout << std::endl;

//        std::cout << "dof set according to incompatible_rows:" << std::endl;
        for(std::set<std::size_t>::iterator it = incompatible_rows.begin(); it != incompatible_rows.end(); ++it)
        {
            for(typename ModelPart::DofsArrayType::iterator it2 = rdof_set.begin(); it2 != rdof_set.end(); ++it2)
            {
                if(it2->EquationId() == *it)
                {
//                    std::cout << "dof id " << it2->Id() << ", EquationId: " << it2->EquationId() << std::endl;
                    incompatible_nodes.insert(it2->Id());
                    break;
                }
            }
        }

        std::cout << "incompatible_nodes:";
        for(std::set<std::size_t>::iterator it = incompatible_nodes.begin(); it != incompatible_nodes.end(); ++it)
            std::cout << " " << *it;
        std::cout << std::endl;
        #endif

//        std::vector<std::size_t> sample_rows = {294, 295, 296, 300, 301, 302, 306, 307, 308, 312, 313, 314, 315, 316, 317, 318, 319, 320};
//        for(typename ModelPart::DofsArrayType::iterator it = rdof_set.begin(); it != rdof_set.end(); ++it)
//        {
//            std::vector<std::size_t>::iterator it2 = std::find(sample_rows.begin(), sample_rows.end(), it->EquationId());
//            if(it2 != sample_rows.end())
//            {
//                std::cout << "row " << *it2 << " node: " << it->Id() << std::endl;
//            }
//        }

//        std::vector<std::size_t> sample_nodes = {76450, 67455};
//        std::vector<std::size_t> sample_nodes = {1000, 7399, 4285};
//        for(typename ModelPart::DofsArrayType::iterator it = rdof_set.begin(); it != rdof_set.end(); ++it)
//        {
//            std::vector<std::size_t>::iterator it2 = std::find(sample_nodes.begin(), sample_nodes.end(), it->Id());
//            if(it2 != sample_nodes.end())
//            {
//                std::cout << "node " << it->Id() << " has row " << it->EquationId() << std::endl;
//            }
//        }
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
        std::size_t n = rA.size1();
        double diagonal_ave = 0.0;
        std::size_t diagonal_count = 0;
        std::size_t num_negative_rows = 0;
        std::set<std::size_t> zero_rows = GetZeroRows(rA, 1.0e-12);

//        std::vector<std::size_t> sample_rows = {96600, 96601, 96602, 106893, 106894, 106895};
//        std::vector<std::size_t> sample_rows = {6999, 21807, 7000, 7349, 21814, 7350};
//        std::vector<std::size_t> sample_rows = {2704, 2705, 2706, 14203, 14204, 14205, 9367, 9368, 9369};
//        for(std::size_t i = 0; i < sample_rows.size(); ++i)
//        {
//            const std::size_t& row = sample_rows[i];
//            std::cout << "row " << sample_rows[i] << " has lhs(d) = " << rA(row, row) << ", rhs = " << rB(row) << std::endl;
//        }

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
            rB(*it) = 0.0;
            std::cout << " " << *it;
        }
        std::cout << std::endl;

        // print the zero rows
        std::cout << zero_rows.size() << " zero rows detected; set diagonal to " << diagonal_ave
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
    bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
    {
        KRATOS_THROW_ERROR(std::logic_error, "ERROR: This solver can be used for single RHS only", "");
        return false;
    }

    /// Return information about this object.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "DiagonalFitSolver";
        return buffer.str();
    }
    
    /**
     * Print information about this object.
     */
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "DiagonalFitSolver finished.";
    }

    /**
     * Print object's data.
     */
    void PrintData(std::ostream& rOStream) const
    {
    }

private:

    typename BaseType::Pointer mpLinearSolver;

    template<typename TMatrixType>
    std::set<std::size_t> GetZeroRows(const TMatrixType& rA, const double& tol)
    {
        std::set<std::size_t> zero_rows;

        for(std::size_t i = 0; i < rA.size1(); ++i)
        {
            if( fabs(rA(i, i)) < tol )
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

};
// Class DiagonalFitSolver

/**
 * input stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::istream& operator >>(
        std::istream& rIStream,
        DiagonalFitSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis) {
    return rIStream;
}

/**
 * output stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator <<(
        std::ostream& rOStream,
        const DiagonalFitSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.

#undef CHECK_INCOMPATIBLE_ROWS

#endif // KRATOS_MULTITHREADED_SOLVERS_APPLICATION_DIAGONAL_FIT_SOLVER_H_INCLUDED  defined
