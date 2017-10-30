/*
see multithreaded_solvers_application/LICENSE.txt
*/

#if !defined(KRATOS_MULTITHREADED_SOLVERS_APPLICATION_VARIABLE_SOLVER_H_INCLUDED )
#define  KRATOS_MULTITHREADED_SOLVERS_APPLICATION_VARIABLE_SOLVER_H_INCLUDED

// #define BOOST_NUMERIC_BINDINGS_SUPERLU_PRINT

// System includes
#include <stdio.h>
#include <cstdlib>
#include <cmath>

// External includes


// Project includes
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "includes/ublas_interface.h"
#include "linear_solvers/linear_solver.h"

#define ENABLE_PROFILING

namespace ublas = boost::numeric::ublas;

namespace Kratos
{
template< class TSparseSpaceType, class TDenseSpaceType,
          class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class VariableSolver : public DirectSolver< TSparseSpaceType, TDenseSpaceType, TReordererType>
{
public:
    /**
     * Counted pointer
     */
    typedef boost::shared_ptr<VariableSolver> Pointer;

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    VariableSolver()
    {
    }

    /**
     * Destructor
     */
    virtual ~VariableSolver() {}

    void AddSolver(typename BaseType::Pointer pSolver)
    {
        mSolverList.push_back(pSolver);
    }
    
    virtual bool AdditionalPhysicalDataIsNeeded()
    {
        bool is_needed = false;
        for(unsigned int i = 0; i < mSolverList.size(); ++i)
            is_needed = is_needed || (mSolverList[i]->AdditionalPhysicalDataIsNeeded());
        return is_needed;
    }

    virtual void ProvideAdditionalData(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB,
        typename ModelPart::DofsArrayType& rdof_set,
        ModelPart& r_model_part
    )
    {
        if(AdditionalPhysicalDataIsNeeded())
        {
            for(unsigned int i = 0; i < mSolverList.size(); ++i)
            {
                if(mSolverList[i]->AdditionalPhysicalDataIsNeeded())
                    mSolverList[i]->ProvideAdditionalData(rA, rX, rB, rdof_set, r_model_part);
            }
        }
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
        if(mSolverList.size() == 0)
            KRATOS_THROW_ERROR(std::logic_error, "Error: The variable solver does not contain any solver. Please add a solver for it to work", "")
        else
        {
            int cnt = 0;
            while(cnt < mSolverList.size())
            {
                bool isSolved = mSolverList[cnt]->Solve(rA, rX, rB);
                
                if(!isSolved)
                {
                    std::cout << "The solver " << *(mSolverList[cnt]) << " does not converge, switch to the next solver" << std::endl;
                    ++cnt;
                }
                else
                {
                    std::cout << "The solver " << *(mSolverList[cnt]) << " completed" << std::endl;
                    break;
                }
            }
            if(cnt == mSolverList.size())
                KRATOS_THROW_ERROR(std::logic_error, "None of the provided solver converged for the current problem. Consider adding more solvers", "")
        }
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
    }

    /// Return information about this object.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Variable solver";
        return buffer.str();
    }
    
    /**
     * Print information about this object.
     */
    void  PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Variable solver finished.";
    }

    /**
     * Print object's data.
     */
    void  PrintData(std::ostream& rOStream) const
    {
    }

private:

    std::vector<typename BaseType::Pointer> mSolverList;

    /**
     * Assignment operator.
     */
    VariableSolver& operator=(const VariableSolver& Other);

    /**
     * Copy constructor.
     */
//    VariableSolver(const ParallelSuperLUSolver& Other);

}; // Class VariableSolver

}  // namespace Kratos.

#endif // KRATOS_MULTITHREADED_SOLVERS_APPLICATION_VARIABLE_SOLVER_H_INCLUDED  defined 


