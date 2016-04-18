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
* Date:                $Date: 26 Oct 2014 $			                     *
* Revision:            $Revision: 1.4 $ 				                 *
*========================================================================*
* International Center of Numerical Methods in Engineering - CIMNE	     *
* Barcelona - Spain 							                         *
*========================================================================*
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


