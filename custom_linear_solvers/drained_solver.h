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
//   Date:                $Date: 24 Jul 2020 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_MULTITHREADED_SOLVERS_APPLICATION_DRAINED_SOLVER_H_INCLUDED )
#define KRATOS_MULTITHREADED_SOLVERS_APPLICATION_DRAINED_SOLVER_H_INCLUDED


// System includes
#include <cmath>


// External includes


// Project includes
#include "includes/define.h"
#include "linear_solvers/linear_solver.h"
#include "custom_utilities/multithreaded_solvers_math_utils.h"


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

/// Short class definition.
/** Detail class definition.
THis solver solve the linear system from drained analysis by only solve the u-part and leave the p-part untouched
*/
template<class TSparseSpaceType, class TDenseSpaceType,
         class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class DrainedSolver : public LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of  DrainedSolver
    KRATOS_CLASS_POINTER_DEFINITION( DrainedSolver );

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> LinearSolverType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef std::size_t  SizeType;

    typedef std::size_t  IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DrainedSolver() {}

    DrainedSolver(typename BaseType::Pointer pSolver)
    : BaseType(), mpSolver(pSolver)
    {}

    /// Copy constructor.
    DrainedSolver(const  DrainedSolver& Other)
    : BaseType(Other) {}

    /// Destructor.
    virtual ~DrainedSolver() {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    DrainedSolver& operator=(const  DrainedSolver& Other)
    {
        BaseType::operator=(Other);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    virtual bool AdditionalPhysicalDataIsNeeded()
    {
        return true;
    }

    virtual void ProvideAdditionalData(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB,
        typename ModelPart::DofsArrayType& rdof_set,
        ModelPart& r_model_part
    )
    {
        //count pressure dofs
        unsigned int n_pressure_dofs = 0;
        unsigned int tot_active_dofs = 0;
        unsigned int system_size = TSparseSpaceType::Size(rB);
        for (ModelPart::DofsArrayType::iterator it = rdof_set.begin(); it != rdof_set.end(); ++it)
            if (it->EquationId() < system_size)
            {
                ++tot_active_dofs;
                if ( (it->GetVariable().Key() == WATER_PRESSURE) || (it->GetVariable().Key() == AIR_PRESSURE) )
                    ++n_pressure_dofs;
            }
        if (tot_active_dofs != rA.size1() )
            KRATOS_THROW_ERROR (std::logic_error,"total system size does not coincide with the free dof map","");

        KRATOS_WATCH(tot_active_dofs)
        KRATOS_WATCH(n_pressure_dofs)

        //resize arrays as needed
        unsigned int other_dof_size = tot_active_dofs - n_pressure_dofs;
        mpressure_indices.resize(n_pressure_dofs, false);
        mother_indices.resize(other_dof_size, false);
        mglobal_to_local_indexing.resize(tot_active_dofs, false);
        mis_pressure_block.resize(tot_active_dofs, false);
        //construct aux_lists as needed
        //"other_counter[i]" i will contain the position in the global system of the i-th NON-pressure node
        //"pressure_counter[i]" will contain the in the global system of the i-th NON-pressure node
        //
        //mglobal_to_local_indexing[i] will contain the position in the local blocks of the
        unsigned int pressure_counter = 0;
        unsigned int other_counter = 0;
        unsigned int global_pos;
        madof_set.clear();
        for (ModelPart::DofsArrayType::iterator it = rdof_set.begin(); it != rdof_set.end(); ++it)
        {
            global_pos = it->EquationId();
            if (global_pos < system_size)
            {
                if ( (it->GetVariable().Key() == WATER_PRESSURE) || (it->GetVariable().Key() == AIR_PRESSURE) )
                {
                    mpressure_indices[pressure_counter] = global_pos;
                    mglobal_to_local_indexing[global_pos] = pressure_counter;
                    mis_pressure_block[global_pos] = true;
                    ++pressure_counter;
                }
                else
                {
                    mother_indices[other_counter] = global_pos;
                    mglobal_to_local_indexing[global_pos] = other_counter;
                    mis_pressure_block[global_pos] = false;
                    madof_set.push_back(*it);
                    ++other_counter;
                }
            }
        }

        // madof_set.Sort();
        KRATOS_WATCH(madof_set.size())

        mpSolver->ProvideAdditionalData(rA, rX, rB, madof_set, r_model_part);
    }

    /** Normal solve method.
    Solves the linear system Ax=b and puts the result on SystemVector& rX.
    rX is also th initial guess for iterative methods.
    @param rA. System matrix
    @param rX. Solution vector. it's also the initial
    guess for iterative linear solvers.
    @param rB. Right hand side vector.
    */
    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
        if(this->IsNotConsistent(rA, rX, rB))
            return false;

        KRATOS_WATCH(rA.size1())
        KRATOS_WATCH(rA.size2())

        std::cout << "Extract non-pressure block begin" << std::endl;
        double start = OpenMPUtils::GetCurrentTime();
        // SparseMatrixType A;
        // MultithreadedSolversMathUtils::ExtractBlock(rA,
        //     mother_indices, mglobal_to_local_indexing, mis_other_block,
        //     A);
        SparseMatrixType A;
        MultithreadedSolversMathUtils::FillBlockMatrices(rA,
            mother_indices, mpressure_indices, mglobal_to_local_indexing, mis_pressure_block,
            A);
        std::cout << "Extract non-pressure block completed..." << OpenMPUtils::GetCurrentTime() - start << " s" << std::endl;

        KRATOS_WATCH(MultithreadedSolversMathUtils::ComputeFrobeniusNorm(A))
        KRATOS_WATCH(MultithreadedSolversMathUtils::DiagonalNorm(A))

        // TODO: check the diagonal dominance of matrix A
        #ifdef CHECK_DIAGONAL_DOMINANCE
        MultithreadedSolversMathUtils::CheckDiagonalDominance(A);
        #endif

        VectorType ru, u, p;
        TSparseSpaceType::Resize(ru, mother_indices.size());
        TSparseSpaceType::Resize(u, mother_indices.size());
        TSparseSpaceType::Resize(p, mpressure_indices.size());

        // Get the initial u
        MultithreadedSolversMathUtils::GetPart(rX, mother_indices, u);

        // Get ru
        MultithreadedSolversMathUtils::GetPart(rB, mother_indices, ru);
        KRATOS_WATCH(norm_2(ru))

        mpSolver->Initialize(A, u, ru);

        mpSolver->Solve(A, u, ru);
//        KRATOS_WATCH(u)
        KRATOS_WATCH(norm_2(u))
        std::cout << "Solve for u completed" << std::endl;

        // Write back the result to solution vector
        MultithreadedSolversMathUtils::WritePart(rX, mother_indices, u);

        // Write back the result to solution vector
        TSparseSpaceType::SetToZero(p);
        MultithreadedSolversMathUtils::WritePart(rX, mpressure_indices, p);

        return true;
    }

    /**
    Multi solve method for solving a set of linear systems with same coefficient matrix.
    Solves the linear system Ax=b and puts the result on SystemVector& rX.
    rX is also th initial guess for iterative methods.
    @param rA. System matrix
    @param rX. Solution vector. it's also the initial
    guess for iterative linear solvers.
    @param rB. Right hand side vector.
    */
    bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Multisolve is not yet supported for", typeid(*this).name())
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
        buffer << "Linear solver for drained analysis (solid solver: " << mpSolver->Info() << ")";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& OStream) const
    {
        OStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& OStream) const
    {
        BaseType::PrintData(OStream);
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

    typename BaseType::Pointer mpSolver;

    std::vector<SizeType> mpressure_indices;
    std::vector<SizeType> mother_indices;
    std::vector<SizeType> mglobal_to_local_indexing;
    std::vector<int> mis_pressure_block;

    typename ModelPart::DofsArrayType madof_set;

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


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

}; // Class  DrainedSolver

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::istream& operator >> (std::istream& IStream, DrainedSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis)
{
    return IStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator << (std::ostream& OStream, const  DrainedSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis)
{
    rThis.PrintInfo(OStream);
    OStream << std::endl;
    rThis.PrintData(OStream);

    return OStream;
}
///@}


}  // namespace Kratos.

#undef CHECK_DIAGONAL_DOMINANCE

#endif //  KRATOS_MULTITHREADED_SOLVERS_APPLICATION_DRAINED_SOLVER_H_INCLUDED  defined

