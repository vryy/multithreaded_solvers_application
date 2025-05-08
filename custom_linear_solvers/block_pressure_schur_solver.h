/*
see multithreaded_solvers_application/LICENSE.txt
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 7 Nov 2014 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BLOCK_PRESSURE_SCHUR_SOLVER_H_INCLUDED )
#define KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BLOCK_PRESSURE_SHUR_SOLVER_H_INCLUDED


// System includes
#include <cmath>


// External includes


// Project includes
#include "includes/define.h"
#include "linear_solvers/linear_solver.h"
#include "custom_linear_solvers/block_2phase_schur_solver.h"


#define STRINGIFY(name) #name
// #define CHECK_DIAGONAL_DOMINANCE


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
THis solver solve the coupled u-p problem by using a staggerred scheme:
[A  B1][u] = [ru]
[B2  C][p]   [rp]

=> u = A^(-1) * (ru -B1 * p)                                (1)
   (C - B2 * A^(-1) * B1) * p = (rp - B2 * A^(-1) * ru)     (2)

Firstly p in (2) is solved approximately by approximating A^(-1) by diag(A)^(-1) (remarks: A needs to be diagonal-dominant for this to be a good approximation)
Then u in (1) is solved exactly
*/
template<class TSparseSpaceType, class TDenseSpaceType,
         class TModelPartType,
         class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class BlockPressureSchurSolver : public Block2PhaseSchurSolver<TSparseSpaceType, TDenseSpaceType, TModelPartType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of  BlockPressureSchurSolver
    KRATOS_CLASS_POINTER_DEFINITION( BlockPressureSchurSolver );

    typedef Block2PhaseSchurSolver<TSparseSpaceType, TDenseSpaceType, TModelPartType, TReordererType> BaseType;

    typedef typename BaseType::LinearSolverType LinearSolverType;

    typedef typename BaseType::SparseMatrixType SparseMatrixType;

    typedef typename BaseType::VectorType VectorType;

    typedef typename BaseType::DenseMatrixType DenseMatrixType;

    typedef typename BaseType::SizeType SizeType;

    typedef typename BaseType::IndexType IndexType;

    typedef typename BaseType::DataType DataType;

    typedef typename BaseType::ValueType ValueType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BlockPressureSchurSolver() : BaseType() {}

    BlockPressureSchurSolver(typename LinearSolverType::Pointer pSolver)
    : BaseType(pSolver)
    {}

    BlockPressureSchurSolver(typename LinearSolverType::Pointer pSolverA,
        typename LinearSolverType::Pointer pSolverS)
    : BaseType(pSolverA, pSolverS)
    {}

    /// Copy constructor.
     BlockPressureSchurSolver(const  BlockPressureSchurSolver& Other)
     : BaseType(Other) {}

    /// Destructor.
    ~BlockPressureSchurSolver() override {}

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    BlockPressureSchurSolver& operator=(const  BlockPressureSchurSolver& Other)
    {
        BaseType::operator=(Other);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    bool AdditionalPhysicalDataIsNeeded() override
    {
        return true;
    }

    void ProvideAdditionalData(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB,
        typename TModelPartType::DofsArrayType& rdof_set,
        TModelPartType& r_model_part
    ) override
    {
        //count pressure dofs
        unsigned int n_pressure_dofs = 0;
        unsigned int tot_active_dofs = 0;
        unsigned int system_size = TSparseSpaceType::Size(rB);
        for (auto it = rdof_set.begin(); it != rdof_set.end(); ++it)
            if (it->EquationId() < system_size)
            {
                ++tot_active_dofs;
                if ( (it->GetVariable().Key() == WATER_PRESSURE) || (it->GetVariable().Key() == AIR_PRESSURE) )
                    ++n_pressure_dofs;
            }
        if (tot_active_dofs != rA.size1() )
            KRATOS_THROW_ERROR (std::logic_error,"total system size does not coincide with the free dof map","");

        // KRATOS_WATCH(tot_active_dofs)
        // KRATOS_WATCH(n_pressure_dofs)

        //resize arrays as needed
        unsigned int other_dof_size = tot_active_dofs - n_pressure_dofs;
        BaseType::mpressure_indices.resize(n_pressure_dofs, false);
        BaseType::mother_indices.resize(other_dof_size, false);
        BaseType::mglobal_to_local_indexing.resize(tot_active_dofs, false);
        BaseType::mis_pressure_block.resize(tot_active_dofs, false);
        //construct aux_lists as needed
        //"other_counter[i]" i will contain the position in the global system of the i-th NON-pressure node
        //"pressure_counter[i]" will contain the in the global system of the i-th NON-pressure node
        //
        //mglobal_to_local_indexing[i] will contain the position in the local blocks of the
        unsigned int pressure_counter = 0;
        unsigned int other_counter = 0;
        unsigned int global_pos;
        BaseType::madof_set.clear();
        BaseType::msdof_set.clear();
        for (auto it = rdof_set.begin(); it != rdof_set.end(); ++it)
        {
            global_pos = it->EquationId();
            if (global_pos < system_size)
            {
                if ( (it->GetVariable().Key() == WATER_PRESSURE) || (it->GetVariable().Key() == AIR_PRESSURE) )
                {
                    BaseType::mpressure_indices[pressure_counter] = global_pos;
                    BaseType::mglobal_to_local_indexing[global_pos] = pressure_counter;
                    BaseType::mis_pressure_block[global_pos] = true;
                    BaseType::msdof_set.push_back(*it);
                    ++pressure_counter;
                }
                else
                {
                    BaseType::mother_indices[other_counter] = global_pos;
                    BaseType::mglobal_to_local_indexing[global_pos] = other_counter;
                    BaseType::mis_pressure_block[global_pos] = false;
                    BaseType::madof_set.push_back(*it);
                    ++other_counter;
                }
            }
        }

        SparseMatrixType DummyA;
        VectorType DummyX, DummyB;
        BaseType::mpSolverS->ProvideAdditionalData(DummyA, DummyX, DummyB, BaseType::msdof_set, r_model_part);
        BaseType::mpSolverA->ProvideAdditionalData(DummyA, DummyX, DummyB, BaseType::madof_set, r_model_part);
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
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Linear solver using Schur complement reduction scheme for displacement-pressure" << std::endl;
        buffer << BaseType::Info();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& OStream) const override
    {
        OStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& OStream) const override
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
    ///@name Un accessible methods
    ///@{


    ///@}

}; // Class  BlockPressureSchurSolver

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

}  // namespace Kratos.

#undef CHECK_DIAGONAL_DOMINANCE
#undef STRINGIFY

#endif //  KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BLOCK_PRESSURE_SHUR_SOLVER_H_INCLUDED  defined
