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
//   Date:                $Date: 26 Aug 2014 $
//   Revision:            $Revision: 1.2 $
//
//

#if !defined(KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BLOCK_PRESSURE_INDEX_BASED_SCHUR_PRECONDITIONER_H_INCLUDED)
#define KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BLOCK_PRESSURE_INDEX_BASED_SCHUR_PRECONDITIONER_H_INCLUDED




// System includes



// External includes
#include <boost/smart_ptr.hpp>
#include <boost/numeric/ublas/vector.hpp>


// Project includes
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "linear_solvers/linear_solver.h"
#include "custom_preconditioners/block_2phase_schur_preconditioner.h"


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

///@name  Preconditioners
///@{

/// BlockPressureIndexBasedSchurPreconditioner class.
/**
REF: White, Borja
*/
template<class TSparseSpaceType, class TDenseSpaceType>
class BlockPressureIndexBasedSchurPreconditioner : public Block2PhaseSchurPreconditioner<TSparseSpaceType, TDenseSpaceType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BlockPressureIndexBasedSchurPreconditioner
    KRATOS_CLASS_POINTER_DEFINITION (BlockPressureIndexBasedSchurPreconditioner);

    typedef Block2PhaseSchurPreconditioner<TSparseSpaceType, TDenseSpaceType> BaseType;

    typedef Preconditioner<TSparseSpaceType, TDenseSpaceType> PreconditionerType;

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType> LinearSolverType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef std::size_t  SizeType;

    typedef std::size_t  IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BlockPressureIndexBasedSchurPreconditioner(
        typename PreconditionerType::Pointer prec_A,
        typename PreconditionerType::Pointer prec_S,
        const std::string& SchurComputeMode,
        const std::size_t& block_size1,
        const std::size_t& block_size2)
    : BaseType(prec_A, prec_S, SchurComputeMode)
    , mblock_size1(block_size1), mblock_size2(block_size2)
    {}

    BlockPressureIndexBasedSchurPreconditioner(
        typename PreconditionerType::Pointer prec_A,
        typename PreconditionerType::Pointer prec_S,
        const std::string& SchurComputeMode,
        const std::string& InverseOption,
        const std::size_t& block_size1,
        const std::size_t& block_size2)
    : BaseType(prec_A, prec_S, SchurComputeMode, InverseOption)
    , mblock_size1(block_size1), mblock_size2(block_size2)
    {}

    BlockPressureIndexBasedSchurPreconditioner(
        typename PreconditionerType::Pointer prec_A,
        typename PreconditionerType::Pointer prec_S,
        const std::string& SchurComputeMode,
        typename LinearSolverType::Pointer solver_S,
        const std::size_t& block_size1,
        const std::size_t& block_size2)
    : BaseType(prec_A, prec_S, SchurComputeMode, solver_S)
    , mblock_size1(block_size1), mblock_size2(block_size2)
    {}

    BlockPressureIndexBasedSchurPreconditioner(
        typename PreconditionerType::Pointer prec_A,
        typename PreconditionerType::Pointer prec_S,
        const std::string& SchurComputeMode,
        typename LinearSolverType::Pointer solver_S,
        const std::string& InverseOption,
        const std::size_t& block_size1,
        const std::size_t& block_size2)
    : BaseType(prec_A, prec_S, SchurComputeMode, solver_S, InverseOption)
    , mblock_size1(block_size1), mblock_size2(block_size2)
    {}

    /// Copy constructor.
    BlockPressureIndexBasedSchurPreconditioner(const BlockPressureIndexBasedSchurPreconditioner& Other)
    : BaseType(Other), mblock_size1(Other.mblock_size1), mblock_size2(Other.mblock_size2)
    {}


    /// Destructor.
    virtual ~BlockPressureIndexBasedSchurPreconditioner()
    {
    }


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    BlockPressureIndexBasedSchurPreconditioner& operator=(const BlockPressureIndexBasedSchurPreconditioner& Other)
    {
        BaseType::operator=(Other);
        mblock_size1 = Other.mblock_size1;
        mblock_size2 = Other.mblock_size2;
        return *this;
    }


    ///@}
    ///@name Operations
    ///@{

    virtual void Initialize(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
        //count pressure dofs
        unsigned int n_nodes = rA.size1() / (mblock_size1 + mblock_size2);
        unsigned int n_pressure_dofs = n_nodes * mblock_size2;
        unsigned int tot_active_dofs = rA.size1();
        unsigned int system_size = TSparseSpaceType::Size(rB);

        KRATOS_WATCH(tot_active_dofs)
        KRATOS_WATCH(n_pressure_dofs)

        //resize arrays as needed
        unsigned int other_dof_size = tot_active_dofs - n_pressure_dofs;
        BaseType::mpressure_indices.resize(n_pressure_dofs, false);
        BaseType::mother_indices.resize(other_dof_size, false);
        BaseType::mglobal_to_local_indexing.resize(tot_active_dofs, false);
        BaseType::mis_pressure_block.resize(tot_active_dofs, false);
        //construct aux_lists as needed
        //"other_counter[i]" i will contain the position in the global system of the i-th NON-pressure node
        //"pressure_counter[i]" will contain the in the global system of the i-th pressure node
        //
        //mglobal_to_local_indexing[i] will contain the position in the local blocks of the
        unsigned int pressure_counter = 0;
        unsigned int other_counter = 0;
        unsigned int global_pos;
        for (unsigned int i = 0; i < tot_active_dofs; ++i)
        {
            global_pos = i;
            if ( global_pos % (mblock_size1 + mblock_size2) >= mblock_size1 )
            {
                BaseType::mpressure_indices[pressure_counter] = global_pos;
                BaseType::mglobal_to_local_indexing[global_pos] = pressure_counter;
                BaseType::mis_pressure_block[global_pos] = true;
                ++pressure_counter;
            }
            else
            {
                BaseType::mother_indices[other_counter] = global_pos;
                BaseType::mglobal_to_local_indexing[global_pos] = other_counter;
                BaseType::mis_pressure_block[global_pos] = false;
                ++other_counter;
            }
        }

        BaseType::Initialize(rA, rX, rB);
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
        buffer << BaseType::Info() << " for displacement-pressure system";
        return buffer.str();
    }


    /// Print information about this object.
    virtual void  PrintInfo(std::ostream& OStream) const
    {
        OStream << Info();
    }


    virtual void PrintData(std::ostream& OStream) const
    {
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

    std::size_t mblock_size1;
    std::size_t mblock_size2;

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
    ///@name Unaccessible methods
    ///@{


    ///@}

}; // Class BlockPressureIndexBasedSchurPreconditioner

///@}


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType>
inline std::istream& operator >> (std::istream& IStream,
                                  BlockPressureIndexBasedSchurPreconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
{
    return IStream;
}


/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType>
inline std::ostream& operator << (std::ostream& OStream,
                                  const BlockPressureIndexBasedSchurPreconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
{
    rThis.PrintInfo(OStream);
    OStream << std::endl;
    rThis.PrintData(OStream);


    return OStream;
}
///@}


}  // namespace Kratos.


#endif // KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BLOCK_PRESSURE_INDEX_BASED_SCHUR_PRECONDITIONER_H_INCLUDED  defined

