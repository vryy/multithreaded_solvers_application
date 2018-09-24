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
//   Date:                $Date: 7 Nov 2014 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BLOCK_2_PHASE_INDEX_BASED_SCHUR_SOLVER_H_INCLUDED )
#define KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BLOCK_2_PHASE_INDEX_BASED_SCHUR_SOLVER_H_INCLUDED


// System includes
#include <cmath>


// External includes
#include "boost/smart_ptr.hpp"


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
         class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class Block2PhaseIndexBasedSchurSolver : public Block2PhaseSchurSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of  Block2PhaseIndexBasedSchurSolver
    KRATOS_CLASS_POINTER_DEFINITION( Block2PhaseIndexBasedSchurSolver );

    typedef Block2PhaseSchurSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename BaseType::LinearSolverType LinearSolverType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
    
    typedef std::size_t  SizeType;
    
    typedef std::size_t  IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Block2PhaseIndexBasedSchurSolver() : BaseType() {}

    Block2PhaseIndexBasedSchurSolver(
        typename LinearSolverType::Pointer pSolver
    ) : BaseType(pSolver), mblock_size1(3), mblock_size2(1)
    {}

    Block2PhaseIndexBasedSchurSolver(
        typename LinearSolverType::Pointer pSolver,
        const unsigned int& block_size1,
        const unsigned int& block_size2
    ) : BaseType(pSolver), mblock_size1(block_size1), mblock_size2(block_size2)
    {}

    /// Copy constructor.
     Block2PhaseIndexBasedSchurSolver(const  Block2PhaseIndexBasedSchurSolver& Other) : BaseType(Other) {}

    /// Destructor.
    virtual ~ Block2PhaseIndexBasedSchurSolver() {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
     Block2PhaseIndexBasedSchurSolver& operator=(const  Block2PhaseIndexBasedSchurSolver& Other)
    {
        BaseType::operator=(Other);
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
        buffer << "Linear solver using Schur complement reduction scheme for displacement-pressure" << std::endl;
        buffer << BaseType::Info();
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

    unsigned int mblock_size1;
    unsigned int mblock_size2;

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

}; // Class  Block2PhaseIndexBasedSchurSolver

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::istream& operator >> (std::istream& IStream, Block2PhaseIndexBasedSchurSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis)
{
    return IStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator << (std::ostream& OStream, const  Block2PhaseIndexBasedSchurSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis)
{
    rThis.PrintInfo(OStream);
    OStream << std::endl;
    rThis.PrintData(OStream);

    return OStream;
}
///@}


}  // namespace Kratos.

#undef CHECK_DIAGONAL_DOMINANCE
#undef STRINGIFY

#endif //  KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BLOCK_2_PHASE_INDEX_BASED_SCHUR_SOLVER_H_INCLUDED  defined 

