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
//   Date:                $Date: 3 July 2015 $
//   Revision:            $Revision: 1.3 $
//
//

#if !defined(KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BLOCK_JACOBI_NODAL_BASED_PRESSURE_PRECONDITIONER_H_INCLUDED)
#define KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BLOCK_JACOBI_NODAL_BASED_PRESSURE_PRECONDITIONER_H_INCLUDED


// System includes



// External includes
#include <boost/foreach.hpp>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>

// Project includes
#include "includes/define.h"
#include "linear_solvers/linear_solver.h"
#include "utilities/openmp_utils.h"
#include "block_jacobi_preconditioner.h"

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

/// BlockJacobiNodalBasedPressurePreconditioner class.
/**
THis is similar to BlockJacobiNodalBasedPreconditioner but introduce pressure to a separate block
*/
template<class TSparseSpaceType, class TDenseSpaceType, class TModelPartType>
class BlockJacobiNodalBasedPressurePreconditioner : public BlockJacobiPreconditioner<TSparseSpaceType, TDenseSpaceType, TModelPartType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BlockJacobiNodalBasedPressurePreconditioner
    KRATOS_CLASS_POINTER_DEFINITION (BlockJacobiNodalBasedPressurePreconditioner);

    typedef Preconditioner<TSparseSpaceType, TDenseSpaceType, TModelPartType> BaseType;

    typedef BlockJacobiPreconditioner<TSparseSpaceType, TDenseSpaceType, TModelPartType> SuperType;

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TModelPartType> LinearSolverType;

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

    /// Default constructor
    /// Provide the nodes to the preconditioner. The nodes is arranged in list of list. For example:
    /// nodes = [[1, 2, 3], [4, 5, 6]]
    BlockJacobiNodalBasedPressurePreconditioner(boost::python::list& pyListNodes)
    {
        /* extract the nodal indices from python list */
        typedef boost::python::stl_input_iterator<boost::python::list> iterator_value_type;
        BOOST_FOREACH(const iterator_value_type::value_type& nodes_list,
                      std::make_pair(iterator_value_type(pyListNodes), // begin
                      iterator_value_type() ) )                        // end
        {
            std::set<IndexType> nodal_indices;
            typedef boost::python::stl_input_iterator<int> sub_iterator_value_type;
            BOOST_FOREACH(const sub_iterator_value_type::value_type& node_id,
                      std::make_pair(sub_iterator_value_type(nodes_list), // begin
                      sub_iterator_value_type() ) )                       // end
            {
                nodal_indices.insert(static_cast<IndexType>(node_id));
            }
            mNodalIndices.push_back(nodal_indices);
        }
    }

    /// Copy constructor.
    BlockJacobiNodalBasedPressurePreconditioner(const BlockJacobiNodalBasedPressurePreconditioner& rOther)
    : SuperType(rOther)
    {
    }

    /// Destructor.
    ~BlockJacobiNodalBasedPressurePreconditioner() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    BlockJacobiNodalBasedPressurePreconditioner& operator=(const BlockJacobiNodalBasedPressurePreconditioner& rOther)
    {
        SuperType::operator=(rOther);
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
        /* build the map from node id to the block index */
        std::map<unsigned int, unsigned int> node_to_block_map;
        for(unsigned int block_id = 0; block_id < mNodalIndices.size(); ++block_id)
            for(auto it = mNodalIndices[block_id].begin(); it != mNodalIndices[block_id].end(); ++it)
                node_to_block_map[*it] = block_id;

        /* clear existing data */
        SuperType::mBlockIndices.clear();
        SuperType::mBlockDofs.clear();

        /* build the block indices */
        unsigned int tot_active_dofs = 0;
        unsigned int system_size = rA.size1();
        unsigned int node_id;
        unsigned int block_id;
        unsigned int num_blocks = mNodalIndices.size() + 1;
        SuperType::mBlockIndices.resize(num_blocks);
        SuperType::mBlockDofs.resize(num_blocks);

        std::vector<unsigned int> local_equation_id(num_blocks);
        for(unsigned int i = 0; i < num_blocks; ++i)
            local_equation_id[i] = 0;

        for(auto it = rdof_set.begin(); it != rdof_set.end(); ++it)
            if(it->EquationId() < system_size)
            {
                ++tot_active_dofs;
                node_id = it->Id();
                if((it->GetVariable().Key() == WATER_PRESSURE) || (it->GetVariable().Key() == AIR_PRESSURE))
                    block_id = num_blocks - 1;
                else
                    block_id = node_to_block_map[node_id];
                SuperType::mBlockIndices[block_id].insert(it->EquationId());

                typename TModelPartType::DofType ident_dof = *it;
                ident_dof.SetEquationId(local_equation_id[block_id]++);
                SuperType::mBlockDofs[block_id].push_back(ident_dof);
            }

        /* size check */
        if(tot_active_dofs != rA.size1() )
            KRATOS_ERROR << "total system size does not equal total active dofs";
        std::cout << "Build the block indices completed, tot_active_dofs = " << tot_active_dofs << std::endl;

        /* provide additional data for each sub-preconditioner */
        SuperType::ProvideAdditionalData(rA, rX, rB, rdof_set, r_model_part);
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

    virtual std::string Name() const {return "BlockJacobiNodalBasedPressurePreconditioner";}

    /// Return information about this object.
    std::string Info() const override
    {
        return SuperType::Info();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& OStream) const override
    {
        OStream << Info();
    }

    void PrintData(std::ostream& OStream) const override
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

    std::vector<std::set<IndexType> > mNodalIndices;

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

}; // Class BlockJacobiNodalBasedPressurePreconditioner

///@}


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

}  // namespace Kratos.

#endif // KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BLOCK_JACOBI_NODAL_BASED_PRECONDITIONER_H_INCLUDED  defined
