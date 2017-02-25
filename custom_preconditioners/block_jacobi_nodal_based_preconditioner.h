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
//   Date:                $Date: 29 Jun 2015 $
//   Revision:            $Revision: 1.3 $
//
//

#if !defined(KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BLOCK_JACOBI_NODAL_BASED_PRECONDITIONER_H_INCLUDED)
#define KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BLOCK_JACOBI_NODAL_BASED_PRECONDITIONER_H_INCLUDED


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

/// BlockJacobiNodalBasedPreconditioner class.
/**   */
template<class TSparseSpaceType, class TDenseSpaceType>
class BlockJacobiNodalBasedPreconditioner : public BlockJacobiPreconditioner<TSparseSpaceType, TDenseSpaceType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BlockJacobiNodalBasedPreconditioner
    KRATOS_CLASS_POINTER_DEFINITION (BlockJacobiNodalBasedPreconditioner);

    typedef Preconditioner<TSparseSpaceType, TDenseSpaceType> BaseType;

    typedef BlockJacobiPreconditioner<TSparseSpaceType, TDenseSpaceType> SuperType;
    
    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType> LinearSolverType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef std::size_t  SizeType;
    
    typedef std::size_t  IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    /// Provide the nodes to the preconditioner. The nodes is arranged in list of list. For example:
    /// nodes = [[1, 2, 3], [4, 5, 6]]
    BlockJacobiNodalBasedPreconditioner(boost::python::list& pyListNodes)
    {
        /* extract the nodal indices from python list */
        typedef boost::python::stl_input_iterator<boost::python::list> iterator_value_type;
        BOOST_FOREACH(const iterator_value_type::value_type& nodes_list, 
                      std::make_pair(iterator_value_type(pyListNodes), // begin
                      iterator_value_type() ) )                        // end
        {
            std::set<unsigned int> nodal_indices;
            typedef boost::python::stl_input_iterator<int> sub_iterator_value_type;
            BOOST_FOREACH(const sub_iterator_value_type::value_type& node_id, 
                      std::make_pair(sub_iterator_value_type(nodes_list), // begin
                      sub_iterator_value_type() ) )                       // end
            {
                nodal_indices.insert(node_id);
            }
            mNodalIndices.push_back(nodal_indices);
        }
    }

    /// Copy constructor.
    BlockJacobiNodalBasedPreconditioner(const BlockJacobiNodalBasedPreconditioner& rOther)
    : SuperType(rOther)
    {
    }

    /// Destructor.
    virtual ~BlockJacobiNodalBasedPreconditioner()
    {
    }


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    BlockJacobiNodalBasedPreconditioner& operator=(const BlockJacobiNodalBasedPreconditioner& rOther)
    {
        SuperType::operator=(rOther);
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
        /* build the map from node id to the block index */
        std::map<unsigned int, unsigned int> node_to_block_map;
        for(unsigned int block_id = 0; block_id < mNodalIndices.size(); ++block_id)
            for(std::set<unsigned int>::iterator it = mNodalIndices[block_id].begin(); it != mNodalIndices[block_id].end(); ++it)
                node_to_block_map[*it] = block_id;

        /* clear existing data */
        SuperType::mBlockIndices.clear();
        SuperType::mBlockDofs.clear();

        /* build the block indices */
        unsigned int tot_active_dofs = 0;
        unsigned int system_size = rA.size1();
        unsigned int node_id;
        unsigned int block_id;
        unsigned int num_blocks = mNodalIndices.size();
        SuperType::mBlockIndices.resize(num_blocks);
        SuperType::mBlockDofs.resize(num_blocks);
        
        std::vector<unsigned int> local_equation_id(num_blocks);
        for(unsigned int i = 0; i < num_blocks; ++i)
            local_equation_id[i] = 0;

        for(ModelPart::DofsArrayType::iterator it = rdof_set.begin(); it != rdof_set.end(); ++it)
            if(it->EquationId() < system_size)
            {
                ++tot_active_dofs;
                node_id = it->Id();
                block_id = node_to_block_map[node_id];
                SuperType::mBlockIndices[block_id].insert(it->EquationId());

                ModelPart::DofType ident_dof = *it;
                ident_dof.SetEquationId(local_equation_id[block_id]++);
                SuperType::mBlockDofs[block_id].push_back(ident_dof);
            }

        /* size check */
        if(tot_active_dofs != rA.size1() )
            KRATOS_THROW_ERROR(std::logic_error, "total system size does not equal total active dofs", "");
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

    virtual std::string Name() const {return "BlockJacobiNodalBasedPreconditioner";}

    /// Return information about this object.
    virtual std::string Info() const
    {
        return SuperType::Info();
    }


    /// Print information about this object.
    virtual void PrintInfo(std::ostream& OStream) const
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

    std::vector<std::set<unsigned int> > mNodalIndices;

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

}; // Class BlockJacobiNodalBasedPreconditioner

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
                                  BlockJacobiNodalBasedPreconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
{
    return IStream;
}


/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType>
inline std::ostream& operator << (std::ostream& OStream,
                                  const BlockJacobiNodalBasedPreconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
{
    rThis.PrintInfo(OStream);
    OStream << std::endl;
    rThis.PrintData(OStream);


    return OStream;
}
///@}


}  // namespace Kratos.


#endif // KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BLOCK_JACOBI_NODAL_BASED_PRECONDITIONER_H_INCLUDED  defined 

