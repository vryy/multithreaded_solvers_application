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
//   Date:                $Date: 1 July 2015 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(MULTITHREADED_SOLVERS_APP_DEFLATED_SUBDOMAIN_NODAL_BASED_CG_SOLVER_H_INCLUDED )
#define  MULTITHREADED_SOLVERS_APP_DEFLATED_SUBDOMAIN_NODAL_BASED_CG_SOLVER_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <vector>
#include <set>

// External includes
#include <boost/foreach.hpp>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>

// Project includes
#include "includes/define.h"
#include "linear_solvers/iterative_solver.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"
#include "utilities/deflation_utils.h"
#include "deflated_cg_solver_2.h"

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
 */
template<class TSparseSpaceType, class TDenseSpaceType,
         class TModelPartType,
         class TPreconditionerType = Preconditioner<TSparseSpaceType, TDenseSpaceType, TModelPartType>,
         class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class DeflatedSubdomainNodalBasedCGSolver : public DeflatedCGSolver2<TSparseSpaceType, TDenseSpaceType, TModelPartType, TPreconditionerType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DeflatedSubdomainNodalBasedCGSolver
    KRATOS_CLASS_POINTER_DEFINITION(DeflatedSubdomainNodalBasedCGSolver);

    typedef DeflatedCGSolver2<TSparseSpaceType, TDenseSpaceType, TModelPartType, TPreconditionerType, TReordererType> BaseType;

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TModelPartType, TReordererType> LinearSolverType;

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
    DeflatedSubdomainNodalBasedCGSolver(ValueType NewMaxTolerance, unsigned int NewMaxIterationsNumber, typename TPreconditionerType::Pointer pNewPreconditioner, boost::python::list& pyListNodes)
    : BaseType(NewMaxTolerance, NewMaxIterationsNumber, pNewPreconditioner)
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
    DeflatedSubdomainNodalBasedCGSolver(const DeflatedSubdomainNodalBasedCGSolver& Other) : BaseType(Other)
    {
    }

    /// Destructor.
    ~DeflatedSubdomainNodalBasedCGSolver() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.

    DeflatedSubdomainNodalBasedCGSolver & operator=(const DeflatedSubdomainNodalBasedCGSolver& Other)
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
        /* clear the internal data */
        mw.clear();

        /* build the map from node id to the block index */
        std::map<unsigned int, unsigned int> node_to_block_map;
        for(unsigned int block_id = 0; block_id < mNodalIndices.size(); ++block_id)
            for(auto it = mNodalIndices[block_id].begin(); it != mNodalIndices[block_id].end(); ++it)
                node_to_block_map[*it] = block_id;

        /* populate the deflation space */
        unsigned int system_size = rA.size1();
        mw.resize(system_size, 0);
        unsigned int tot_active_dofs = 0;
        unsigned int node_id, block_id;
        for(auto it = rdof_set.begin(); it != rdof_set.end(); ++it)
            if(it->EquationId() < system_size)
            {
                ++tot_active_dofs;
                node_id = it->Id();
                block_id = node_to_block_map[node_id];
                mw[it->EquationId()] = block_id;
            }

        /* size check */
        if(tot_active_dofs != rA.size1())
            KRATOS_THROW_ERROR(std::logic_error, "total system size does not equal total active dofs", "");

        /* call super class ProvideAdditionalData */
        BaseType::ProvideAdditionalData(rA, rX, rB, rdof_set, r_model_part);
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

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Deflated subdomain nodal-based Conjugate gradient linear solver with " << BaseType::GetPreconditioner()->Info();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
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

    void ConstructW(SparseMatrixType& rA, std::vector<int>& w, SparseMatrixType& deflatedA) const override
    {
        /* export the deflation space */
        w.resize(mw.size());
        std::copy(mw.begin(), mw.end(), w.begin());

        SizeType full_size = rA.size1();
        SizeType reduced_size = mNodalIndices.size();

        /* construct non-zero structure of deflation matrix */
        // Non-zero structure of deflatedA
        std::vector<std::set<SizeType> > deflatedANZ(reduced_size);

        // Loop over non-zero structure of A and build non-zero structure of deflatedA
        typename SparseMatrixType::iterator1 a_iterator = rA.begin1();
        for(IndexType i = 0; i < full_size; ++i)
        {
            #ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            for(typename SparseMatrixType::iterator2 row_iterator = a_iterator.begin() ;
                                        row_iterator != a_iterator.end() ; ++row_iterator)
            {
            #else
            for (typename SparseMatrixType::iterator2 row_iterator = begin(a_iterator,
                    boost::numeric::ublas::iterator1_tag());
                    row_iterator != end(a_iterator, boost::numeric::ublas::iterator1_tag()); ++row_iterator )
            {
            #endif
                deflatedANZ[w[a_iterator.index1()]].insert(w[row_iterator.index2()]);
            }
            ++a_iterator;
        }

        // Count the number of non-zeros in deflatedA
        SizeType NZ = 0;
        for (IndexType i = 0; i < reduced_size; ++i)
            NZ += deflatedANZ[i].size();

        // Resize as needed the deflatedA
        deflatedA.resize(reduced_size, reduced_size, NZ);

        // Insert the non-zero structure into deflatedA
        for(IndexType i = 0 ; i < reduced_size ; ++i)
        {
            for(auto j = deflatedANZ[i].begin() ; j != deflatedANZ[i].end() ; ++j)
            {
                deflatedA.push_back(i, *j, 0.00);
            }
        }

        std::cout << "Non-zero structure of deflatedA is built, number of nonzeros = " << NZ << std::endl;
    }

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
    std::vector<int> mw;

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

}; // Class DeflatedSubdomainNodalBasedCGSolver

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

} // namespace Kratos.

#endif // MULTITHREADED_SOLVERS_APP_DEFLATED_SUBDOMAIN_NODAL_BASED_CG_SOLVER_H_INCLUDED  defined
