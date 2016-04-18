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

#if !defined(KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BLOCK_JACOBI_PRECONDITIONER_H_INCLUDED)
#define KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BLOCK_JACOBI_PRECONDITIONER_H_INCLUDED


// System includes



// External includes
#include <boost/smart_ptr.hpp>
#include <boost/numeric/ublas/vector.hpp>


// Project includes
#include "includes/define.h"
#include "linear_solvers/linear_solver.h"
#include "utilities/openmp_utils.h"


#define PROBE_BLOCK_PROPERTIES


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

/// BlockJacobiPreconditioner class.
/**   */
template<class TSparseSpaceType, class TDenseSpaceType>
class BlockJacobiPreconditioner : public Preconditioner<TSparseSpaceType, TDenseSpaceType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BlockJacobiPreconditioner
    KRATOS_CLASS_POINTER_DEFINITION (BlockJacobiPreconditioner);

    typedef Preconditioner<TSparseSpaceType, TDenseSpaceType> BaseType;

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
    BlockJacobiPreconditioner()
    {
    }

    /// Copy constructor.
    BlockJacobiPreconditioner(const BlockJacobiPreconditioner& rOther)
    {
//        rOther.mBlockIndices = mBlockIndices;
//        rOther.mpPrecs = mpPrecs;
//        rOther.mBlockMats = mBlockMats;
//        rOther.mBlockVecs = mBlockVecs;
    }

    /// Destructor.
    virtual ~BlockJacobiPreconditioner()
    {
    }


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    BlockJacobiPreconditioner& operator=(const BlockJacobiPreconditioner& rOther)
    {
//        rOther.mBlockIndices = mBlockIndices;
//        rOther.mpPrecs = mpPrecs;
//        rOther.mBlockMats = mBlockMats;
//        rOther.mBlockVecs = mBlockVecs;
        return *this;
    }


    ///@}
    ///@name Operations
    ///@{

    void AddPreconditioner(typename BaseType::Pointer pPrec)
    {
        mpPrecs.push_back(pPrec);
    }

    void SetPreconditioner(unsigned int i, typename BaseType::Pointer pPrec)
    {
        mpPrecs[i] = pPrec;
    }

    unsigned int NumberOfBlocks() const
    {
        return mBlockIndices.size();
    }

    typename BaseType::Pointer GetPreconditioner(unsigned int i) const
    {
        return mpPrecs[i];
    }

    virtual void Initialize(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
        /* check if the preconditioner has been set */
        unsigned int num_blocks = mBlockIndices.size();
        if(mpPrecs.size() != num_blocks)
            KRATOS_THROW_ERROR(std::logic_error, "The number of sub-preconditioners is not equal to number of blocks", "")

        /* initialize the preconditioner */
        for(unsigned int i = 0; i < mpPrecs.size(); ++i)
        {
            mpPrecs[i]->Initialize(mBlockMats[i], mBlockVecs[i], mBlockVecs[i]); // should we do initialize for rX as well?
        }
    }

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
        /* check if the preconditioner has been set */
        unsigned int num_blocks = mBlockIndices.size();
        if(mpPrecs.size() != num_blocks)
            KRATOS_THROW_ERROR(std::logic_error, "The number of sub-preconditioners is not equal to number of blocks", "")

        /* fill the block matrices and vectors */
        std::cout << "Fill blocks begin" << std::endl;
        double start = OpenMPUtils::GetCurrentTime();
        FillBlockMatrices(rA, mBlockIndices, mBlockMats);
        FillBlockVectors(rB, mBlockIndices, mBlockVecs);
        std::cout << "Fill block matrices and vectors completed..." << OpenMPUtils::GetCurrentTime() - start << " s" << std::endl;

        #ifdef PROBE_BLOCK_PROPERTIES
        for(unsigned int i = 0; i < num_blocks; ++i)
        {
            std::cout << "Frobenius norm of block " << i << ": " << ComputeFrobeniusNorm(mBlockMats[i]) << std::endl;
            std::cout << "Largest (abs) diagonal entry of block " << i << ": " << ExtractLargestDiagonalEntry(mBlockMats[i]) << std::endl;
            std::cout << "Smallest (abs) diagonal entry of block " << i << ": " << ExtractSmallestDiagonalEntry(mBlockMats[i]) << std::endl;
            std::cout << "Average (abs) of diagonal entries of block " << i << ": " << ComputeAverageAbsDiagonal(mBlockMats[i]) << std::endl;
        }
        #endif

        /* provide additional data for each preconditioner */
        for(unsigned int i = 0; i < mpPrecs.size(); ++i)
            mpPrecs[i]->ProvideAdditionalData(mBlockMats[i], mBlockVecs[i], mBlockVecs[i], mBlockDofs[i], r_model_part);
    }


    /** calculate preconditioned_u = A^{-1} * mu; preconditioned_p = C^{-1} * mp (for left-preconditioning)
        @param rX  Unknows of preconditioner suystem
    */
    virtual VectorType& ApplyLeft(VectorType& rX)
    {
        VectorType U;
        for(unsigned int i = 0; i < mBlockIndices.size(); ++i)
        {
            U.resize(mBlockIndices[i].size());
            ExtractVector(mBlockIndices[i], rX, U);
            mpPrecs[i]->ApplyLeft(U);
            AssignVector(mBlockIndices[i], rX, U);
        }
        return rX;
    }

    /** calculate preconditioned_u = A^{-1} * mu; preconditioned_p = C^{-1} * mp (for right-preconditioning)
        @param rX  Unknows of preconditioner suystem
    */
    virtual VectorType& ApplyRight(VectorType& rX)
    {
        VectorType U;
        for(unsigned int i = 0; i < mBlockIndices.size(); ++i)
        {
            U.resize(mBlockIndices[i].size());
            ExtractVector(mBlockIndices[i], rX, U);
            mpPrecs[i]->ApplyRight(U);
            AssignVector(mBlockIndices[i], rX, U);
        }
        return rX;
    }

    virtual VectorType& ApplyInverseRight(VectorType& rX)
    {
        VectorType U;
        for(unsigned int i = 0; i < mBlockIndices.size(); ++i)
        {
            U.resize(mBlockIndices[i].size());
            ExtractVector(mBlockIndices[i], rX, U);
            mpPrecs[i]->ApplyInverseRight(U);
            AssignVector(mBlockIndices[i], rX, U);
        }
        return rX;
    }

    virtual VectorType& ApplyTransposeLeft(VectorType& rX)
    {
        VectorType U;
        for(unsigned int i = 0; i < mBlockIndices.size(); ++i)
        {
            U.resize(mBlockIndices[i].size());
            ExtractVector(mBlockIndices[i], rX, U);
            mpPrecs[i]->ApplyTransposeLeft(U);
            AssignVector(mBlockIndices[i], rX, U);
        }
        return rX;
    }

    virtual VectorType& ApplyTransposeRight(VectorType& rX)
    {
        VectorType U;
        for(unsigned int i = 0; i < mBlockIndices.size(); ++i)
        {
            U.resize(mBlockIndices[i].size());
            ExtractVector(mBlockIndices[i], rX, U);
            mpPrecs[i]->ApplyTransposeRight(U);
            AssignVector(mBlockIndices[i], rX, U);
        }
        return rX;
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

    virtual std::string Name() const {return "BlockJacobiPreconditioner";}

    /// Return information about this object.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << this->Name() << ", list of sub-preconditioners:";
        for(unsigned int i = 0; i < mpPrecs.size(); ++i)
            buffer << " {" << mpPrecs[i]->Info() << "}";
        return buffer.str();
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

    std::vector<std::set<SizeType> > mBlockIndices;

    std::vector<ModelPart::DofsArrayType> mBlockDofs;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    //extract sub-vector rU from the total vector rX
    void ExtractVector(const std::set<SizeType>& rIndices, const VectorType& rX, VectorType& rU) const
    {
        unsigned int i = 0;
        for(std::set<SizeType>::iterator it = rIndices.begin(); it != rIndices.end(); ++it)
            rU[i++] = rX[*it];
    }

    //assign the values of sub-vector to part of the total vector
    void AssignVector(const std::set<SizeType>& rIndices, VectorType& rX, const VectorType& rU) const
    {
        unsigned int i = 0;
        for(std::set<SizeType>::iterator it = rIndices.begin(); it != rIndices.end(); ++it)
            rX[*it] = rU[i++];
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

    std::vector<typename BaseType::Pointer> mpPrecs;

    std::vector<SparseMatrixType> mBlockMats;

    std::vector<VectorType> mBlockVecs;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /// fill the block matrices from the block indices
    void FillBlockMatrices(const SparseMatrixType& rA,
                           const std::vector<std::set<SizeType> >& rBlockIndices,
                           std::vector<SparseMatrixType>& rBlockMats) const
    {
        KRATOS_TRY

        //get access to rA data
        unsigned int system_size = rA.size1();
        const SizeType* index1 = rA.index1_data().begin();
        const SizeType* index2 = rA.index2_data().begin();
        const double*   values = rA.value_data().begin();

        rBlockMats.clear();
        rBlockMats.resize(rBlockIndices.size());
        for(unsigned int i = 0; i < rBlockIndices.size(); ++i)
        {
            // clear and allocate the block matrix
            rBlockMats[i].clear();
            TSparseSpaceType::Resize(rBlockMats[i], rBlockIndices[i].size(), rBlockIndices[i].size());

            // extract the block indices and create the map from global indexing to locel indexing
            std::vector<SizeType> block_indices(rBlockIndices[i].begin(), rBlockIndices[i].end()); // make a copy of the block indices
            std::map<unsigned int, unsigned int> global_to_local_map;
            std::map<unsigned int, bool> is_block_local;

            for(unsigned int j = 0; j < system_size; ++j)
                is_block_local[j] = false;

            for(unsigned int j = 0; j < block_indices.size(); ++j)
            {
                global_to_local_map[block_indices[j]] = j;
                is_block_local[block_indices[j]] = true;
            }
            KRATOS_WATCH(block_indices.size())

            // fill the matrix by push_back the values by row
            unsigned int local_row_id = 0, local_col_id;
            for(unsigned int j = 0; j < block_indices.size(); ++j)
            {
                unsigned int row = block_indices[j];
                unsigned int row_idx_begin = index1[row];
                unsigned int row_idx_end   = index1[row + 1];
                for(unsigned int k = row_idx_begin; k < row_idx_end; ++k)
                {
                    unsigned int col = index2[k];
                    if(is_block_local[col])
                    {
                        local_col_id = global_to_local_map[col];
                        rBlockMats[i].push_back(local_row_id, local_col_id, values[k]);
                    }
                }
                ++local_row_id;
            }

            // complete the fill by calling complete_index1_data
            rBlockMats[i].complete_index1_data();
        }

        KRATOS_CATCH ("")
    }

    /// fill the block vectors from the block indices
    void FillBlockVectors(const VectorType& rV,
                           const std::vector<std::set<SizeType> >& rBlockIndices,
                           std::vector<VectorType>& rBlockVecs) const
    {
        rBlockVecs.clear();
        rBlockVecs.resize(rBlockIndices.size());
        for(unsigned int i = 0; i < rBlockIndices.size(); ++i)
        {
            // clear and allocate the vector
            rBlockVecs[i].clear();
            rBlockVecs[i].resize(rBlockIndices[i].size());
            // fill the vector
            unsigned int j = 0;
            for(std::set<SizeType>::const_iterator it = rBlockIndices[i].begin(); it != rBlockIndices[i].end(); ++it)
                rBlockVecs[i][j++] = rV[*it];
        }
    }

    /// fast function to compute the Frobenius norm of the sparse matrix
    double ComputeFrobeniusNorm(SparseMatrixType& rA)
    {
        int n = rA.size1();
        const std::size_t* ia = rA.index1_data().begin();
        const std::size_t* ja = rA.index2_data().begin();
        const double*	   a  = rA.value_data().begin();

        double norm = 0.0;
        for(int i = 0; i < n; ++i)
        {
            int nz = ia[i + 1] - ia[i];
            for(int j = 0; j < nz; ++j)
                norm += pow(a[ia[i] + j], 2);
        }

        return sqrt(norm);
    }

    /// compute the largest diagonal entry of matrix A
    double ExtractLargestDiagonalEntry(SparseMatrixType& rA)
    {
        int n = rA.size1();
        double max = -static_cast<double>(INT_MAX);
        for(unsigned int i = 0; i < n; ++i)
            if(fabs(rA(i, i)) > max)
                max = fabs(rA(i, i));
        return max;
    }

    /// compute the smallest diagonal entry of matrix A
    double ExtractSmallestDiagonalEntry(SparseMatrixType& rA)
    {
        int n = rA.size1();
        double min = static_cast<double>(INT_MAX);
        for(unsigned int i = 0; i < n; ++i)
            if(fabs(rA(i, i)) < min)
                min = fabs(rA(i, i));
        return min;
    }

    /// compute the average value (abs) of the diagonal
    double ComputeAverageAbsDiagonal(SparseMatrixType& rA)
    {
        int n = rA.size1();
        double sum = 0.0;
        for(unsigned int i = 0; i < n; ++i)
            sum += fabs(rA(i, i));
        return sum / n;
    }
    
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

}; // Class BlockJacobiPreconditioner

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
                                  BlockJacobiPreconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
{
    return IStream;
}


/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType>
inline std::ostream& operator << (std::ostream& OStream,
                                  const BlockJacobiPreconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
{
    rThis.PrintInfo(OStream);
    OStream << std::endl;
    rThis.PrintData(OStream);


    return OStream;
}
///@}


}  // namespace Kratos.

#undef PROBE_BLOCK_PROPERTIES

#endif // KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BLOCK_JACOBI_PRESSURE_PRECONDITIONER_H_INCLUDED  defined 

