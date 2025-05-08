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
//   Date:                $Date: 28 Jun 2015 $
//   Revision:            $Revision: 1.3 $
//
//

#if !defined(KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BLOCK_JACOBI_PRESSURE_PRECONDITIONER_H_INCLUDED)
#define KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BLOCK_JACOBI_PRESSURE_PRECONDITIONER_H_INCLUDED


// System includes



// External includes
#include <boost/smart_ptr.hpp>
#include <boost/numeric/ublas/vector.hpp>


// Project includes
#include "includes/define.h"
#include "linear_solvers/linear_solver.h"
#include "utilities/openmp_utils.h"

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

/// BlockJacobiPressurePreconditioner class.
/**
This preconditioner treats the displacement and pressure differently in each block. It has the form
P=[A 0
   0 C]
where A and C will be preconditioned by normal preconditioner
*/
template<class TSparseSpaceType, class TDenseSpaceType, class TModelPartType>
class BlockJacobiPressurePreconditioner : public Preconditioner<TSparseSpaceType, TDenseSpaceType, TModelPartType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BlockJacobiPressurePreconditioner
    KRATOS_CLASS_POINTER_DEFINITION (BlockJacobiPressurePreconditioner);

    typedef Preconditioner<TSparseSpaceType, TDenseSpaceType, TModelPartType> BaseType;

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

    /// Default constructor.
    BlockJacobiPressurePreconditioner(typename BaseType::Pointer prec_A,
                                      typename BaseType::Pointer prec_C)
    {
        mprec_A = prec_A;
        mprec_C = prec_C;
    }

    /// Copy constructor.
    BlockJacobiPressurePreconditioner(const BlockJacobiPressurePreconditioner& Other)
    {
        mprec_A = Other.mprec_A;
        mprec_C = Other.mprec_C;
    }

    /// Destructor.
    ~BlockJacobiPressurePreconditioner() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    BlockJacobiPressurePreconditioner& operator=(const BlockJacobiPressurePreconditioner& Other)
    {
        mprec_A = Other.mprec_A;
        mprec_C = Other.mprec_C;
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    void Initialize(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        std::cout << "Fill blocks begin" << std::endl;
        double start = OpenMPUtils::GetCurrentTime();
        FillBlockMatrices(rA, mA, mB1, mB2, mC, mu, mp);
        std::cout << "Fill blocks completed..." << OpenMPUtils::GetCurrentTime() - start << " s" << std::endl;

        //this is rather slow
//        KRATOS_WATCH(norm_frobenius(mA))
//        KRATOS_WATCH(norm_frobenius(mB1))
//        KRATOS_WATCH(norm_frobenius(mB2))
//        KRATOS_WATCH(norm_frobenius(mC))

//        KRATOS_WATCH(ComputeFrobeniusNorm(mA))
//        KRATOS_WATCH(ComputeFrobeniusNorm(mB1))
//        KRATOS_WATCH(ComputeFrobeniusNorm(mB2))
//        KRATOS_WATCH(ComputeFrobeniusNorm(mC))

//        mprec_A->Initialize(mA, mu, mru);
        mprec_A->Initialize(mA, mu, rB); //take rB as temporary, but it should not be

//        mprec_C->Initialize(mS, mp, mrp);
        mprec_C->Initialize(mC, mp, rB); //take rB as temporary, but it should not be

        std::cout << "Block jacobi preconditioner is initialized" << std::endl;

        //debugging: try to rebuild the preconditioner
//        std::cout.precision(15);
//        int n = rB.size();
//        VectorType x(n);
//        typename TDenseSpaceType::MatrixType P(n, n);
//        for(int i = 0; i < n; ++i)
//        {
//            noalias(x) = ZeroVector(n);
//            x(i) = 1.0;
//            ApplyLeft(x);
//            noalias(column(P, i)) = x;
//        }
//        typename TDenseSpaceType::MatrixType rG = prod(P, rA);
//        KRATOS_WATCH(P)
//        KRATOS_WATCH(rA)
//        KRATOS_WATCH(rG)
//        exit(0);
    }

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
        //"other_counter[i]" i will contain the position in the global system of the i-th NON-pressure node
        //"pressure_counter[i]" will contain the in the global system of the i-th NON-pressure node
        //
        //mglobal_to_local_indexing[i] will contain the position in the local blocks of the
        unsigned int system_size = rA.size1();
        unsigned int pressure_counter = 0;
        unsigned int other_counter = 0;
        unsigned int global_pos;
        for(auto it = rdof_set.begin(); it != rdof_set.end(); ++it)
        {
            global_pos = it->EquationId();
            if(global_pos < system_size)
            {
                if( (it->GetVariable().Key() == WATER_PRESSURE) || (it->GetVariable().Key() == AIR_PRESSURE) )
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
                    ++other_counter;
                }
            }
        }

//        if(mprec_A->AdditionalPhysicalDataIsNeeded())
//            mprec_A->ProvideAdditionalData(rA, rX, rB, rdof_set, r_model_part); //TODO: the parameters need to be customized to the block A

//        if(mprec_C->AdditionalPhysicalDataIsNeeded())
//            mprec_C->ProvideAdditionalData(rA, rX, rB, rdof_set, r_model_part); //TODO: the parameters need to be customized to the block S
    }



    /** calculate preconditioned_u = A^{-1} * mu; preconditioned_p = C^{-1} * mp (for left-preconditioning)
        @param rX  Unknows of preconditioner suystem
    */
    VectorType& ApplyLeft(VectorType& rX) override
    {
        GetUPart(rX, mu);
        mprec_A->ApplyLeft(mu);
        WriteUPart(rX, mu);

        GetPPart(rX, mp);
        mprec_C->ApplyLeft(mp);
        WritePPart(rX, mp);

        return rX;
    }

    /** calculate preconditioned_u = A^{-1} * mu; preconditioned_p = C^{-1} * mp (for right-preconditioning)
        @param rX  Unknows of preconditioner suystem
    */
    VectorType& ApplyRight(VectorType& rX) override
    {
        GetUPart(rX, mu);
        mprec_A->ApplyRight(mu);
        WriteUPart(rX, mu);

        GetPPart(rX, mp);
        mprec_C->ApplyRight(mp);
        WritePPart(rX, mp);

        return rX;
    }

    VectorType& ApplyInverseRight(VectorType& rX) override
    {
        GetUPart(rX, mu);
        mprec_A->ApplyInverseRight(mu);
        WriteUPart(rX, mu);

        GetPPart(rX, mp);
        mprec_C->ApplyInverseRight(mp);
        WritePPart(rX, mp);

        return rX;
    }

    VectorType& ApplyTransposeLeft(VectorType& rX) override
    {
        GetUPart(rX, mu);
        mprec_A->ApplyTransposeLeft(mu);
        WriteUPart(rX, mu);

        GetPPart(rX, mp);
        mprec_C->ApplyTransposeLeft(mp);
        WritePPart(rX, mp);

        return rX;
    }

    VectorType& ApplyTransposeRight(VectorType& rX) override
    {
        GetUPart(rX, mu);
        mprec_A->ApplyTransposeRight(mu);
        WriteUPart(rX, mu);

        GetPPart(rX, mp);
        mprec_C->ApplyTransposeRight(mp);
        WritePPart(rX, mp);

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

    /// Return information about this object.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "BlockJacobiPressurePreconditioner";
        buffer << ", P_A = {" << mprec_A->Info() << "}";
        buffer << ", P_C = {" << mprec_C->Info() << "}";
        return buffer.str();
    }

    /// Print information about this object.
    void  PrintInfo(std::ostream& OStream) const override
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

    std::map<unsigned int, unsigned int> mpressure_indices;
    std::map<unsigned int, unsigned int> mother_indices;
    std::map<unsigned int, unsigned int> mglobal_to_local_indexing;
    std::map<unsigned int, bool> mis_pressure_block;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    //this function extracts from a vector which has the size of the
    //overall r, the part that corresponds to u-dofs
    void GetUPart (const VectorType& rtot, VectorType& ru) const
    {
        if (ru.size() != mother_indices.size() )
            ru.resize (mother_indices.size(), false);
        #pragma omp parallel for
        for (unsigned int i = 0; i < static_cast<int>(ru.size()); ++i)
        {
            auto it = mother_indices.find(i);
            if (it == mother_indices.end())
                KRATOS_ERROR << "Index " << i << " does not exist";
            ru[i] = rtot[it->second];
        }
    }

    //this function extracts from a vector which has the size of the
    //overall r, the part that corresponds to p-dofs
    void GetPPart (const VectorType& rtot, VectorType& rp) const
    {
        if (rp.size() != mpressure_indices.size() )
            rp.resize (mpressure_indices.size(), false);
        #pragma omp parallel for
        for (unsigned int i = 0; i < static_cast<int>(rp.size()); ++i)
        {
            auto it = mpressure_indices.find(i);
            if (it == mpressure_indices.end())
                KRATOS_ERROR << "Index " << i << " does not exist";
            rp[i] = rtot[it->second];
        }
    }

    void WriteUPart (VectorType& rtot, const VectorType& ru) const
    {
        #pragma omp parallel for
        for (unsigned int i = 0; i < static_cast<int>(ru.size()); ++i)
        {
            auto it = mother_indices.find(i);
            if (it == mother_indices.end())
                KRATOS_ERROR << "Index " << i << " does not exist";
            rtot[it->second] = ru[i];
        }
    }

    void WritePPart (VectorType& rtot, const VectorType& rp) const
    {
        #pragma omp parallel for
        for (unsigned int i = 0; i < static_cast<int>(rp.size()); ++i)
        {
            auto it = mpressure_indices.find(i);
            if (it == mpressure_indices.end())
                KRATOS_ERROR << "Index " << i << " does not exist";
            rtot[it->second] = rp[i];
        }
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

    typename BaseType::Pointer mprec_A;
    typename BaseType::Pointer mprec_C;

    SparseMatrixType mA;
    SparseMatrixType mB1;
    SparseMatrixType mB2;
    SparseMatrixType mC;

    VectorType mp;
    VectorType mu;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    ///this function generates the subblocks of matrix J
    ///as J = ( A  B1 ) u
    ///       ( B2 C  ) p
    void FillBlockMatrices (const SparseMatrixType& rA,
        SparseMatrixType& A, SparseMatrixType& B1, SparseMatrixType& B2, SparseMatrixType& C,
        VectorType& u, VectorType& p) const
    {
        KRATOS_TRY

        //get access to J data
        const SizeType* index1 = rA.index1_data().begin();
        const SizeType* index2 = rA.index2_data().begin();
        const double*   values = rA.value_data().begin();

        A.clear();
        B1.clear();
        B2.clear();
        C.clear();

        //do allocation
        TSparseSpaceType::Resize(A,  mother_indices.size(), mother_indices.size());
        TSparseSpaceType::Resize(B1, mother_indices.size(), mpressure_indices.size());
        TSparseSpaceType::Resize(B2, mpressure_indices.size(), mother_indices.size());
        TSparseSpaceType::Resize(C,  mpressure_indices.size(), mpressure_indices.size());

        TSparseSpaceType::Resize(p, mpressure_indices.size());
        TSparseSpaceType::Set(p, 0.00);
        TSparseSpaceType::Resize(u, mother_indices.size());
        TSparseSpaceType::Set(u, 0.00);

        //allocate the blocks by push_back
        for (unsigned int i = 0; i < rA.size1(); ++i)
        {
            unsigned int row_begin = index1[i];
            unsigned int row_end   = index1[i + 1];
            unsigned int local_row_id = mglobal_to_local_indexing.at(i);

            if ( mis_pressure_block.at(i) == false) //either A or B1
            {
                for (unsigned int j = row_begin; j < row_end; ++j)
                {
                    unsigned int col_index = index2[j];
                    double value = values[j];
                    unsigned int local_col_id = mglobal_to_local_indexing.at(col_index);
                    if (mis_pressure_block.at(col_index) == false) //A block
                        A.push_back ( local_row_id, local_col_id, value);
                    else //B1 block
                        B1.push_back ( local_row_id, local_col_id, value);
                }
            }
            else //either B2 or C
            {
                for (unsigned int j = row_begin; j < row_end; ++j)
                {
                    unsigned int col_index = index2[j];
                    double value = values[j];
                    unsigned int local_col_id = mglobal_to_local_indexing.at(col_index);
                    if (mis_pressure_block.at(col_index) == false) //B2 block
                        B2.push_back ( local_row_id, local_col_id, value);
                    else //C block
                        C.push_back ( local_row_id, local_col_id, value);
                }
            }
        }

        KRATOS_CATCH ("")
    }

    double ComputeFrobeniusNorm(const SparseMatrixType& rA)const
    {
        int n = rA.size1();
        const std::size_t* ia = rA.index1_data().begin();
        const std::size_t* ja = rA.index2_data().begin();
        const double*      a  = rA.value_data().begin();

        double norm = 0.0;
//        KRATOS_WATCH(rA.size1())
//        KRATOS_WATCH(rA.size2())
        for(int i = 0; i < n; ++i)
        {
//            std::cout << "i: " << i << ", ia[i]: " << ia[i] << std::endl;
            int nz = ia[i + 1] - ia[i];
            for(int j = 0; j < nz; ++j)
                norm += pow(a[ia[i] + j], 2);
        }
        return std::sqrt(norm);
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

}; // Class BlockJacobiPressurePreconditioner

///@}


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

}  // namespace Kratos.

#endif // KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BLOCK_JACOBI_PRESSURE_PRECONDITIONER_H_INCLUDED  defined
