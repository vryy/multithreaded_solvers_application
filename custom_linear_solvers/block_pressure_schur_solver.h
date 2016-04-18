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


#if !defined(KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BLOCK_PRESSURE_SCHUR_SOLVER_H_INCLUDED )
#define KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BLOCK_PRESSURE_SHUR_SOLVER_H_INCLUDED


// System includes
#include <cmath>


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "linear_solvers/linear_solver.h"


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
class BlockPressureSchurSolver : public LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of  BlockPressureSchurSolver
    KRATOS_CLASS_POINTER_DEFINITION( BlockPressureSchurSolver );

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
    
    typedef std::size_t  SizeType;
    
    typedef std::size_t  IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BlockPressureSchurSolver() {}

    BlockPressureSchurSolver(
        typename BaseType::Pointer pSolver
    ) : BaseType()
    {
        mpSolver = pSolver;
    }

    /// Copy constructor.
     BlockPressureSchurSolver(const  BlockPressureSchurSolver& Other) : BaseType(Other) {}

    /// Destructor.
    virtual ~ BlockPressureSchurSolver() {}


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
                    ++other_counter;
                }
            }
        }
        
        if(mpSolver->AdditionalPhysicalDataIsNeeded())
            mpSolver->ProvideAdditionalData(rA, rX, rB, rdof_set, r_model_part);
    }
    
    virtual void Initialize(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
        std::cout << "Fill blocks begin" << std::endl;
        double start = OpenMPUtils::GetCurrentTime();
        FillBlockMatrices(rA, mA, mB1, mB2, mC);
        std::cout << "Fill blocks completed..." << OpenMPUtils::GetCurrentTime() - start << " s" << std::endl;
        
        //this is rather slow
//        KRATOS_WATCH(norm_frobenius(mA))
//        KRATOS_WATCH(norm_frobenius(mB1))
//        KRATOS_WATCH(norm_frobenius(mB2))
//        KRATOS_WATCH(norm_frobenius(mC))
        
        KRATOS_WATCH(ComputeFrobeniusNorm(mA))
        KRATOS_WATCH(ComputeFrobeniusNorm(mB1))
        KRATOS_WATCH(ComputeFrobeniusNorm(mB2))
        KRATOS_WATCH(ComputeFrobeniusNorm(mC))
        
        // TODO: check the diagonal dominance of matrix A
        #ifdef CHECK_DIAGONAL_DOMINANCE
        CheckDiagonalDominance(mA);
        #endif

        mpSolver->Initialize(mA, mu, rB); //take rB as temporary, but it should not be
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
        
        Initialize(rA, rX, rB);
        
        VectorType ru, rp, u, p;

        // Get the initial u & p
        GetUPart(rX, u);
        GetPPart(rX, p);
        
        // Get ru, rp
        GetUPart(rB, ru);
        GetPPart(rB, rp);
        
        // Extract inverse of diagonal of mA
        int n = mA.size1();
        VectorType invDiagBlockA(n);
        for(unsigned int i = 0; i < n; ++i)
            invDiagBlockA[i] = 1.0 / mA(i, i);
//        KRATOS_WATCH(invDiagBlockA)
//        KRATOS_WATCH(ru)
//        KRATOS_WATCH(rp)
        
        // Solve for p
        VectorType rhs_p = ru;
        VectorScale(rhs_p, invDiagBlockA);
        rhs_p = rp - prod(mB2, rhs_p);
//        KRATOS_WATCH(rhs_p)
        
        SparseMatrixType lhs_p = mB1;
        RowScale(lhs_p, invDiagBlockA);
        lhs_p = mC - SparseMatrixType(prod(mB2, lhs_p));
        
        mpSolver->Solve(lhs_p, p, rhs_p);
//        KRATOS_WATCH(p)
        
        // Solve for u
        noalias(ru) = ru - prod(mB1, p);
        mpSolver->Solve(mA, u, ru);
//        KRATOS_WATCH(u)
        
        // Write back the result to solution vector
        WriteUPart(rX, u);
        WritePPart(rX, p);
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
        buffer << "Linear solver using Schur complement reduction scheme, linear solver = " << mpSolver->Info();
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

    std::vector<SizeType> mpressure_indices;
    std::vector<SizeType> mother_indices;
    std::vector<int> mglobal_to_local_indexing;
    std::vector<int> mis_pressure_block;
    
    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    //this function extracts from a vector which has the size of the
    //overall r, the part that corresponds to u-dofs
    void GetUPart (const VectorType& rtot, VectorType& ru)
    {
        if (ru.size() != mother_indices.size() )
            ru.resize (mother_indices.size(), false);
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(ru.size()); ++i)
            ru[i] = rtot[mother_indices[i]];
    }

    //this function extracts from a vector which has the size of the
    //overall r, the part that corresponds to p-dofs
    void GetPPart (const VectorType& rtot, VectorType& rp)
    {
        if (rp.size() != mpressure_indices.size() )
            rp.resize (mpressure_indices.size(), false);
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(rp.size()); ++i)
            rp[i] = rtot[mpressure_indices[i]];
    }

    void WriteUPart (VectorType& rtot, const VectorType& ru)
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(ru.size()); ++i)
            rtot[mother_indices[i]] = ru[i];
    }

    void WritePPart (VectorType& rtot, const VectorType& rp)
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(rp.size()); ++i)
            rtot[mpressure_indices[i]] = rp[i];
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

    typename BaseType::Pointer mpSolver;
    
    SparseMatrixType mA;
    SparseMatrixType mB1;
    SparseMatrixType mB2;
    SparseMatrixType mC;
    SparseMatrixType mS;
    
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
    void FillBlockMatrices (SparseMatrixType& rA, SparseMatrixType& A, SparseMatrixType& B1, SparseMatrixType& B2, SparseMatrixType& C)
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
	    
	    TSparseSpaceType::Resize(mp, mpressure_indices.size());
	    TSparseSpaceType::Set(mp, 0.00);
	    TSparseSpaceType::Resize(mu, mother_indices.size());
	    TSparseSpaceType::Set(mu, 0.00);

        //allocate the blocks by push_back
        for (unsigned int i = 0; i < rA.size1(); ++i)
        {
            unsigned int row_begin = index1[i];
            unsigned int row_end   = index1[i + 1];
            unsigned int local_row_id = mglobal_to_local_indexing[i];

            if ( mis_pressure_block[i] == false) //either A or B1
            {
                for (unsigned int j = row_begin; j < row_end; ++j)
                {
                    unsigned int col_index = index2[j];
                    double value = values[j];
                    unsigned int local_col_id = mglobal_to_local_indexing[col_index];
                    if (mis_pressure_block[col_index] == false) //A block
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
                    unsigned int local_col_id = mglobal_to_local_indexing[col_index];
                    if (mis_pressure_block[col_index] == false) //B2 block
                        B2.push_back ( local_row_id, local_col_id, value);
                    else //C block
                        C.push_back ( local_row_id, local_col_id, value);
                }
            }
        }

        KRATOS_CATCH ("")
    }
    
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
    
    void VectorScale(VectorType& rX, VectorType& rD)
    {
        int n = rX.size();

        #pragma omp parallel for
        for(int i = 0; i < n; ++i)
            rX(i) *= rD(i);
    }
    
    void RowScale(SparseMatrixType& rA, VectorType& rDL)
    {
        int n = rA.size1();
        std::size_t*   ia = rA.index1_data().begin();
        std::size_t*   ja = rA.index2_data().begin();
        double*         a = rA.value_data().begin();
        
        #pragma omp parallel for
        for(int i = 0; i < n; ++i)
        {
            int nz = ia[i + 1] - ia[i];
            for(int j = 0; j < nz; ++j)
                a[ia[i] + j] *= rDL(i);
        }
    }
    
    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{

    void CheckDiagonalDominance(SparseMatrixType& rA)
    {
        int n = rA.size1();
        const std::size_t* ia = rA.index1_data().begin();
        const std::size_t* ja = rA.index2_data().begin();
        const double*	   a  = rA.value_data().begin();
        
        double diag_norm, off_diag_norm;
        int nz;
        for(int i = 0; i < n; ++i)
        {
            nz = ia[i + 1] - ia[i];
            off_diag_norm = 0.0;
            for(int j = 0; j < nz; ++j)
            {
                if(ja[ia[i] + j] == i)
                    diag_norm = fabs(a[ia[i] + j]);
                else
                    off_diag_norm += fabs(a[ia[i] + j]);
            }
            if(diag_norm < off_diag_norm)
                std::cout << "Matrix A is not diagonal dominant at row " << i << ", the ratio is " << diag_norm / off_diag_norm << std::endl; 
        }
    }

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


/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::istream& operator >> (std::istream& IStream, BlockPressureSchurSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis)
{
    return IStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator << (std::ostream& OStream, const  BlockPressureSchurSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis)
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

#endif //  KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BLOCK_PRESSURE_SHUR_SOLVER_H_INCLUDED  defined 

