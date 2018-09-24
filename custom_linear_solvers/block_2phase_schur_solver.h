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


#if !defined(KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BLOCK_2_PHASE_SCHUR_SOLVER_H_INCLUDED )
#define KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BLOCK_2_PHASE_SCHUR_SOLVER_H_INCLUDED


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
THis solver solve the coupled u-p problem by using one-step staggerred scheme:
[A  B1][u] = [ru]
[B2  C][p]   [rp]

=> u = A^(-1) * (ru -B1 * p)                                (1)
   (C - B2 * A^(-1) * B1) * p = (rp - B2 * A^(-1) * ru)     (2)
   
Firstly p in (2) is solved approximately by approximating A^(-1) by diag(A)^(-1) (remarks: A needs to be diagonal-dominant for this to be a good approximation)
Then u in (1) is solved exactly
*/
template<class TSparseSpaceType, class TDenseSpaceType,
         class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class Block2PhaseSchurSolver : public LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of  Block2PhaseSchurSolver
    KRATOS_CLASS_POINTER_DEFINITION( Block2PhaseSchurSolver );

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
    Block2PhaseSchurSolver() {}

    Block2PhaseSchurSolver(
        typename BaseType::Pointer pSolver
    ) : BaseType(), mpSolver(pSolver)
    {}

    /// Copy constructor.
    Block2PhaseSchurSolver(const  Block2PhaseSchurSolver& Other) : BaseType(Other) {}

    /// Destructor.
    virtual ~Block2PhaseSchurSolver() {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    Block2PhaseSchurSolver& operator=(const  Block2PhaseSchurSolver& Other)
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
        std::cout << "Initialize completed" << std::endl;

        KRATOS_WATCH(rA.size1())
        KRATOS_WATCH(rA.size2())

        VectorType ru, rp, u, p;

        // Get the initial u & p
        GetUPart(rX, u);
        GetPPart(rX, p);

        // Get ru, rp
        GetUPart(rB, ru);
        GetPPart(rB, rp);

        KRATOS_WATCH(u.size())
        KRATOS_WATCH(p.size())
        KRATOS_WATCH(ru.size())
        KRATOS_WATCH(rp.size())

        // Extract inverse of diagonal of mA
        std::size_t n = mA.size1();
        KRATOS_WATCH(n)
        VectorType invDiagBlockA(n);
        for(std::size_t i = 0; i < n; ++i)
            invDiagBlockA[i] = 1.0 / mA(i, i);
//        KRATOS_WATCH(invDiagBlockA)
//        KRATOS_WATCH(ru)
//        KRATOS_WATCH(rp)

        // Solve for p
        VectorType tmp_u = ru;
        KRATOS_WATCH(norm_2(tmp_u))
//        KRATOS_WATCH(tmp_u)
        KRATOS_WATCH(norm_2(invDiagBlockA))
        VectorScale(tmp_u, invDiagBlockA);
//        KRATOS_WATCH(tmp_u)
        KRATOS_WATCH(norm_2(tmp_u))
        Vector aux1(mB2.size1());
        TSparseSpaceType::Mult(mB2, tmp_u, aux1);
        KRATOS_WATCH(norm_2(aux1))
        VectorType rhs_p = rp - aux1;
        KRATOS_WATCH(norm_2(rhs_p))
//        KRATOS_WATCH(rhs_p.size())
//        KRATOS_WATCH(rp.size())
//        KRATOS_WATCH(mB2.size1())
//        KRATOS_WATCH(mB2.size2())
//        KRATOS_WATCH(rhs_p.size())
        std::cout << "Create RHS for p completed" << std::endl;

        SparseMatrixType tmpS = mB1;
        RowScale(tmpS, invDiagBlockA);
        DenseMatrixType S(mC.size1(), mC.size2());
        MatrixMult(S, mB2, tmpS);
        S *= -1.0;
        noalias(S) += mC;
        std::cout << "Compute Schur completed" << std::endl;
        KRATOS_WATCH(norm_frobenius(S))

        SparseMatrixType sS = S;
        KRATOS_WATCH(ComputeFrobeniusNorm(sS))
        mpSolver->Solve(sS, p, rhs_p);
//        KRATOS_WATCH(p)
        KRATOS_WATCH(norm_2(p))
        std::cout << "Solve for p completed" << std::endl;

        // Solve for u
        Vector aux2(ru.size());
        TSparseSpaceType::Mult(mB1, p, aux2);
        noalias(ru) -= aux2;
        KRATOS_WATCH(norm_2(ru))
        mpSolver->Solve(mA, u, ru);
//        KRATOS_WATCH(u)
        KRATOS_WATCH(norm_2(u))
        std::cout << "Solve for u completed" << std::endl;

        // Write back the result to solution vector
        WriteUPart(rX, u);
        WritePPart(rX, p);

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
        buffer << "Linear solver using Schur complement reduction scheme for 2 phases, linear solver = " << mpSolver->Info();
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
//        #pragma omp parallel for
        for (std::size_t i = 0; i < ru.size(); ++i)
            ru[i] = rtot[mother_indices[i]];
    }

    //this function extracts from a vector which has the size of the
    //overall r, the part that corresponds to p-dofs
    void GetPPart (const VectorType& rtot, VectorType& rp)
    {
        if (rp.size() != mpressure_indices.size() )
            rp.resize (mpressure_indices.size(), false);
//        #pragma omp parallel for
        for (std::size_t i = 0; i < rp.size(); ++i)
            rp[i] = rtot[mpressure_indices[i]];
    }

    void WriteUPart (VectorType& rtot, const VectorType& ru)
    {
//        #pragma omp parallel for
        for (std::size_t i = 0; i < ru.size(); ++i)
            rtot[mother_indices[i]] = ru[i];
    }

    void WritePPart (VectorType& rtot, const VectorType& rp)
    {
//        #pragma omp parallel for
        for (std::size_t i = 0; i < rp.size(); ++i)
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
//                        A.push_back ( local_row_id, local_col_id, value);
                        A( local_row_id, local_col_id ) = value;
                    else //B1 block
//                        B1.push_back ( local_row_id, local_col_id, value);
                        B1( local_row_id, local_col_id ) = value;
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
//                        B2.push_back ( local_row_id, local_col_id, value);
                        B2( local_row_id, local_col_id ) = value;
                    else //C block
//                        C.push_back ( local_row_id, local_col_id, value);
                        C( local_row_id, local_col_id ) = value;
                }
            }
        }

        A.complete_index1_data();
        B1.complete_index1_data();
        B2.complete_index1_data();
        C.complete_index1_data();

        KRATOS_CATCH ("")
    }

    double ComputeFrobeniusNorm(SparseMatrixType& rA)
    {
        std::size_t n = rA.size1();
        const std::size_t* ia = rA.index1_data().begin();
        const std::size_t* ja = rA.index2_data().begin();
        const double*	   a  = rA.value_data().begin();
        
        double norm = 0.0;
        for(std::size_t i = 0; i < n; ++i)
        {
            int nz = ia[i + 1] - ia[i];
            for(int j = 0; j < nz; ++j)
                norm += pow(a[ia[i] + j], 2);
        }
        return sqrt(norm);
    }

    void VectorScale(VectorType& rX, const VectorType& rD)
    {
        std::size_t n = rX.size();

//        #pragma omp parallel for
        for(std::size_t i = 0; i < n; ++i)
            rX(i) *= rD(i);
    }

    void RowScale(SparseMatrixType& rA, const VectorType& rDL)
    {
        std::size_t n = rA.size1();
        std::size_t*   ia = rA.index1_data().begin();
        std::size_t*   ja = rA.index2_data().begin();
        double*         a = rA.value_data().begin();
        
//        #pragma omp parallel for
        for(std::size_t i = 0; i < n; ++i)
        {
            int nz = ia[i + 1] - ia[i];
            for(int j = 0; j < nz; ++j)
                a[ia[i] + j] *= rDL(i);
        }
    }

    void MatrixMult(DenseMatrixType& rC, const SparseMatrixType& rA, const DenseMatrixType& rB)
    {
        std::size_t n = rA.size1();
        std::size_t m = rB.size2();
        const std::size_t*   ia = rA.index1_data().begin();
        const std::size_t*   ja = rA.index2_data().begin();
        const double*         a = rA.value_data().begin();
        
//        #pragma omp parallel for
        for(std::size_t i = 0; i < n; ++i)
        {
            int nz = ia[i + 1] - ia[i];
            for(int k = 0; k < nz; ++k)
            {
                double v = a[ia[i] + k];
                double col = ja[ia[i] + k];
                for(std::size_t j = 0; j < m; ++j)
                {
                    rC(i, j) += v*rB(col, j);
                }
            }
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

}; // Class  Block2PhaseSchurSolver

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::istream& operator >> (std::istream& IStream, Block2PhaseSchurSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis)
{
    return IStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator << (std::ostream& OStream, const  Block2PhaseSchurSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis)
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

#endif //  KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BLOCK_2_PHASE_SCHUR_SOLVER_H_INCLUDED  defined 

