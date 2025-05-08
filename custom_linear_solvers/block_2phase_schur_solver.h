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


#if !defined(KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BLOCK_2_PHASE_SCHUR_SOLVER_H_INCLUDED )
#define KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BLOCK_2_PHASE_SCHUR_SOLVER_H_INCLUDED


// System includes
#include <cmath>


// External includes


// Project includes
#include "includes/define.h"
#include "linear_solvers/linear_solver.h"
#include "custom_utilities/multithreaded_solvers_math_utils.h"


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
         class TModelPartType,
         class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class Block2PhaseSchurSolver : public LinearSolver<TSparseSpaceType, TDenseSpaceType, TModelPartType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of  Block2PhaseSchurSolver
    KRATOS_CLASS_POINTER_DEFINITION( Block2PhaseSchurSolver );

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TModelPartType, TReordererType> BaseType;

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
    Block2PhaseSchurSolver() {}

    Block2PhaseSchurSolver(typename BaseType::Pointer pSolver)
    : BaseType(), mpSolverA(pSolver), mpSolverS(pSolver)
    {}

    Block2PhaseSchurSolver(typename BaseType::Pointer pSolverA,
        typename BaseType::Pointer pSolverS)
    : BaseType(), mpSolverA(pSolverA), mpSolverS(pSolverS)
    {}

    /// Copy constructor.
    Block2PhaseSchurSolver(const  Block2PhaseSchurSolver& Other)
    : BaseType(Other) {}

    /// Destructor.
    ~Block2PhaseSchurSolver() override {}


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
        mpSolverA->ProvideAdditionalData(rA, rX, rB, rdof_set, r_model_part);
        mpSolverS->ProvideAdditionalData(rA, rX, rB, rdof_set, r_model_part);
    }

    void Initialize(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        std::cout << "Fill blocks begin" << std::endl;
        double start = OpenMPUtils::GetCurrentTime();
        MultithreadedSolversMathUtils::FillBlockMatrices(rA,
            mother_indices, mpressure_indices, mglobal_to_local_indexing, mis_pressure_block,
            mA, mB1, mB2, mC);
        std::cout << "Fill blocks completed..." << OpenMPUtils::GetCurrentTime() - start << " s" << std::endl;

        TSparseSpaceType::Resize(mp, mpressure_indices.size());
        TSparseSpaceType::Set(mp, 0.00);
        TSparseSpaceType::Resize(mu, mother_indices.size());
        TSparseSpaceType::Set(mu, 0.00);

        //this is rather slow
//        KRATOS_WATCH(norm_frobenius(mA))
//        KRATOS_WATCH(norm_frobenius(mB1))
//        KRATOS_WATCH(norm_frobenius(mB2))
//        KRATOS_WATCH(norm_frobenius(mC))

        KRATOS_WATCH(MultithreadedSolversMathUtils::ComputeFrobeniusNorm(mA))
        KRATOS_WATCH(MultithreadedSolversMathUtils::ComputeFrobeniusNorm(mB1))
        KRATOS_WATCH(MultithreadedSolversMathUtils::ComputeFrobeniusNorm(mB2))
        KRATOS_WATCH(MultithreadedSolversMathUtils::ComputeFrobeniusNorm(mC))

        // TODO: check the diagonal dominance of matrix A
        #ifdef CHECK_DIAGONAL_DOMINANCE
        MultithreadedSolversMathUtils::CheckDiagonalDominance(mA);
        #endif

        mpSolverA->Initialize(mA, mu, rB); //take rB as temporary, but it should not be
        mpSolverS->Initialize(mC, mp, rB); //take rB as temporary, but it should not be
    }

    /** Normal solve method.
    Solves the linear system Ax=b and puts the result on SystemVector& rX.
    rX is also th initial guess for iterative methods.
    @param rA. System matrix
    @param rX. Solution vector. it's also the initial
    guess for iterative linear solvers.
    @param rB. Right hand side vector.
    */
    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        if(this->IsNotConsistent(rA, rX, rB))
            return false;

        Initialize(rA, rX, rB);
        #ifdef ENABLE_MORE_OUTPUT
        std::cout << "Initialize completed" << std::endl;
        KRATOS_WATCH(rA.size1())
        KRATOS_WATCH(rA.size2())
        #endif

        VectorType ru, rp, u, p;

        // Get the initial u & p
        MultithreadedSolversMathUtils::GetPart(rX, mother_indices, u);
        MultithreadedSolversMathUtils::GetPart(rX, mpressure_indices, p);

        // Get ru, rp
        MultithreadedSolversMathUtils::GetPart(rB, mother_indices, ru);
        MultithreadedSolversMathUtils::GetPart(rB, mpressure_indices, rp);

        // KRATOS_WATCH(u.size())
        // KRATOS_WATCH(p.size())
        // KRATOS_WATCH(ru.size())
        // KRATOS_WATCH(rp.size())

        // Extract inverse of diagonal of mA
        std::size_t n = mA.size1();
        // KRATOS_WATCH(n)
        VectorType invDiagBlockA(n);
        for(std::size_t i = 0; i < n; ++i)
            invDiagBlockA[i] = 1.0 / mA(i, i);
//        KRATOS_WATCH(invDiagBlockA)
//        KRATOS_WATCH(ru)
//        KRATOS_WATCH(rp)

        // Solve for p
        VectorType tmp_u = ru;
        // KRATOS_WATCH(norm_2(tmp_u))
//        KRATOS_WATCH(tmp_u)
        // KRATOS_WATCH(norm_2(invDiagBlockA))
        MultithreadedSolversMathUtils::VectorScale(tmp_u, invDiagBlockA);
//        KRATOS_WATCH(tmp_u)
        // KRATOS_WATCH(norm_2(tmp_u))
        Vector aux1(mB2.size1());
        TSparseSpaceType::Mult(mB2, tmp_u, aux1);
        // KRATOS_WATCH(norm_2(aux1))
        VectorType rhs_p = rp - aux1;
        // KRATOS_WATCH(norm_2(rhs_p))
//        KRATOS_WATCH(rhs_p.size())
//        KRATOS_WATCH(rp.size())
//        KRATOS_WATCH(mB2.size1())
//        KRATOS_WATCH(mB2.size2())
//        KRATOS_WATCH(rhs_p.size())
        #ifdef ENABLE_MORE_OUTPUT
        std::cout << "Create RHS for p completed" << std::endl;
        #endif

        SparseMatrixType tmpS = mB1;
        // KRATOS_WATCH(__LINE__)
        MultithreadedSolversMathUtils::RowScale(tmpS, invDiagBlockA);
        // KRATOS_WATCH(__LINE__)
        DenseMatrixType S(mC.size1(), mC.size2());
        // KRATOS_WATCH(__LINE__)
        MultithreadedSolversMathUtils::MatrixMult(S, mB2, static_cast<DenseMatrixType>(tmpS));
        // KRATOS_WATCH(__LINE__)
        S *= -1.0;
        // KRATOS_WATCH(__LINE__)
        noalias(S) += mC;
        #ifdef ENABLE_MORE_OUTPUT
        std::cout << "Compute Schur completed" << std::endl;
        KRATOS_WATCH(norm_frobenius(S))
        #endif

        SparseMatrixType sS = S;
        #ifdef ENABLE_MORE_OUTPUT
        KRATOS_WATCH(MultithreadedSolversMathUtils::ComputeFrobeniusNorm(sS))
        #endif
        mpSolverS->Solve(sS, p, rhs_p);
//        KRATOS_WATCH(p)
        #ifdef ENABLE_MORE_OUTPUT
        KRATOS_WATCH(norm_2(p))
        std::cout << "Solve for p completed" << std::endl;
        #endif

        // Solve for u
        Vector aux2(ru.size());
        TSparseSpaceType::Mult(mB1, p, aux2);
        noalias(ru) -= aux2;
        #ifdef ENABLE_MORE_OUTPUT
        KRATOS_WATCH(norm_2(ru))
        #endif
        mpSolverA->Solve(mA, u, ru);
//        KRATOS_WATCH(u)
        #ifdef ENABLE_MORE_OUTPUT
        KRATOS_WATCH(norm_2(u))
        std::cout << "Solve for u completed" << std::endl;
        #endif

        // Write back the result to solution vector
        MultithreadedSolversMathUtils::WritePart(rX, mother_indices, u);
        MultithreadedSolversMathUtils::WritePart(rX, mpressure_indices, p);

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
    bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB) override
    {
        KRATOS_ERROR << "Multisolve is not yet supported for " << typeid(*this).name();
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
        buffer << "Linear solver using Schur complement reduction scheme for 2 phases";
        buffer << ", A linear solver = " << mpSolverA->Info();
        buffer << ", S linear solver = " << mpSolverS->Info();
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

    typename BaseType::Pointer mpSolverA;
    typename BaseType::Pointer mpSolverS;

    std::vector<SizeType> mpressure_indices;
    std::vector<SizeType> mother_indices;
    std::vector<SizeType> mglobal_to_local_indexing;
    std::vector<int> mis_pressure_block;

    typename TModelPartType::DofsArrayType madof_set;
    typename TModelPartType::DofsArrayType msdof_set;

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

}; // Class  Block2PhaseSchurSolver

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
#undef ENABLE_MORE_OUTPUT

#endif //  KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BLOCK_2_PHASE_SCHUR_SOLVER_H_INCLUDED  defined
