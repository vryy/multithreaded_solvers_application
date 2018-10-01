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

#if !defined(KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BLOCK_2_PHASE_SCHUR_PRECONDITIONER_H_INCLUDED)
#define KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BLOCK_2_PHASE_SCHUR_PRECONDITIONER_H_INCLUDED




// System includes



// External includes
#include <boost/smart_ptr.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/progress.hpp>


// Project includes
#include "includes/define.h"
#include "linear_solvers/linear_solver.h"
#include "utilities/openmp_utils.h"

// #define STRINGIFY(name) #name
// #define SCHUR_DIAGONAL          0
// #define SCHUR_DIAGONAL_LUMPING  1
// #define SCHUR_PRECONDITIONER    2
// #define SCHUR_PRESSURE          3
// #define SCHUR_SOLVER            4

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

/// Block2PhaseSchurPreconditioner class.
/**
REF: White, Borja
*/
template<class TSparseSpaceType, class TDenseSpaceType>
class Block2PhaseSchurPreconditioner : public Preconditioner<TSparseSpaceType, TDenseSpaceType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Block2PhaseSchurPreconditioner
    KRATOS_CLASS_POINTER_DEFINITION (Block2PhaseSchurPreconditioner);

    typedef Preconditioner<TSparseSpaceType, TDenseSpaceType> BaseType;

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType> LinearSolverType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef boost::shared_ptr<SparseMatrixType> SparseMatrixPointerType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef std::size_t  SizeType;

    typedef std::size_t  IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Block2PhaseSchurPreconditioner(
        typename BaseType::Pointer prec_A,
        typename BaseType::Pointer prec_S,
        const std::string& SchurComputeMode
    ) : BaseType(), mprec_A(prec_A), mprec_S(prec_S), mSchurComputeMode(SchurComputeMode)
    , mInverseOption("SL")
    , mpS(SparseMatrixPointerType(new SparseMatrixType()))
    {
        if(mSchurComputeMode == std::string("SCHUR_SOLVER"))
            KRATOS_THROW_ERROR(std::logic_error, "In this Schur complement compute mode, a solver must be provided", "")
    }

    Block2PhaseSchurPreconditioner(
        typename BaseType::Pointer prec_A,
        typename BaseType::Pointer prec_S,
        const std::string& SchurComputeMode,
        const std::string& InverseOption
    ) : BaseType(), mprec_A(prec_A), mprec_S(prec_S), mSchurComputeMode(SchurComputeMode)
    , mInverseOption(InverseOption)
    , mpS(SparseMatrixPointerType(new SparseMatrixType()))
    {
        if(mSchurComputeMode == std::string("SCHUR_SOLVER"))
            KRATOS_THROW_ERROR(std::logic_error, "In this Schur complement compute mode, a solver must be provided", "")
    }

    Block2PhaseSchurPreconditioner(
        typename BaseType::Pointer prec_A,
        typename BaseType::Pointer prec_S,
        const std::string& SchurComputeMode,
        typename LinearSolverType::Pointer solver_S
    ) : BaseType(), mprec_A(prec_A), mprec_S(prec_S), mSchurComputeMode(SchurComputeMode)
    , msolver_S(solver_S), mInverseOption("SL")
    , mpS(SparseMatrixPointerType(new SparseMatrixType()))
    {
        if(mSchurComputeMode != std::string("SCHUR_SOLVER"))
            KRATOS_THROW_ERROR(std::logic_error, "Only SCHUR_SOLVER is supported", "")
    }

    Block2PhaseSchurPreconditioner(
        typename BaseType::Pointer prec_A,
        typename BaseType::Pointer prec_S,
        const std::string& SchurComputeMode,
        typename LinearSolverType::Pointer solver_S,
        const std::string& InverseOption
    ) : BaseType(), mprec_A(prec_A), mprec_S(prec_S), mSchurComputeMode(SchurComputeMode)
    , msolver_S(solver_S), mInverseOption(InverseOption)
    , mpS(SparseMatrixPointerType(new SparseMatrixType()))
    {
        if(mSchurComputeMode != std::string("SCHUR_SOLVER"))
            KRATOS_THROW_ERROR(std::logic_error, "Only SCHUR_SOLVER is supported", "")
    }

    /// Copy constructor.
    Block2PhaseSchurPreconditioner(const Block2PhaseSchurPreconditioner& Other)
    {
        mprec_A = Other.mprec_A;
        mprec_S = Other.mprec_S;
        mSchurComputeMode = Other.mSchurComputeMode;
        mInverseOption = Other.mInverseOption;
    }


    /// Destructor.
    virtual ~Block2PhaseSchurPreconditioner()
    {
    }


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    Block2PhaseSchurPreconditioner& operator=(const Block2PhaseSchurPreconditioner& Other)
    {
        mprec_A = Other.mprec_A;
        mprec_S = Other.mprec_S;
        mSchurComputeMode = Other.mSchurComputeMode;
        mInverseOption = Other.mInverseOption;
        return *this;
    }


    ///@}
    ///@name Operations
    ///@{

    void SetSchurMatrix(SparseMatrixPointerType pS)
    {
        mpS = pS;
        mSchurComputeMode = "SCHUR_GIVEN";
    }

    virtual void Initialize(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
        std::cout << "Fill blocks begin" << std::endl;
        double start = OpenMPUtils::GetCurrentTime();
        MultithreadedSolversMathUtils::FillBlockMatrices(rA,
            mother_indices, mpressure_indices, mglobal_to_local_indexing, mis_pressure_block,
            mA, mB1, mB2, mC);
        KRATOS_WATCH(MultithreadedSolversMathUtils::ComputeFrobeniusNorm(mA))
        KRATOS_WATCH(MultithreadedSolversMathUtils::ComputeFrobeniusNorm(mB1))
        KRATOS_WATCH(MultithreadedSolversMathUtils::ComputeFrobeniusNorm(mB2))
        KRATOS_WATCH(MultithreadedSolversMathUtils::ComputeFrobeniusNorm(mC))
        // KRATOS_WATCH(mA)
        // KRATOS_WATCH(mB1)
        // KRATOS_WATCH(mB2)
        // KRATOS_WATCH(mC)
        // WriteMatrixMarketMatrix("A.mm", mA, false);
        // WriteMatrixMarketMatrix("B1.mm", mB1, false);
        // WriteMatrixMarketMatrix("B2.mm", mB2, false);
        // WriteMatrixMarketMatrix("C.mm", mC, false);
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

//        KRATOS_WATCH(MultithreadedSolversMathUtils::ComputeFrobeniusNorm(mA))
//        KRATOS_WATCH(MultithreadedSolversMathUtils::ComputeFrobeniusNorm(mB1))
//        KRATOS_WATCH(MultithreadedSolversMathUtils::ComputeFrobeniusNorm(mB2))
//        KRATOS_WATCH(MultithreadedSolversMathUtils::ComputeFrobeniusNorm(mC))

//        mprec_A->Initialize(mA, mu, mru);
        mprec_A->Initialize(mA, mu, rB); //take rB as temporary, but it should not be
        std::cout << "mprec_A is initialized" << std::endl;
        // KRATOS_WATCH(mSchurComputeMode)
        // KRATOS_WATCH(mpS)
        // KRATOS_WATCH(mpS->size1())
        // KRATOS_WATCH(mpS->size2())
        // // KRATOS_WATCH((*mpS)(0, 0))
        // KRATOS_WATCH(*mpS)

        if (mSchurComputeMode != std::string("SCHUR_GIVEN"))
        {
            this->CalculateSchurComplement(*mpS, mA, mB1, mB2, mC);
            std::cout << "Schur complement is computed" << std::endl;
        }
        KRATOS_WATCH(MultithreadedSolversMathUtils::ComputeFrobeniusNorm(*mpS))
        // WriteMatrixMarketMatrix("S.mm", *mpS, false);

//        mprec_S->Initialize(*mpS, mp, mrp);
        mprec_S->Initialize(*mpS, mp, rB); //take rB as temporary, but it should not be
        std::cout << "mprec_S is initialized" << std::endl;

        std::cout << "Block preconditioner is initialized" << std::endl;

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

    virtual bool AdditionalPhysicalDataIsNeeded()
    {
        return false;
    }


    virtual void ProvideAdditionalData(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB,
        typename ModelPart::DofsArrayType& rdof_set,
        ModelPart& r_model_part
    )
    {}


    virtual void Mult(SparseMatrixType& rA, VectorType& rX, VectorType& rY)
    {
        TSparseSpaceType::Mult(rA, rX, rY);
        ApplyLeft(rY);
    }


    /** calculate preconditioned_u = A^{-1} * mu; preconditioned_p = S^{-1} * (mp - B2*preconditioned_u)
        @param rX  Unknowns of preconditioner suystem
    */
    virtual VectorType& ApplyLeft(VectorType& rX)
    {
        MultithreadedSolversMathUtils::GetPart(rX, mother_indices, mu);
        MultithreadedSolversMathUtils::GetPart(rX, mpressure_indices, mp);

        if (mInverseOption == std::string("SL")) // lower triangle inversion
        {
            mprec_A->ApplyLeft(mu); // mu <- A^-1 u

            VectorType paux(mp.size());

            TSparseSpaceType::Mult(mB2, mu, paux); // paux <- B2 A^-1 u

            TSparseSpaceType::UnaliasedAdd(mp, -1.0, paux); // mp <- p - B2 A^-1 u

            mprec_S->ApplyLeft(mp); // mp <- S^-1(p - B2 A^-1 u)
        }
        else if (mInverseOption == std::string("SU")) // upper triangle inversion
        {
            mprec_S->ApplyLeft(mp); // mp <- S^-1 p

            VectorType uaux(mu.size());

            TSparseSpaceType::Mult(mB1, mp, uaux); // paux <- B1 S^-1 p

            TSparseSpaceType::UnaliasedAdd(mu, -1.0, uaux); // mu <- u - B1 S^-1 p

            mprec_A->ApplyLeft(mu); // mu <- A^-1(u - B1 S^-1 p)
        }
        else if (mInverseOption == std::string("SF")) // full inversion
        {
            mprec_A->ApplyLeft(mu); // mu <- A^-1 u

            VectorType paux(mp.size());

            TSparseSpaceType::Mult(mB2, mu, paux); // paux <- B2 A^-1 u

            TSparseSpaceType::UnaliasedAdd(mp, -1.0, paux); // mp <- p - B2 A^-1 u

            mprec_S->ApplyLeft(mp); // mp <- S^-1(p - B2 A^-1 u)

            VectorType uaux(mu.size());

            TSparseSpaceType::Mult(mB1, mp, uaux); // uaux <- B1 S^-1(p - B2 A^-1 u)

            mprec_A->ApplyLeft(uaux); // uax <- A^-1 B1 S^-1(p - B2 A^-1 u)

            TSparseSpaceType::UnaliasedAdd(mu, -1.0, uaux); // mu <- A^-1 u - A^-1 B1 S^-1(p - B2 A^-1 u)
        }
        else
        {
            KRATOS_THROW_ERROR(std::logic_error, "Unknown option", mInverseOption)
        }

        MultithreadedSolversMathUtils::WritePart(rX, mother_indices, mu);
        MultithreadedSolversMathUtils::WritePart(rX, mpressure_indices, mp);

        return rX;
    }



    ///@}
    ///@name Access
    ///@{

    void SetSchurComputeMode(std::string SchurComputeMode)
    {
        mSchurComputeMode = SchurComputeMode;
        KRATOS_WATCH("SchurComputeMode has been changed");
        KRATOS_WATCH(mSchurComputeMode);
    }

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
        buffer << "Block2PhaseSchurPreconditioner";
        buffer << ", SchurComputeMode = " << mSchurComputeMode;
        buffer << ", InverseOption = " << mInverseOption;
        buffer << ", P_A = {" << mprec_A->Info() << "}";
        buffer << ", P_S = {" << mprec_S->Info() << "}";
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

    std::vector<SizeType> mpressure_indices;
    std::vector<SizeType> mother_indices;
    std::vector<SizeType> mglobal_to_local_indexing;
    std::vector<int> mis_pressure_block;

    typename BaseType::Pointer mprec_A;
    typename BaseType::Pointer mprec_S;
    typename LinearSolverType::Pointer msolver_S;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    void CalculateSchurComplement(
        SparseMatrixType& S,
        SparseMatrixType& A,
        SparseMatrixType& B1,
        SparseMatrixType& B2,
        SparseMatrixType& C)
    {
        KRATOS_TRY

        TSparseSpaceType::Resize(S, mpressure_indices.size(), mpressure_indices.size());

        if(mSchurComputeMode == std::string("SCHUR_DIAGONAL"))
        {
            //allocate the Schur complement
            // ConstructSystemMatrix(S, B1, B2, C);
            // calculate the Schur by diagonal approximant of A
            Vector approxA(mother_indices.size());
            ComputeDiagonalByExtracting(A, approxA);
            //fill the Schur complement
            // CalculateSchurComplementByDiagonalApproximation(S, A, B1, B2, C, approxA);
            CalculateSchurComplementByDiagonalApproximationSimple(S, A, B1, B2, C, approxA);
        }
        else if(mSchurComputeMode == std::string("SCHUR_DIAGONAL_LUMPING"))
        {
            //allocate the Schur complement
            // ConstructSystemMatrix(S, B1, B2, C);
            // calculate the Schur by diagonal approximant of A
            Vector approxA(mother_indices.size());
            ComputeDiagonalByLumping(A, approxA);
            //fill the Schur complement
            // CalculateSchurComplementByDiagonalApproximation(S, A, B1, B2, C, approxA);
            CalculateSchurComplementByDiagonalApproximationSimple(S, A, B1, B2, C, approxA);
        }
        else if(mSchurComputeMode == std::string("SCHUR_PRECONDITIONER"))
        {
            // calculate the Schur using the preconditioner (note that this may be very expensive, use it to test the algorithm only. Another point is that this may generate dense matrix which is not optimal for a sparse solver)
            CalculateSchurComplementByPreconditioner(S, A, B1, B2, C, mprec_A);
        }
        else if(mSchurComputeMode == std::string("SCHUR_SOLVER"))
        {
            // calculate the Schur using the solver (note that this may be very expensive, use it to test the algorithm only. Another point is that this may generate dense matrix which is not optimal for a sparse solver)
            CalculateSchurComplementBySolver(S, A, B1, B2, C, msolver_S);
        }
        else
        {
            KRATOS_THROW_ERROR(std::logic_error, "Unknown Schur compute mode: ", mSchurComputeMode);
        }

        KRATOS_CATCH ("")
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
    SparseMatrixPointerType mpS;

    VectorType mp;
    VectorType mu;

    std::string mSchurComputeMode;
    std::string mInverseOption;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /**
     * Compute the Schur complement A = L - D*Inv(Diag(K))*G. The multiplication
     * is performed in random order, so each row will be stored in a temporary
     * variable, ordered and copied in input matrix A.
     */
    void CalculateSchurComplementByDiagonalApproximation (
        SparseMatrixType& A,
        SparseMatrixType& K,
        SparseMatrixType& rG,
        SparseMatrixType& rD,
        SparseMatrixType& rL,
        VectorType& diagK
    ) const
    {
        // Retrieve matrices

        // Compute Inv(Diag(K))
        VectorType& rIDiagS = diagK;

        typedef boost::numeric::ublas::vector<int> IndexVector;
        typedef typename SparseMatrixType::iterator1 OuterIt;
        typedef typename SparseMatrixType::iterator2 InnerIt;
        typedef typename boost::numeric::ublas::matrix_row< SparseMatrixType > RowType;

        int DiagSize = int (diagK.size()); // to avoid comparison between int & unsigned int
        #pragma omp parallel for
        for ( int i = 0; i < DiagSize; ++i)
            rIDiagS[i] = 1.0 / diagK[i];
        OpenMPUtils::PartitionVector Partition;
        int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::DivideInPartitions (A.size1(),NumThreads,Partition);
        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();
            VectorType CurrentRow(K.size2());

            for (unsigned int i = 0; i < rL.size1(); ++i) CurrentRow[i] = 0.0;
            boost::shared_ptr< IndexVector > pNext ( new IndexVector (rL.size1()) );
            IndexVector& Next = *pNext; // Keeps track of which columns were filled
            for (unsigned int m=0; m < rL.size1(); ++m) Next[m] = -1;
            SizeType NumTerms = 0; // Full positions in a row
            boost::shared_ptr< std::vector<unsigned int> > pUsedCols(new std::vector<unsigned int>);
            std::vector<unsigned int>& UsedCols = *pUsedCols;
            UsedCols.reserve (rL.size1());
            for ( int RowIndex = Partition[k] ; RowIndex != Partition[k+1] ; ++RowIndex )
            {
                RowType RowD (rD,RowIndex);
                RowType RowL (rL,RowIndex);
                int head = -2;
                SizeType Length = 0;
                // Write L in A
                for ( typename RowType::iterator ItL = RowL.begin(); ItL != RowL.end(); ++ItL )
                {
                    CurrentRow (ItL.index() ) = *ItL;
                    if ( Next[ItL.index()] == -1)
                    {
                        Next[ItL.index()] = head;
                        head = ItL.index();
                        ++Length;
                    }
                }
                // Substract D*Inv(Diag(S))*G
                for ( typename RowType::iterator ItD = RowD.begin(); ItD != RowD.end(); ++ItD )
                {
                    RowType RowG (rG,ItD.index() );
                    for ( typename RowType::iterator ItG = RowG.begin(); ItG != RowG.end(); ++ItG )
                    {
                        CurrentRow[ItG.index()] -= (*ItD) * rIDiagS[ItD.index()] * (*ItG);
                        if ( Next[ItG.index()] == -1)
                        {
                            Next[ItG.index()] = head;
                            head = ItG.index();
                            ++Length;
                        }
                    }
                }
                // Identify full terms for ordering
                for ( SizeType i = 0; i < Length; ++i)
                {
                    if ( Next[head] != -1 )
                    {
                        UsedCols.push_back (head);
                        ++NumTerms;
                    }
                    int temp = head;
                    head = Next[head];
                    // Clear 'Next' for next iteration
                    Next[temp] = -1;
                }
                // Sort Column indices
                SortCols (UsedCols,NumTerms);
                // Fill matrix row, then clean temporary variables.
                RowType RowA (A,RowIndex);
                SizeType n = 0;
                unsigned int Col;
                for ( typename RowType::iterator ItA = RowA.begin(); ItA != RowA.end(); ++ItA)
                {
                    Col = UsedCols[n++];
                    *ItA = CurrentRow[Col];
                    CurrentRow[Col] = 0;
                }
                NumTerms = 0;
                UsedCols.resize (0,false);
            }
        }

        //KRATOS_WATCH(896)
        //add stabilization matrix L
        /*				const SizeType* L_index1 = rL.index1_data().begin();
        				const SizeType* L_index2 = rL.index2_data().begin();
        				const double*	   L_values = rL.value_data().begin();
        				for (unsigned int i=0; i<rL.size1(); i++)
        				{
        					unsigned int row_begin = L_index1[i];
        					unsigned int row_end   = L_index1[i+1];
        					diagA[i] = 0.0;
        					for (unsigned int j=row_begin; j<row_end; j++)
        					{
        						unsigned int col = L_index2[j];
        						rS(i,col) += L_values[j];
        					}
        				}*/
    }

    /**
     * Compute the Schur complement S = C - B2*Inv(Diag(A))*B1.
     */
    void CalculateSchurComplementByDiagonalApproximationSimple (
        SparseMatrixType& S,
        SparseMatrixType& A,
        SparseMatrixType& B1,
        SparseMatrixType& B2,
        SparseMatrixType& C,
        VectorType& diagA
    ) const
    {
        VectorType invDiagA(diagA.size());
        for(std::size_t i = 0; i < diagA.size(); ++i)
            invDiagA[i] = 1.0 / diagA(i);
        // KRATOS_WATCH(invDiagA)
        SparseMatrixType tmpS = B1;
        MultithreadedSolversMathUtils::RowScale(tmpS, invDiagA);
        MultithreadedSolversMathUtils::MatrixMult(S, B2, tmpS);
        S *= -1.0;
        noalias(S) += C;
    }

    void CalculateSchurComplementByPreconditioner(
        SparseMatrixType& S,
        SparseMatrixType& A,
        SparseMatrixType& B1,
        SparseMatrixType& B2,
        SparseMatrixType& C,
        typename BaseType::Pointer prec_A_Schur) const
    {
        typename TDenseSpaceType::MatrixType Tmp(TSparseSpaceType::Size1(B2), TSparseSpaceType::Size2(B1));
        boost::progress_display show_progress(TSparseSpaceType::Size2(B1));
        for(IndexType i = 0; i < TSparseSpaceType::Size2(B1); ++i)
        {
            VectorType b1(TSparseSpaceType::Size1(B1));

            TSparseSpaceType::GetColumn(i, B1, b1);

            b1 = prec_A_Schur->ApplyLeft(b1);

            VectorType b2(TSparseSpaceType::Size1(B2));

            TSparseSpaceType::Mult(B2, b1, b2);

            noalias(column(Tmp, i)) = b2;
            ++show_progress;
        }
        std::cout << "Schur complement using preconditioner done" << std::endl;
        noalias(S) = C - Tmp;

//        KRATOS_WATCH(A)
//        KRATOS_WATCH(B1)
//        KRATOS_WATCH(B2)
//        KRATOS_WATCH(C)
//        KRATOS_WATCH(S)
    }

    void CalculateSchurComplementBySolver(
        SparseMatrixType& S,
        SparseMatrixType& A,
        SparseMatrixType& B1,
        SparseMatrixType& B2,
        SparseMatrixType& C,
        typename LinearSolverType::Pointer solver_S) const
    {
        std::cout << "compute Schur complement using solver" <<  std::endl;
        DenseMatrixType aux1(A.size1(), B1.size2());
        DenseMatrixType aux2 = B1;
        solver_S->Solve(A, aux1, aux2);
        std::cout << "solver_S->Solve completed" <<  std::endl;
        noalias(S) = C - prod(B2, aux1);
        std::cout << "Schur complement using solver done" << std::endl;
    }

    /// Helper function for Schur complement functions
    void SortCols (std::vector<unsigned int>& ColList, SizeType& NumCols) const
    {
        bool swap = true;
        unsigned int d = NumCols;
        int temp;
        while ( swap || d > 1 )
        {
            swap = false;
            d = (d+1) / 2;
            for ( unsigned int i=0; i< (NumCols - d); ++i)
                if ( ColList[i+d] < ColList[i] )
                {
                    temp = ColList[i+d];
                    ColList[i+d] = ColList[i];
                    ColList[i] = temp;
                    swap = true;
                }
        }
    }

    /// Identify non-zero tems in the Schur complement
    void ConstructSystemMatrix(
        SparseMatrixType& A,
        SparseMatrixType& rG,
        SparseMatrixType& rD,
        SparseMatrixType& rL
    ) const
    {
        typedef boost::numeric::ublas::vector<int> IndexVector;
        typedef OpenMPUtils::PartitionVector PartitionVector;
        typedef typename SparseMatrixType::iterator1 OuterIt;
        typedef typename SparseMatrixType::iterator2 InnerIt;
        typedef typename boost::numeric::ublas::matrix_row< SparseMatrixType > RowType;

        PartitionVector Partition;
        int NumThreads = OpenMPUtils::GetNumThreads();

        OpenMPUtils::DivideInPartitions(A.size1(),NumThreads,Partition);

        for (int k = 0 ; k < NumThreads ; ++k)
        {
            // This code is serial, the pragma is here to ensure that each
            // row block is assigned to the processor that will fill it
            #pragma omp parallel
            if ( OpenMPUtils::ThisThread() == k)
            {
                boost::shared_ptr< IndexVector > pNext( new IndexVector(rL.size1() ) );
                IndexVector& Next = *pNext; // Keeps track of which columns were filled
                for (unsigned int m = 0; m < rL.size1(); m++) Next[m] = -1;

                SizeType NumTerms = 0; // Full positions in a row
                boost::shared_ptr< std::vector<unsigned int> > pUsedCols( new std::vector<unsigned int>);
                std::vector<unsigned int>& UsedCols = *pUsedCols;
                UsedCols.reserve(rL.size1());

                for ( int RowIndex = Partition[k] ; RowIndex != Partition[k+1] ; ++RowIndex )
                {
                    RowType RowD(rD,RowIndex);
                    RowType RowL(rL,RowIndex);

                    int head = -2;
                    SizeType Length = 0;

                    // Terms filled by L
                    for ( typename RowType::iterator ItL = RowL.begin(); ItL != RowL.end(); ++ItL )
                    {
                        if ( Next[ItL.index()] == -1)
                        {
                            Next[ItL.index()] = head;
                            head = ItL.index();
                            Length++;
                        }
                    }

                    // Additional terms due to D*Inv(Diag(S))*G
                    for ( typename RowType::iterator ItD = RowD.begin(); ItD != RowD.end(); ++ItD )
                    {
                        RowType RowG(rG,ItD.index());

                        for ( typename RowType::iterator ItG = RowG.begin(); ItG != RowG.end(); ++ItG )
                        {
                            if ( Next[ItG.index()] == -1)
                            {
                                Next[ItG.index()] = head;
                                head = ItG.index();
                                Length++;
                            }
                        }
                    }

                    // Identify full terms for ordering
                    for ( SizeType i = 0; i < Length; ++i)
                    {
                        if ( Next[head] != -1 )
                        {
                            UsedCols.push_back(head);
                            NumTerms++;
                        }

                        int temp = head;
                        head = Next[head];

                        // Clear 'Next' for next iteration
                        Next[temp] = -1;
                    }

                    // Sort Column indices
                    SortCols(UsedCols,NumTerms);

                    // Store row in matrix, clean temporary variables
                    for ( unsigned int i = 0; i < NumTerms; ++i)
                    {
                        A.push_back(RowIndex,UsedCols[i],0);
                    }
                    NumTerms = 0;
                    UsedCols.resize(0,false);
                }
            }
        }
    }

    void ComputeDiagonalByLumping (SparseMatrixType& A, VectorType& diagA) const
    {
        if (diagA.size() != A.size1() )
            TSparseSpaceType::Resize(diagA, A.size1());
        //get access to A data
        const SizeType* index1 = A.index1_data().begin();
//        const SizeType* index2 = A.index2_data().begin();
        const double*	   values = A.value_data().begin();

        #pragma omp parallel for
        for (SizeType i = 0; i < A.size1(); ++i)
        {
            unsigned int row_begin = index1[i];
            unsigned int row_end   = index1[i+1];
            double temp = 0.0;
            for (unsigned int j = row_begin; j < row_end; ++j)
                temp += values[j] * values[j];

            diagA[i] = sqrt(temp);
        }
    }

    void ComputeDiagonalByExtracting (SparseMatrixType& A, VectorType& diagA) const
    {
        if (diagA.size() != A.size1() )
            TSparseSpaceType::Resize(diagA, A.size1());
        //get access to A data
        const SizeType* index1 = A.index1_data().begin();
        const SizeType* index2 = A.index2_data().begin();
        const double*	   values = A.value_data().begin();

        #pragma omp parallel for
        for (SizeType i = 0; i < A.size1(); ++i)
        {
            unsigned int row_begin = index1[i];
            unsigned int row_end   = index1[i+1];

            unsigned int jj = 0;
            bool has_diagonal = false;
            for (unsigned int j = row_begin; j < row_end; ++j)
            {
                if(index2[j] == i)
                {
                    jj = j;
                    has_diagonal = true;
                    break;
                }
            }

            if( has_diagonal )
                diagA[i] = values[jj];
            else
                diagA[i] = 0.0;
        }
    }

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{

    // std::string SchurComputeModeNames(int SchurComputeMode) const
    // {
    //     switch(SchurComputeMode)
    //     {
    //         case SCHUR_DIAGONAL:
    //             return std::string("Schur_Diagonal");
    //         case SCHUR_DIAGONAL_LUMPING:
    //             return std::string("Schur_DiagonalLumping");
    //         case SCHUR_PRECONDITIONER:
    //             return std::string("Schur_Preconditioner");
    //         case SCHUR_PRESSURE:
    //             return std::string("Schur_Pressure");
    //         default:
    //             return std::string("None");
    //     }
    // }

    ///@}
    ///@name Unaccessible methods
    ///@{


    ///@}

}; // Class Block2PhaseSchurPreconditioner

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
                                  Block2PhaseSchurPreconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
{
    return IStream;
}


/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType>
inline std::ostream& operator << (std::ostream& OStream,
                                  const Block2PhaseSchurPreconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
{
    rThis.PrintInfo(OStream);
    OStream << std::endl;
    rThis.PrintData(OStream);


    return OStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BLOCK_2_PHASE_SCHUR_PRECONDITIONER_H_INCLUDED  defined
