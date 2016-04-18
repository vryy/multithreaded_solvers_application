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

#if !defined(KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BLOCK_PRESSURE_PRECONDITIONER_H_INCLUDED)
#define KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BLOCK_PRESSURE_PRECONDITIONER_H_INCLUDED




// System includes



// External includes
#include <boost/smart_ptr.hpp>
#include <boost/numeric/ublas/vector.hpp>


// Project includes
#include "includes/define.h"
#include "linear_solvers/linear_solver.h"
#include "utilities/openmp_utils.h"

#define STRINGIFY(name) #name
#define SCHUR_DIAGONAL          0
#define SCHUR_DIAGONAL_LUMPING  1
#define SCHUR_PRECONDITIONER    2
#define SCHUR_PRESSURE          3
#define SCHUR_SOLVER            4

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

/// BlockPressurePreconditioner class.
/**
REF: White, Borja
*/
template<class TSparseSpaceType, class TDenseSpaceType>
class BlockPressurePreconditioner : public Preconditioner<TSparseSpaceType, TDenseSpaceType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BlockPressurePreconditioner
    KRATOS_CLASS_POINTER_DEFINITION (BlockPressurePreconditioner);

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

    /// Default constructor.
    BlockPressurePreconditioner(
        typename BaseType::Pointer prec_A,
        typename BaseType::Pointer prec_S,
        std::string SchurComputeMode
    )
    {
        mprec_A = prec_A;
        mprec_S = prec_S;
        mSchurComputeMode = SchurComputeMode;
        if(mSchurComputeMode == std::string("SCHUR_SOLVER"))
            KRATOS_THROW_ERROR(std::logic_error, "In this Schur complement compute mode, a solver must be provided", "")
    }

    BlockPressurePreconditioner(
        typename BaseType::Pointer prec_A,
        typename BaseType::Pointer prec_S,
        std::string SchurComputeMode,
        typename LinearSolverType::Pointer solver_S
    )
    {
        mprec_A = prec_A;
        mprec_S = prec_S;
        mSchurComputeMode = SchurComputeMode;
        msolver_S = solver_S;
    }

    /// Copy constructor.
    BlockPressurePreconditioner(const BlockPressurePreconditioner& Other)
    {
        mprec_A = Other.mprec_A;
        mprec_S = Other.mprec_S;
        mSchurComputeMode = Other.mSchurComputeMode;
    }


    /// Destructor.
    virtual ~BlockPressurePreconditioner()
    {
    }


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    BlockPressurePreconditioner& operator=(const BlockPressurePreconditioner& Other)
    {
        mprec_A = Other.mprec_A;
        mprec_S = Other.mprec_S;
        return *this;
    }


    ///@}
    ///@name Operations
    ///@{

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

//        KRATOS_WATCH(ComputeFrobeniusNorm(mA))
//        KRATOS_WATCH(ComputeFrobeniusNorm(mB1))
//        KRATOS_WATCH(ComputeFrobeniusNorm(mB2))
//        KRATOS_WATCH(ComputeFrobeniusNorm(mC))

//        mprec_A->Initialize(mA, mu, mru);
        mprec_A->Initialize(mA, mu, rB); //take rB as temporary, but it should not be
        std::cout << "mprec_A is initialized" << std::endl;

        CalculateSchurComplement(mS, mA, mB1, mB2, mC);
        std::cout << "Schur complement is computed" << std::endl;

//        mprec_S->Initialize(mS, mp, mrp);
        mprec_S->Initialize(mS, mp, rB); //take rB as temporary, but it should not be
        std::cout << "mprec_B is initialized" << std::endl;

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
        
//        if(mprec_A->AdditionalPhysicalDataIsNeeded())
//            mprec_A->ProvideAdditionalData(rA, rX, rB, rdof_set, r_model_part); //TODO: the parameters need to be customized to the block A
//            
//        if(mprec_S->AdditionalPhysicalDataIsNeeded())
//            mprec_S->ProvideAdditionalData(rA, rX, rB, rdof_set, r_model_part); //TODO: the parameters need to be customized to the block S
    }


    virtual void Mult(SparseMatrixType& rA, VectorType& rX, VectorType& rY)
    {
        TSparseSpaceType::Mult(rA, rX, rY);
        ApplyLeft(rY);
    }


    /** calculate preconditioned_u = A^{-1} * mu; preconditioned_p = S^{-1} * (mp - B2*preconditioned_u)
        @param rX  Unknows of preconditioner suystem
    */
    virtual VectorType& ApplyLeft(VectorType& rX)
    {
        GetUPart(rX, mu);

        mprec_A->ApplyLeft(mu);

        WriteUPart(rX, mu);

        GetPPart(rX, mp);

        VectorType paux(mp.size());
        
        TSparseSpaceType::Mult(mB2, mu, paux);

        TSparseSpaceType::UnaliasedAdd(mp, -1.0, paux);
        
        mprec_S->ApplyLeft(mp);

        WritePPart(rX, mp);
        
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
        buffer << "BlockPressurePreconditioner";
        buffer << ", SchurComputeMode = " << mSchurComputeMode;
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

    void CalculateSchurComplement(
        SparseMatrixType& S,
        SparseMatrixType& A,
        SparseMatrixType& B1,
        SparseMatrixType& B2,
        SparseMatrixType& C)
    {
        KRATOS_TRY

        TSparseSpaceType::Resize(S, mpressure_indices.size(), mpressure_indices.size());

        //allocate the Schur complement
        ConstructSystemMatrix(S, B1, B2, C); //this is only effective for Schur_Diagonal & Schur_DiagonalLumping mode

        if(mSchurComputeMode == std::string("SCHUR_DIAGONAL"))
        {
            // calculate the Schur by diagonal approximant of A
            Vector approxA(mother_indices.size());
            ComputeDiagonalByExtracting(A, approxA);
            //fill the Schur complement
            CalculateSchurComplementByDiagonalApproximation(S, A, B1, B2, C, approxA);
        }
        else if(mSchurComputeMode == std::string("SCHUR_DIAGONAL_LUMPING"))
        {
            // calculate the Schur by diagonal approximant of A
            Vector approxA(mother_indices.size());
            ComputeDiagonalByLumping(A, approxA);
            //fill the Schur complement
            CalculateSchurComplementByDiagonalApproximation(S, A, B1, B2, C, approxA);
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

    typename BaseType::Pointer mprec_A;
    typename BaseType::Pointer mprec_S;
    typename LinearSolverType::Pointer msolver_S;
    
    SparseMatrixType mA;
    SparseMatrixType mB1;
    SparseMatrixType mB2;
    SparseMatrixType mC;
    SparseMatrixType mS;
    
    VectorType mp;
    VectorType mu;

    std::string mSchurComputeMode;
    
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
    )
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
            std::size_t NumTerms = 0; // Full positions in a row
            boost::shared_ptr< std::vector<unsigned int> > pUsedCols(new std::vector<unsigned int>);
            std::vector<unsigned int>& UsedCols = *pUsedCols;
            UsedCols.reserve (rL.size1());
            for ( int RowIndex = Partition[k] ; RowIndex != Partition[k+1] ; ++RowIndex )
            {
                RowType RowD (rD,RowIndex);
                RowType RowL (rL,RowIndex);
                int head = -2;
                std::size_t Length = 0;
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
                for ( std::size_t i = 0; i < Length; ++i)
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
                std::size_t n = 0;
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
        /*				const std::size_t* L_index1 = rL.index1_data().begin();
        				const std::size_t* L_index2 = rL.index2_data().begin();
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
    
    
    void CalculateSchurComplementByPreconditioner(
        SparseMatrixType& S,
        SparseMatrixType& A,
        SparseMatrixType& B1,
        SparseMatrixType& B2,
        SparseMatrixType& C,
        typename BaseType::Pointer prec_A_Schur)
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
        typename LinearSolverType::Pointer solver_S)
    {
        DenseMatrixType aux1(A.size1(), B1.size2());
        DenseMatrixType aux2 = B1;
        solver_S->Solve(A, aux1, aux2);
        std::cout << "solver_S->Solve completed" <<  std::endl;
        noalias(S) = C - prod(B2, aux1);
        std::cout << "Schur complement using solver done" << std::endl;
    }
    
    
    /// Helper function for Schur complement functions
    void SortCols (std::vector<unsigned int>& ColList, std::size_t& NumCols)
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
    )
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

                std::size_t NumTerms = 0; // Full positions in a row
                boost::shared_ptr< std::vector<unsigned int> > pUsedCols( new std::vector<unsigned int>);
                std::vector<unsigned int>& UsedCols = *pUsedCols;
                UsedCols.reserve(rL.size1());

                for ( int RowIndex = Partition[k] ; RowIndex != Partition[k+1] ; ++RowIndex )
                {
                    RowType RowD(rD,RowIndex);
                    RowType RowL(rL,RowIndex);

                    int head = -2;
                    std::size_t Length = 0;

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
                    for ( std::size_t i = 0; i < Length; ++i)
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
    
    
    void ComputeDiagonalByLumping (SparseMatrixType& A,VectorType& diagA)
    {
        if (diagA.size() != A.size1() )
            TSparseSpaceType::Resize(diagA, A.size1());
        //get access to A data
        const std::size_t* index1 = A.index1_data().begin();
//        const std::size_t* index2 = A.index2_data().begin();
        const double*	   values = A.value_data().begin();

        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(A.size1()); ++i)
        {
            unsigned int row_begin = index1[i];
            unsigned int row_end   = index1[i+1];
            double temp = 0.0;
            for (unsigned int j = row_begin; j < row_end; ++j)
                temp += values[j] * values[j];

            diagA[i] = sqrt(temp);
        }
    }
    
    
    void ComputeDiagonalByExtracting (SparseMatrixType& A,VectorType& diagA)
    {
        if (diagA.size() != A.size1() )
            TSparseSpaceType::Resize(diagA, A.size1());
        //get access to A data
        const std::size_t* index1 = A.index1_data().begin();
        const std::size_t* index2 = A.index2_data().begin();
        const double*	   values = A.value_data().begin();

        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(A.size1()); ++i)
        {
            unsigned int row_begin = index1[i];
            unsigned int row_end   = index1[i+1];
            
            unsigned int jj = 0;
            bool has_diagonal = false;
            for (unsigned int j=row_begin; j<row_end; ++j)
                if(index2[j] == i)
                {
                    jj = j;
                    has_diagonal = true;
                    break;
                }

            if( has_diagonal )
                diagA[i] = values[jj];
            else
                diagA[i] = 0.0;
        }
    }
    
    double ComputeFrobeniusNorm(SparseMatrixType& rA)
    {
        int n = rA.size1();
        const std::size_t* ia = rA.index1_data().begin();
        const std::size_t* ja = rA.index2_data().begin();
        const double*	   a  = rA.value_data().begin();
        
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
        return sqrt(norm);
    }
    
    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{

    std::string SchurComputeModeNames(int SchurComputeMode) const
    {
        switch(SchurComputeMode)
        {
            case SCHUR_DIAGONAL:
                return std::string("Schur_Diagonal");
            case SCHUR_DIAGONAL_LUMPING:
                return std::string("Schur_DiagonalLumping");
            case SCHUR_PRECONDITIONER:
                return std::string("Schur_Preconditioner");
            case SCHUR_PRESSURE:
                return std::string("Schur_Pressure");
            default:
                return std::string("None");
        }
    }

    ///@}
    ///@name Unaccessible methods
    ///@{


    ///@}

}; // Class BlockPressurePreconditioner

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
                                  BlockPressurePreconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
{
    return IStream;
}


/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType>
inline std::ostream& operator << (std::ostream& OStream,
                                  const BlockPressurePreconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
{
    rThis.PrintInfo(OStream);
    OStream << std::endl;
    rThis.PrintData(OStream);


    return OStream;
}
///@}


}  // namespace Kratos.


#endif // KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BLOCK_PRESSURE_PRECONDITIONER_H_INCLUDED  defined 

