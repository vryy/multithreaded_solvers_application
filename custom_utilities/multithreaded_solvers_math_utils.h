//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         finite_cell_application/LICENSE.txt
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Hoang-Giang Bui
//  Date:            25 Sep 2018
//


#if !defined(KRATOS_MULTITHREADED_SOLVERS_MATH_UTILITY_H_INCLUDED )
#define  KRATOS_MULTITHREADED_SOLVERS_MATH_UTILITY_H_INCLUDED



// System includes
#include <string>


// External includes


// Project includes
#include "includes/define.h"
#include "utilities/openmp_utils.h"


namespace Kratos
{
///@addtogroup FiniteCellApplication
///@{

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
/** class for auxilliary routines
*/
class MultithreadedSolversMathUtils
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MultithreadedSolversMathUtils
    KRATOS_CLASS_POINTER_DEFINITION(MultithreadedSolversMathUtils);

    typedef typename Element::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;

    typedef typename NodeType::CoordinatesArrayType CoordinatesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MultithreadedSolversMathUtils() {}

    /// Destructor.
    virtual ~MultithreadedSolversMathUtils() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    //this function extracts from a vector which has the size of the
    //overall r, the part that corresponds to u-dofs
    template<typename TVectorType, typename TIndexType = std::size_t>
    static void GetPart (const TVectorType& rtot, const std::vector<TIndexType>& indices, TVectorType& ru)
    {
        const TIndexType size = indices.size();
        if (ru.size() != size )
            ru.resize (size, false);
        #pragma omp parallel for
        for (TIndexType i = 0; i < size; ++i)
            ru[i] = rtot[indices[i]];
    }

    //this function assign to a vector the values that corresponds to u-dofs
    template<typename TVectorType, typename TIndexType = std::size_t>
    static void WritePart (TVectorType& rtot, const std::vector<TIndexType>& indices, const TVectorType& ru)
    {
        const TIndexType size = indices.size();
        #pragma omp parallel for
        for (TIndexType i = 0; i < size; ++i)
            rtot[indices[i]] = ru[i];
    }

    ///this function extracts the subblock A of matrix J
    ///as J = ( A  B1 ) u
    ///       ( B2 C  ) p
    template<typename TSparseMatrixType, typename TIndexType = std::size_t, typename TValueType = double>
    static void FillBlockMatrices (const TSparseMatrixType& rA,
        const std::vector<TIndexType>& first_indices,
        const std::vector<TIndexType>& second_indices,
        const std::vector<TIndexType>& global_to_local_indexing,
        const std::vector<int>& is_second_block,
        TSparseMatrixType& A)
    {
        KRATOS_TRY

        //get access to J data
        const TIndexType* index1 = rA.index1_data().begin();
        const TIndexType* index2 = rA.index2_data().begin();
        const TValueType* values = rA.value_data().begin();

        A.clear();

        //do allocation
        A.resize(first_indices.size(), first_indices.size(), false);

        //allocate the blocks by push_back
        for (TIndexType i = 0; i < rA.size1(); ++i)
        {
            TIndexType row_begin = index1[i];
            TIndexType row_end   = index1[i + 1];
            TIndexType local_row_id = global_to_local_indexing[i];

            if ( is_second_block[i] == false) //either A or B1
            {
                for (TIndexType j = row_begin; j < row_end; ++j)
                {
                    TIndexType col_index = index2[j];
                    TValueType value = values[j];
                    TIndexType local_col_id = global_to_local_indexing[col_index];
                    if (is_second_block[col_index] == false) //A block
//                        A.push_back ( local_row_id, local_col_id, value);
                        A( local_row_id, local_col_id ) = value;
                }
            }
        }

        A.complete_index1_data();

        KRATOS_CATCH ("")
    }

    ///this function generates the subblocks of matrix J
    ///as J = ( A  B1 ) u
    ///       ( B2 C  ) p
    template<typename TSparseMatrixType, typename TIndexType = std::size_t, typename TValueType = double>
    static void FillBlockMatrices (const TSparseMatrixType& rA,
        const std::vector<TIndexType>& first_indices,
        const std::vector<TIndexType>& second_indices,
        const std::vector<TIndexType>& global_to_local_indexing,
        const std::vector<int>& is_second_block,
        TSparseMatrixType& A, TSparseMatrixType& B1, TSparseMatrixType& B2, TSparseMatrixType& C)
    {
        KRATOS_TRY

        //get access to J data
        const TIndexType* index1 = rA.index1_data().begin();
        const TIndexType* index2 = rA.index2_data().begin();
        const TValueType* values = rA.value_data().begin();

        A.clear();
        B1.clear();
        B2.clear();
        C.clear();

        //do allocation
        A.resize(first_indices.size(), first_indices.size(), false);
        B1.resize(first_indices.size(), second_indices.size(), false);
        B2.resize(second_indices.size(), first_indices.size(), false);
        C.resize(second_indices.size(), second_indices.size(), false);

        //allocate the blocks by push_back
        for (TIndexType i = 0; i < rA.size1(); ++i)
        {
            TIndexType row_begin = index1[i];
            TIndexType row_end   = index1[i + 1];
            TIndexType local_row_id = global_to_local_indexing[i];

            if ( is_second_block[i] == false) //either A or B1
            {
                for (TIndexType j = row_begin; j < row_end; ++j)
                {
                    TIndexType col_index = index2[j];
                    TValueType value = values[j];
                    TIndexType local_col_id = global_to_local_indexing[col_index];
                    if (is_second_block[col_index] == false) //A block
//                        A.push_back ( local_row_id, local_col_id, value);
                        A( local_row_id, local_col_id ) = value;
                    else //B1 block
//                        B1.push_back ( local_row_id, local_col_id, value);
                        B1( local_row_id, local_col_id ) = value;
                }
            }
            else //either B2 or C
            {
                for (TIndexType j = row_begin; j < row_end; ++j)
                {
                    TIndexType col_index = index2[j];
                    TValueType value = values[j];
                    TIndexType local_col_id = global_to_local_indexing[col_index];
                    if (is_second_block[col_index] == false) //B2 block
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

    template<typename TSparseMatrixType, typename TIndexType = std::size_t, typename TValueType = double>
    static TValueType ComputeFrobeniusNorm(const TSparseMatrixType& rA)
    {
        const TIndexType n = rA.size1();
        const TIndexType* ia = rA.index1_data().begin();
        const TIndexType* ja = rA.index2_data().begin();
        const TValueType* a  = rA.value_data().begin();

        TValueType norm = 0.0;
        for(TIndexType i = 0; i < n; ++i)
        {
            int nz = ia[i + 1] - ia[i];
            for(int j = 0; j < nz; ++j)
                norm += pow(a[ia[i] + j], 2);
        }
        return sqrt(norm);
    }

    template<typename TVectorType, typename TIndexType = std::size_t>
    static void VectorScale(TVectorType& rX, const TVectorType& rD)
    {
        TIndexType n = rX.size();

        #pragma omp parallel for
        for(TIndexType i = 0; i < n; ++i)
            rX(i) *= rD(i);
    }

    template<typename TSparseMatrixType, typename TVectorType, typename TIndexType = std::size_t, typename TValueType = double>
    static void RowScale(TSparseMatrixType& rA, const TVectorType& rDL)
    {
        const TIndexType n = rA.size1();
        const TIndexType* ia = rA.index1_data().begin();
        const TIndexType* ja = rA.index2_data().begin();
        TValueType* a  = rA.value_data().begin();

        #pragma omp parallel for
        for(TIndexType i = 0; i < n; ++i)
        {
            int nz = ia[i + 1] - ia[i];
            for(int j = 0; j < nz; ++j)
                a[ia[i] + j] *= rDL(i);
        }
    }

    template<typename TSparseMatrixType, typename TDenseMatrixType, typename TIndexType = std::size_t, typename TValueType = double>
    static void MatrixMult(TDenseMatrixType& rC, const TSparseMatrixType& rA, const TDenseMatrixType& rB)
    {
        const TIndexType n = rA.size1();
        const TIndexType m = rB.size2();
        const TIndexType* ia = rA.index1_data().begin();
        const TIndexType* ja = rA.index2_data().begin();
        const TValueType* a  = rA.value_data().begin();

        #pragma omp parallel for
        for(TIndexType i = 0; i < n; ++i)
        {
            int nz = ia[i + 1] - ia[i];
            for(int k = 0; k < nz; ++k)
            {
                TValueType v = a[ia[i] + k];
                TIndexType col = ja[ia[i] + k];
                for(std::size_t j = 0; j < m; ++j)
                {
                    rC(i, j) += v*rB(col, j);
                }
            }
        }
    }

    template<typename TSparseMatrixType, typename TIndexType = std::size_t, typename TValueType = double>
    static void CheckDiagonalDominance(const TSparseMatrixType& rA)
    {
        const TIndexType n = rA.size1();
        const TIndexType* ia = rA.index1_data().begin();
        const TIndexType* ja = rA.index2_data().begin();
        const TValueType* a  = rA.value_data().begin();

        TValueType diag_norm, off_diag_norm;
        int nz;
        for(TIndexType i = 0; i < n; ++i)
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

    template<typename TMatrixType>
    static double DiagonalNorm(const TMatrixType& rA)
    {
        const std::size_t n = rA.size1();
        double diag = 0.0;
        for(std::size_t i = 0; i < n; ++i)
            diag += pow(rA(i, i), 2);
        diag = sqrt(diag);
        return diag;
    }

    template<typename TSparseMatrixType, typename TIndexType = std::size_t, typename TValueType = double>
    static void EstimateMinMaxEigenvaluesByGershgorin(const TSparseMatrixType& rA, double& lambda_min, double& lambda_max)
    {
        const TIndexType n = rA.size1();
        const TIndexType* ia = rA.index1_data().begin();
        const TIndexType* ja = rA.index2_data().begin();
        const TValueType* a  = rA.value_data().begin();

        int nz;
        double diag, off_diag_abs_sum;
        lambda_min = 1.0e99;
        lambda_max = -1.0e99;
        for(TIndexType i = 0; i < n; ++i)
        {
            nz = ia[i + 1] - ia[i];
            diag = 0.0;
            off_diag_abs_sum = 0.0;
            for(int j = 0; j < nz; ++j)
            {
                if(ja[ia[i] + j] == i)
                    diag = a[ia[i] + j];
                else
                    off_diag_abs_sum += fabs(a[ia[i] + j]);
            }

            if (lambda_min > diag - off_diag_abs_sum) lambda_min = diag - off_diag_abs_sum;
            if (lambda_max < diag + off_diag_abs_sum) lambda_max = diag + off_diag_abs_sum;
        }
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

    void Print(GeometryType::Pointer pGeometry) const
    {
        for(std::size_t i = 0; i < pGeometry->size(); ++i)
            std::cout << " (" << (*pGeometry)[i].X()
                     << ", " << (*pGeometry)[i].Y()
                     << ", " << (*pGeometry)[i].Z() << "),";
        std::cout << std::endl;
    }

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "Finite Cell Auxiliary Utility";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
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

    /// Assignment operator.
    MultithreadedSolversMathUtils& operator=(MultithreadedSolversMathUtils const& rOther);

    /// Copy constructor.
    MultithreadedSolversMathUtils(MultithreadedSolversMathUtils const& rOther);

    ///@}

}; // Class MultithreadedSolversMathUtils

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream, MultithreadedSolversMathUtils& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const MultithreadedSolversMathUtils& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MULTITHREADED_SOLVERS_MATH_UTILITY_H_INCLUDED  defined
