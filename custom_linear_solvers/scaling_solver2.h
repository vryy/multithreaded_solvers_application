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


#if !defined(KRATOS_MULTITHREADED_SOLVERS_APPLICATION_SCALING_SOLVER_2_H_INCLUDED )
#define KRATOS_MULTITHREADED_SOLVERS_APPLICATION_SCALING_SOLVER_2_H_INCLUDED


// System includes
#include <cmath>


// External includes
#include "boost/smart_ptr.hpp"
#include "mkl_lapacke.h"


// Project includes
#include "includes/define.h"
#include "includes/matrix_market_interface.h"
#include "linear_solvers/iterative_solver.h"

#ifdef MULTITHREADED_SOLVERS_APP_USE_FEAST
#include "custom_utilities/feast_solver.h"
#endif

#ifdef MULTITHREADED_SOLVERS_APP_USE_ARPACK
#include "custom_utilities/arpack_solver.h"
#endif

#include "hsl.h"
#include "external_includes/condition.hpp"
#ifdef MULTITHREADED_SOLVERS_APPLICATION_USE_PARDISO
#include "custom_utilities/pardiso_condition.h"
#endif
// #define CHECK_EIGENVALUES
// #define EXPORT_MATRIX_BEFORE_SOLVE

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
         class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class ScalingSolver2 : public LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of  ScalingSolver2
    KRATOS_CLASS_POINTER_DEFINITION( ScalingSolver2);

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef std::size_t  SizeType;

    typedef std::size_t  IndexType;

    #if defined(MULTITHREADED_SOLVERS_APP_USE_FEAST) && defined(CHECK_EIGENVALUES_USING_FEAST)
    typedef FeastSolver<TSparseSpaceType, TDenseSpaceType> EigenSolverType;
    #endif

    #if defined(MULTITHREADED_SOLVERS_APP_USE_ARPACK) && defined(CHECK_EIGENVALUES_USING_ARPACK)
    typedef ArpackSolver<TSparseSpaceType, TDenseSpaceType> EigenSolverType;
    #endif

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ScalingSolver2(typename BaseType::Pointer pSolver)
    {
        mpSolver = pSolver;
        mScalingType = 0; //no scaling
        mCheckConditionNunber = false;
    }

    ScalingSolver2(typename BaseType::Pointer pSolver, int ScalingType)
    {
        mpSolver = pSolver;
        mScalingType = ScalingType; // 1: Ruiz scaling strategy
                                    // 2: MC29 scaling strategy
        mCheckConditionNunber = false;
    }

    /// Copy constructor.
    ScalingSolver2(const ScalingSolver2& Other) : BaseType(Other),
        mScalingType(Other.mScalingType),
        mCheckConditionNunber(Other.mCheckConditionNunber),
        mpSolver(Other.mpSolver)
    {}

    /// Destructor.
    virtual ~ScalingSolver2() {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    ScalingSolver2& operator=(const  ScalingSolver2& Other)
    {
        BaseType::operator=(Other);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    void EnableCheckConditionNumber()
    {
        mCheckConditionNunber = true;
    }

    void DisableCheckConditionNumber()
    {
        mCheckConditionNunber = false;
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
        int n = rA.size1();
        std::size_t* ia = rA.index1_data().begin();
        std::size_t* ja = rA.index2_data().begin();

        if(mall_indices_column_pos.size() != 0)
            for(unsigned int i = 0; i < mall_indices_column_pos.size(); ++i)
                mall_indices_column_pos[i].clear();
        mall_indices_column_pos.clear();
        mall_indices_column_pos.resize(n);
        for(int i = 0; i < n; ++i)
        {
            int nz = ia[i + 1] - ia[i];
            for(int j = 0; j < nz; ++j)
                mall_indices_column_pos[ja[ia[i] + j]].push_back(ia[i] + j);
        }

        if(mScalingType == 2 || mScalingType == 3)
        {
            // for MC29 & MC77
            int NZ = 0;
            for(int i = 0; i < n; ++i)
            {
                int nz = ia[i + 1] - ia[i];
                NZ += nz;
                for(int j = 0; j < nz; ++j)
                {
                    mcrow.push_back(i + 1);
                    mccol.push_back(ja[ia[i] + j] + 1);
                }
            }
            mcrow.resize(NZ);
            mccol.resize(NZ);
        }

        if(mpSolver->AdditionalPhysicalDataIsNeeded())
            mpSolver->ProvideAdditionalData(rA, rX, rB, rdof_set, r_model_part);

        KRATOS_WATCH(mall_indices_column_pos.size())
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

        SparseMatrixType Acopy = rA;
        VectorType       Bcopy = rB;

        //GetTimeTable()->Start(Info());
        double start;

        if(mCheckConditionNunber)
        {
            start = OpenMPUtils::GetCurrentTime();
            std::cout << "ComputeConditionNumber before scaling" << std::endl;

            #ifdef MULTITHREADED_SOLVERS_APPLICATION_USE_PARDISO
            PardisoCondition<SparseMatrixType, VectorType> CondSolver;
            double norm1_invA = CondSolver.hager_norm1_inv_A(Acopy);
            KRATOS_WATCH(norm1_invA)
            double norm1_A = 0.0;
            int n = Acopy.size1();
            double* a = Acopy.value_data().begin();
            // compute norm_1 of A
            for(int j = 0; j < n; ++j)
            {
                double col_norm = 0.0;
                int nz = mall_indices_column_pos[j].size();
                for(int i = 0; i < nz; ++i)
                    col_norm += fabs(a[mall_indices_column_pos[j][i]]);
                if(col_norm > norm1_A)
                    norm1_A = col_norm;
            }
            KRATOS_WATCH(norm1_A)
            double Cond = norm1_A * norm1_invA;
            #else
//            double Cond = ComputeConditionNumberDGECON(Acopy);
            double Cond = ComputeConditionNumberMC75(Acopy);
//            double Cond = ComputeConditionNumberHager(Acopy);
            #endif
            KRATOS_WATCH(Cond)
            std::cout << "ComputeConditionNumber completed..." << OpenMPUtils::GetCurrentTime() - start << " s" << std::endl;
        }

        std::cout << "CompuleScalingVector begin" << std::endl;
        start = OpenMPUtils::GetCurrentTime();
        VectorType DL, DR;
        switch(mScalingType)
        {
            case 1:
                CompuleScalingVector(Acopy, DL, DR);
                break;
            case 2:
                CompuleScalingVector_mc29(Acopy, DL, DR);
                break;
            default:
                DR.resize(rX.size());
                DL.resize(Bcopy.size());
//                KRATOS_WATCH(rX.size())
//                KRATOS_WATCH(Bcopy.size())
                std::fill(DL.begin(), DL.end(), 1.0);
                std::fill(DR.begin(), DR.end(), 1.0);
                break;
        }
        std::cout << "CompuleScalingVector completed..." << OpenMPUtils::GetCurrentTime() - start << " s" << std::endl;

        std::cout << "Scaling begin" << std::endl;
        start = OpenMPUtils::GetCurrentTime();
        RowScale(Acopy, DL);
        std::cout << "Row scale completed..." << OpenMPUtils::GetCurrentTime() - start << " s" << std::endl;
        double start1 = OpenMPUtils::GetCurrentTime();
        ColumnScale(Acopy, DR);
        std::cout << "Column scale completed..." << OpenMPUtils::GetCurrentTime() - start1 << " s" << std::endl;
        VectorScale(Bcopy, DL);
        std::cout << "Scaling completed... " << OpenMPUtils::GetCurrentTime() - start << " s" << std::endl;

        if(mCheckConditionNunber)
        {
            start = OpenMPUtils::GetCurrentTime();
            std::cout << "ComputeConditionNumber after scaling" << std::endl;
            #ifdef MULTITHREADED_SOLVERS_APPLICATION_USE_PARDISO
            PardisoCondition<SparseMatrixType, VectorType> CondSolver;
            double norm1_invA = CondSolver.hager_norm1_inv_A(Acopy);
            KRATOS_WATCH(norm1_invA)
            // compute norm_1 of A after scaling
            double norm1_A = 0.0;
            int n = Acopy.size1();
            double* a = Acopy.value_data().begin();
            for(int j = 0; j < n; ++j)
            {
                double col_norm = 0.0;
                int nz = mall_indices_column_pos[j].size();
                for(int i = 0; i < nz; ++i)
                    col_norm += fabs(a[mall_indices_column_pos[j][i]]);
                if(col_norm > norm1_A)
                    norm1_A = col_norm;
            }
            KRATOS_WATCH(norm1_A)
            double Cond = norm1_A * norm1_invA;
            #else
//            double Cond = ComputeConditionNumberDGECON(Acopy);
            double Cond = ComputeConditionNumberMC75(Acopy);
//            double Cond = ComputeConditionNumberHager(Acopy);
            #endif
            KRATOS_WATCH(Cond)
            std::cout << "ComputeConditionNumber completed..." << OpenMPUtils::GetCurrentTime() - start << " s" << std::endl;
        }

        #ifdef CHECK_EIGENVALUES_USING_ARPACK
//        EigenSolverType eigen_solver;
        EigenSolverType eigen_solver(this->GetPreconditioner()); //using preconditioner for inverse iteration
        VectorType lambdas(10);
        std::cout << "Checking largest eigenvalues..." << std::endl;
        eigen_solver.SolveLargest(Acopy, 5, lambdas);
        std::cout << "Checking smallest eigenvalues..." << std::endl;
        eigen_solver.SolveSmallest(Acopy, 5, lambdas);

//        std::cout << "Checking eigenvalues between -1.0 and 1.0 ..." << std::endl;
//        int nlambda0 = 50;
//        VectorType rLamdas;
//        eigen_solver.Solve(Acopy, nlambda0, rLamdas, -1.0, 1.0);
        #endif

        #ifdef EXPORT_MATRIX_BEFORE_SOLVE
        std::stringstream lhs_filename;
        lhs_filename << "A_scaled.mm";
        WriteMatrixMarketMatrix(lhs_filename.str().c_str(), Acopy, false);
        #endif

        bool is_solved = mpSolver->Solve(Acopy, rX, Bcopy);

        VectorScale(rX, DR);

        return is_solved;
    }

    /** Multi solve method for solving a set of linear systems with same coefficient matrix.
    Solves the linear system Ax=b and puts the result on SystemVector& rX.
    rX is also th initial guess for iterative methods.
    @param rA. System matrix
    @param rX. Solution vector. it's also the initial
    guess for iterative linear solvers.
    @param rB. Right hand side vector.
    */
    bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Multisolve is not yet implemented for", typeid(*this).name())
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
        buffer << "Scaling solver with ";
        switch(mScalingType)
        {
            case 1:
                buffer << "Ruiz scaling strategy";
                break;
            case 2:
                buffer << "MC29 scaling strategy";
                break;
            default:
                buffer << "no scaling";
                break;
        }
        buffer << ", wrapped solver: " << mpSolver->Info();
        return  buffer.str();
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

    std::vector<std::vector<int> > mall_indices_column_pos; //2d vector to store indices of entry in value vector for each colunn

    // for MC29 & MC77
    std::vector<int> mcrow; //row coordinates
    std::vector<int> mccol; //column coordinates

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
    int mScalingType;
    bool mCheckConditionNunber;
    typename BaseType::Pointer mpSolver;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

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

    void ColumnScale(SparseMatrixType& rA, VectorType& rDR)
    {
        int n = rA.size1();
        std::size_t*   ia = rA.index1_data().begin();
        std::size_t*   ja = rA.index2_data().begin();
        double*         a = rA.value_data().begin();

        #pragma omp parallel for
        for(int i = 0; i < n; ++i)
        {
            int nz = mall_indices_column_pos[i].size();
            for(int j = 0; j < nz; ++j)
                a[mall_indices_column_pos[i][j]] *= rDR(i);
        }
    }

    void VectorScale(VectorType& rX, VectorType& rD)
    {
        int n = rX.size();

        #pragma omp parallel for
        for(int i = 0; i < n; ++i)
            rX(i) *= rD(i);
    }

    // scaling vector computation using Ruiz algorithm
    void CompuleScalingVector(SparseMatrixType& rA, VectorType& rDL, VectorType& rDR)
    {
        // Create a clone of matrix A
        SparseMatrixType rAc = rA;

        int n = rAc.size1();
        std::size_t*   ia = rAc.index1_data().begin();
        std::size_t*   ja = rAc.index2_data().begin();
        double*         a = rAc.value_data().begin();

        rDL.resize(n);
        rDR.resize(n);

        // firstly initialize all scaling matrix to identity
        for(unsigned int i = 0; i < n; ++i)
        {
            rDL(i) = 1.0;
            rDR(i) = 1.0;
        }

        // iterations
        VectorType DL(n);
        VectorType DR(n);
        bool converged = false;
        int nz;
        double v, norm, conv_r, conv_c;
        unsigned int i, j, cnt = 0;
        double tol = 1.0e-10;
        do
        {
            conv_r = 0.0;
            conv_c = 0.0;
            for(i = 0; i < n; ++i)
            {
                // find the inf-norm of row i
                nz = ia[i + 1] - ia[i];
                norm = 0.0;
                for(int j = 0; j < nz; ++j)
                {
                    v = sqrt(fabs(a[ia[i] + j]));
                    if(v > norm) norm = v;
                }
                DL(i) = 1.0 / norm;
                v = fabs(1 - pow(norm, 2));
                if(v > conv_r) conv_r = v;

                // find the inf-norm of column i
                nz = mall_indices_column_pos[i].size();
                norm = 0.0;
                for(j = 0; j < nz; ++j)
                {
                    v = sqrt(fabs(a[mall_indices_column_pos[i][j]]));
                    if(v > norm) norm = v;
                }
                DR(i) = 1.0 / norm;
                v = fabs(1 - pow(norm, 2));
                if(v > conv_c) conv_c = v;
            }

            // Scale the matrix
            RowScale(rAc, DL);
            ColumnScale(rAc, DR);

            // Articulate to the scaling vector
            VectorScale(rDL, DL);
            VectorScale(rDR, DR);

            // Check the convergence criteria
            converged = (conv_r < tol) && (conv_c < tol);
            std::cout << "iteration " << ++cnt << ", conv_r = " << conv_r << ", conv_c = " << conv_c << std::endl;
        }
        while(converged != true);
    }

    // scaling vector computation using MC29
    int CompuleScalingVector_mc29(SparseMatrixType& rA, VectorType& rDL, VectorType& rDR)
    {
        int n = rA.size1();
        int ne = mcrow.size();
//        std::size_t*   ia = rA.index1_data().begin();
//        std::size_t*   ja = rA.index2_data().begin();
        double*         a = rA.value_data().begin();

        double* r = new double[n];
        double* c = new double[n];
        double* w = new double[5 * n];
        int lp = 0, ifail;
        mc29ad(&n, &n, &ne, a, &mcrow[0], &mccol[0], r, c, w, &lp, &ifail);

        rDL.resize(n);
        rDR.resize(n);
        for(unsigned int i = 0; i < n; ++i)
        {
            rDL(i) = exp(r[i]);
            rDR(i) = exp(c[i]);
        }

        delete [] r;
        delete [] c;
        delete [] w;

        return ifail;
    }

    /*
     * A routine to compute condition number of a matrix using svd decomposition (MKL required). It is very expensive so usage with big matrix is warning.
     */
    double ComputeConditionNumberDGESVD(SparseMatrixType& rA)
    {
        int n = rA.size1();
        MKL_INT lda = n, ldu = n, ldvt = n, info;
        double* superb = new double[n - 1];
        /* Local arrays */
        double* s = new double[n];
        double* u = new double[ldu * n];
        double* vt = new double[ldvt * n];
        double* a = new double[lda * n];

        /* Populate a */
        std::size_t*   ia = rA.index1_data().begin();
        std::size_t*   ja = rA.index2_data().begin();
        double*        va = rA.value_data().begin();
        for(int i = 0; i < n; ++i)
        {
            for(int j = 0; j < n; ++j)
                a[i * n + j] = 0.0;
            int nz = ia[i + 1] - ia[i];
            for(int j = 0; j < nz; ++j)
                a[i * n + ja[ia[i] + j]] = va[ia[i] + j];
        }

        /* Compute SVD */
        info = LAPACKE_dgesvd( LAPACK_ROW_MAJOR, 'A', 'A', n, n, a, lda, s, u, ldu, vt, ldvt, superb );

        if( info != 0 )
        {
            std::cout << "The algorithm computing SVD failed to converge." << std::endl;
            return 0.0;
        }

        double cond = s[0] / s[n-1];

        delete [] s;
        delete [] u;
        delete [] vt;
        delete [] a;

        return cond;
    }

    /*
     * A routine to compute condition number of a matrix using svd decomposition (MKL required). It is very expensive so usage with big matrix is warning.
     */
    double ComputeConditionNumberDGECON(SparseMatrixType& rA)
    {
        int n = rA.size1();
        MKL_INT lda = n, info;
        MKL_INT* ipiv = new MKL_INT[n];
        /* Local arrays */
        double* a = new double[lda * n];

        /* Populate a */
        std::size_t*   ia = rA.index1_data().begin();
        std::size_t*   ja = rA.index2_data().begin();
        double*        va = rA.value_data().begin();
        for(int i = 0; i < n; ++i)
        {
            for(int j = 0; j < n; ++j)
                a[i * n + j] = 0.0;
            int nz = ia[i + 1] - ia[i];
            for(int j = 0; j < nz; ++j)
                a[i * n + ja[ia[i] + j]] = va[ia[i] + j];
        }

        /* Compute 1-norm of A */
        double anorm = 0.0;
        for(int i = 0; i < n; ++i)
        {
            double col_norm = 0.0;
            for(int j = 0; j < n; ++j)
                col_norm += fabs(a[j*n + i]);
            if(col_norm > anorm)
                anorm = col_norm;
        }

        /* Compute LU decomposition of A */
        info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, a, lda, ipiv);
        if( info != 0 )
        {
            std::cout << "LU decomposition for A failed." << std::endl;
            return 0.0;
        }

        /* Compute estimated condition number */
        double rcond;
        info = LAPACKE_dgecon( LAPACK_ROW_MAJOR, 'I', n, a, lda, anorm, &rcond );
        if( info != 0 )
        {
            std::cout << "The algorithm computing condition number failed to converge." << std::endl;
            return 0.0;
        }

        delete [] ipiv;
        delete [] a;

        return 1.0 / rcond;
    }

    /*
     * A routine to compute condition number of a sparse matrix using mc75
     */
    double ComputeConditionNumberMC75(SparseMatrixType& rA)
    {
        int n = rA.size1();
        int ne = mcrow.size();
//        std::size_t*   ia = rA.index1_data().begin();
//        std::size_t*   ja = rA.index2_data().begin();
        double*         a = rA.value_data().begin();

        int la = (unsigned int)(ne * log2(ne));
        double* a_ = new double[la];
        int* ir = new int[la];
        int* ic = new int[la];
        double* cond = new double[2];
        int liw = 19 * n + 7;
        int* iw = new int[liw];
        int lw = 8 * n;
        double* w = new double[lw];
        int* icntl = new int[5];
        int* info = new int[5];

        std::copy(&mcrow[0], &mcrow[0] + ne, ir);
        std::copy(&mccol[0], &mccol[0] + ne, ic);
        std::copy(a, a + ne, a_);

        icntl[0] = 6;
        icntl[1] = 0;
        icntl[2] = 0;
        icntl[3] = 0;
        icntl[4] = 0;
        mc75ad(&n, &ne, &la, a_, ir, ic, cond, &liw, iw, &lw, w, icntl, info);

        double my_cond = cond[1];

        delete [] a;
        delete [] ir;
        delete [] ic;
        delete [] cond;
        delete [] iw;
        delete [] w;
        delete [] icntl;
        delete [] info;

        return my_cond;
    }

    /*
     * A routine to compute condition number of a matrix using Hager algorithm. Ref: W. W. Hager, Condition Estimates, SIAM
     */
    double ComputeConditionNumberHager(SparseMatrixType& rA)
    {
        int n = rA.size1();
        /* Local arrays */
        double* a = new double[n * n];

        /* Populate a */
        std::size_t*   ia = rA.index1_data().begin();
        std::size_t*   ja = rA.index2_data().begin();
        double*        va = rA.value_data().begin();
        for(int i = 0; i < n; ++i)
        {
            // remark: storage is by column
            for(int j = 0; j < n; ++j)
                a[j * n + i] = 0.0;
            int nz = ia[i + 1] - ia[i];
            for(int j = 0; j < nz; ++j)
                a[ja[ia[i] + j] * n + i] = va[ia[i] + j];
        }

        double cond = condition_hager(n, a);

        delete [] a;

        return cond;
    }

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

}; // Class  ScalingSolver2

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::istream& operator >> (std::istream& IStream, ScalingSolver2<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis)
{
    return IStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator << (std::ostream& OStream, const  ScalingSolver2<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis)
{
    rThis.PrintInfo(OStream);
    OStream << std::endl;
    rThis.PrintData(OStream);

    return OStream;
}
///@}


}  // namespace Kratos.

#undef CHECK_EIGENVALUES
#undef EXPORT_MATRIX_BEFORE_SOLVE

#endif //  KRATOS_MULTITHREADED_SOLVERS_APPLICATION_SCALING_SOLVER_2_H_INCLUDED  defined
