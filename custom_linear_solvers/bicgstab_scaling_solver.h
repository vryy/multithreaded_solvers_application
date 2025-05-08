/*
see multithreaded_solvers_application/LICENSE.txt
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 26 Aug 2014 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BICGSTAB_SCALING_SOLVER_H_INCLUDED )
#define KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BICGSTAB_SCALING_SOLVER_H_INCLUDED


// System includes
#include <cmath>

// External includes
#include "mkl_lapacke.h"

#ifdef MULTITHREADED_SOLVERS_APP_USE_HSL
#include "hsl.h"
#endif

// Project includes
#include "includes/define.h"
#include "includes/matrix_market_interface.h"
#include "linear_solvers/iterative_solver.h"

#ifdef MULTITHREADED_SOLVERS_APP_USE_FEAST
#include "custom_eigen_solvers/feast_solver.h"
#endif

#ifdef MULTITHREADED_SOLVERS_APP_USE_ARPACK
#include "custom_eigen_solvers/arpack_solver.h"
#endif

#include "external_includes/condition.hpp"
#include "custom_utilities/pardiso_condition.h"

// #define CHECK_EIGENVALUES
// #define EXPORT_MATRIX_BEFORE_SOLVE
#ifndef DBL_MIN
#define DBL_MIN 2.2250738585072014e-308
#define BicgstabScalingSolver_define_DBL_MIN
#endif

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
class BicgstabScalingSolver : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TModelPartType, TPreconditionerType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of  BicgstabScalingSolver
    KRATOS_CLASS_POINTER_DEFINITION( BicgstabScalingSolver);

    typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TModelPartType, TPreconditionerType, TReordererType> BaseType;

    typedef typename BaseType::SparseMatrixType SparseMatrixType;

    typedef typename BaseType::VectorType VectorType;

    typedef typename BaseType::DenseMatrixType DenseMatrixType;

    typedef typename BaseType::SizeType SizeType;

    typedef typename BaseType::IndexType IndexType;

    typedef typename BaseType::DataType DataType;

    typedef typename BaseType::ValueType ValueType;

    #ifdef CHECK_EIGENVALUES

    #if defined(MULTITHREADED_SOLVERS_APP_USE_FEAST) && defined(CHECK_EIGENVALUES_USING_FEAST)
    typedef FeastSolver<TSparseSpaceType, TDenseSpaceType> EigenSolverType;
    #endif

    #if defined(MULTITHREADED_SOLVERS_APP_USE_ARPACK) && defined(CHECK_EIGENVALUES_USING_ARPACK)
    typedef ArpackSolver<TSparseSpaceType, TDenseSpaceType> EigenSolverType;
    #endif

    #endif

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BicgstabScalingSolver()
    {
        mCheckConditionNunber = false;
    }

    BicgstabScalingSolver(
        double NewTolerance,
        unsigned int NewMaxIterationsNumber,
        typename TPreconditionerType::Pointer pNewPreconditioner
    ) : BaseType(NewTolerance, NewMaxIterationsNumber, pNewPreconditioner)
    {
        mScalingType = 0; //no scaling
        mCheckConditionNunber = false;
    }

    BicgstabScalingSolver(
        double NewTolerance,
        unsigned int NewMaxIterationsNumber,
        typename TPreconditionerType::Pointer pNewPreconditioner,
        int ScalingType
    ) : BaseType(NewTolerance, NewMaxIterationsNumber, pNewPreconditioner)
    {
        mScalingType = ScalingType; // 1: Ruiz scaling strategy
                                    // 2: MC29 scaling strategy
        mCheckConditionNunber = false;
    }

    /// Copy constructor.
     BicgstabScalingSolver(const  BicgstabScalingSolver& Other) : BaseType(Other), mCheckConditionNunber(false) {}

    /// Destructor.
    ~BicgstabScalingSolver() override {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    BicgstabScalingSolver& operator=(const  BicgstabScalingSolver& Other)
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
        int n = rA.size1();
        IndexType* ia = rA.index1_data().begin();
        IndexType* ja = rA.index2_data().begin();

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

        if(BaseType::GetPreconditioner()->AdditionalPhysicalDataIsNeeded())
            BaseType::GetPreconditioner()->ProvideAdditionalData(rA, rX, rB, rdof_set, r_model_part);

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
    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
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
            double Cond = 0.0;
//            Cond = ComputeConditionNumberDGECON(Acopy);
            #ifdef MULTITHREADED_SOLVERS_APP_USE_HSL
            Cond = ComputeConditionNumberMC75(Acopy);
            #endif
//            Cond = ComputeConditionNumberHager(Acopy);
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
            #ifdef MULTITHREADED_SOLVERS_APP_USE_HSL
            case 2:
                CompuleScalingVector_mc29(Acopy, DL, DR);
                break;
            #endif
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
            double Cond = 0.0;
//            Cond = ComputeConditionNumberDGECON(Acopy);
            #ifdef MULTITHREADED_SOLVERS_APP_USE_HSL
            Cond = ComputeConditionNumberMC75(Acopy);
            #endif
//            Cond = ComputeConditionNumberHager(Acopy);
            #endif
            KRATOS_WATCH(Cond)
            std::cout << "ComputeConditionNumber completed..." << OpenMPUtils::GetCurrentTime() - start << " s" << std::endl;
        }

        #ifdef CHECK_EIGENVALUES
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

        BaseType::GetPreconditioner()->Initialize(Acopy, rX, Bcopy);

//        BaseType::GetPreconditioner()->ApplyInverseRight(rX);

//        BaseType::GetPreconditioner()->ApplyLeft(Bcopy);

        int ierr = IterativeSolve(Acopy, rX, Bcopy);

        if(ierr != 0)
            std::cout << "Warning: the iterative solver encountered some problem, error code = " << ierr << std::endl;

        BaseType::GetPreconditioner()->Finalize(rX);

        VectorScale(rX, DR);

        //GetTimeTable()->Stop(Info());

        return (ierr == 0);
    }

    /** Multi solve method for solving a set of linear systems with same coefficient matrix.
    Solves the linear system Ax=b and puts the result on SystemVector& rX.
    rX is also th initial guess for iterative methods.
    @param rA. System matrix
    @param rX. Solution vector. it's also the initial
    guess for iterative linear solvers.
    @param rB. Right hand side vector.
    */
    bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB) override
    {
        //GetTimeTable()->Start(Info());

        BaseType::GetPreconditioner()->Initialize(rA, rX, rB);

        bool is_solved = true;
        VectorType x(TDenseSpaceType::Size1(rX));
        VectorType b(TDenseSpaceType::Size1(rB));
        for(unsigned int i = 0 ; i < TDenseSpaceType::Size2(rX) ; ++i)
        {
            TDenseSpaceType::GetColumn(i, rX, x);
            TDenseSpaceType::GetColumn(i, rB, b);

//            BaseType::GetPreconditioner()->ApplyInverseRight(x);
//            BaseType::GetPreconditioner()->ApplyLeft(b);

            int ierr = IterativeSolve(rA, x, b);
            if(ierr != 0)
                std::cout << "Warning: the iterative solver encountered some problem at solve #" << i << ", error code = " << ierr << std::endl;
            is_solved &= (ierr == 0);

            BaseType::GetPreconditioner()->Finalize(x);
        }

        //GetTimeTable()->Stop(Info());

        return is_solved;
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
        buffer << "Biconjugate-gradient stabilized iterative solver with ";
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
        buffer << ", preconditioner = " << BaseType::GetPreconditioner()->Info();
        return  buffer.str();
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

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    int IterativeSolve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
        std::cout.precision(15);
        ValueType resid, tol = this->GetTolerance();
        unsigned int max_iter = this->GetMaxIterationsNumber();
        int i, j = 1, k, size = rX.size();
        DataType rho_1, rho_2, alpha, beta, omega, norms, normr, normb;
        VectorType p(size), phat(size), s(size), shat(size), t(size), v(size), r(size), rtilde(size);

        normb = TSparseSpaceType::TwoNorm(rB);
        TSparseSpaceType::Mult(rA, rX, r); //r=A*x
        TSparseSpaceType::ScaleAndAdd(1.0, rB, -1.0, r); //r=b-A*x
        TSparseSpaceType::Assign(rtilde, 1.0, r);

        if (fabs(normb) < DBL_MIN)
            normb = 1.0;

        normr = TSparseSpaceType::TwoNorm(r);
        if ((resid = normr / normb) < tol) {
            tol = resid;
            max_iter = 0;
            this->SetIterationsNumber(max_iter);
            BaseType::mBNorm = normb;
            this->SetResidualNorm(resid);
            return 0;
        }

        Kratos::progress_display show_progress( max_iter );
        for (int i = 1; i <= max_iter; ++i) {
            rho_1 = TSparseSpaceType::Dot(rtilde, r);
//          KRATOS_WATCH(normr)
//          KRATOS_WATCH(rho_1)
//          KRATOS_WATCH(rho_2)
            if(fabs(rho_1) < DBL_MIN) {
                KRATOS_WATCH(TSparseSpaceType::TwoNorm(rtilde))
                KRATOS_WATCH(TSparseSpaceType::TwoNorm(r))
                tol = normr / normb;
                max_iter = i;
                this->SetIterationsNumber(max_iter);
                BaseType::mBNorm = normb;
                this->SetResidualNorm(resid);
                return 2;
            }
            if(i == 1)
                TSparseSpaceType::Assign(p, 1.0, r);
            else {
                beta = (rho_1 / rho_2) * (alpha / omega);
                TSparseSpaceType::UnaliasedAdd(p, -omega, v); //p = p - omega * v
                TSparseSpaceType::ScaleAndAdd(1.0, r, beta, p); // p = r + beta * p
            }
            TSparseSpaceType::Assign(phat, 1.0, p);
            this->GetPreconditioner()->ApplyLeft(phat);
            TSparseSpaceType::Mult(rA, phat, v);
            alpha = rho_1 / TSparseSpaceType::Dot(rtilde, v);
            TSparseSpaceType::ScaleAndAdd(1.0, r, -alpha, v, s);
            norms = TSparseSpaceType::TwoNorm(s);
//          KRATOS_WATCH(norms)
            if((resid = norms / normb) < tol) {
                TSparseSpaceType::UnaliasedAdd(rX, alpha, phat);
                tol = resid;
                max_iter = i;
                this->SetIterationsNumber(max_iter);
                BaseType::mBNorm = normb;
                this->SetResidualNorm(resid);
                return 0;
            }
            TSparseSpaceType::Assign(shat, 1.0, s);
//          KRATOS_WATCH(norm_2(shat))
            this->GetPreconditioner()->ApplyLeft(shat);
//          KRATOS_WATCH(norm_2(shat))
//          KRATOS_WATCH(norm_frobenius(rA))
            TSparseSpaceType::Mult(rA, shat, t);
//          KRATOS_WATCH(norm_2(t))
//          KRATOS_WATCH(norm_2(s))
            omega = TSparseSpaceType::Dot(t, s) / TSparseSpaceType::Dot(t, t);
            TSparseSpaceType::UnaliasedAdd(rX, alpha, phat);
            TSparseSpaceType::UnaliasedAdd(rX, omega, shat);
            TSparseSpaceType::ScaleAndAdd(1.0, s, -omega, t, r);

            rho_2 = rho_1;
            normr = TSparseSpaceType::TwoNorm(r);
//          KRATOS_WATCH(alpha)
//          KRATOS_WATCH(omega)
            if((resid = normr / normb) < tol) {
                tol = resid;
                max_iter = i;
                this->SetIterationsNumber(max_iter);
                BaseType::mBNorm = normb;
                this->SetResidualNorm(resid);
                return 0;
            }
            if(fabs(omega) < DBL_MIN) {
                tol = normr / normb;
                max_iter = i;
                this->SetIterationsNumber(max_iter);
                BaseType::mBNorm = normb;
                this->SetResidualNorm(resid);
                return 3;
            }

            ++show_progress;
        }

        tol = resid;
        this->SetIterationsNumber(max_iter);
        BaseType::mBNorm = normb;
        this->SetResidualNorm(resid);
        return 1;
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

    #ifdef MULTITHREADED_SOLVERS_APP_USE_HSL
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
    #endif

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

    #ifdef MULTITHREADED_SOLVERS_APP_USE_HSL
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
    #endif

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

}; // Class  BicgstabScalingSolver

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

}  // namespace Kratos.

#undef CHECK_EIGENVALUES
#undef EXPORT_MATRIX_BEFORE_SOLVE

#undef BicgstabScalingSolver_define_DBL_MIN
#undef DBL_MIN

#endif //  KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BICGSTAB_SOLVER_H_INCLUDED  defined
