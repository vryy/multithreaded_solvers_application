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


#if !defined(KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BICGSTAB_BLOCK_PRESSURE_SOLVER_H_INCLUDED )
#define KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BICGSTAB_BLOCK_PRESSURE_SOLVER_H_INCLUDED


// System includes
#include <cmath>


// External includes


// Project includes
#include "includes/variables.h"
#include "linear_solvers/iterative_solver.h"

#ifdef MULTITHREADED_SOLVERS_APP_USE_FEAST
#include "custom_eigen_solvers/feast_solver.h"
#endif

#ifdef MULTITHREADED_SOLVERS_APP_USE_ARPACK
#include "custom_eigen_solvers/arpack_solver.h"
#endif

//#define CHECK_EIGENVALUES
#ifndef DBL_MIN
#define DBL_MIN 2.2250738585072014e-308
#define BicgstabBlockPressureSolver_define_DBL_MIN
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
class BicgstabBlockPressureSolver : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TModelPartType, TPreconditionerType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of  BicgstabBlockPressureSolver
    KRATOS_CLASS_POINTER_DEFINITION( BicgstabBlockPressureSolver);

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
    BicgstabBlockPressureSolver() {}

    BicgstabBlockPressureSolver(
        ValueType NewTolerance,
        unsigned int NewMaxIterationsNumber,
        typename TPreconditionerType::Pointer pNewPreconditioner,
        DataType a1,
        DataType a2,
        DataType b1,
        DataType b2
    ) : BaseType(NewTolerance, NewMaxIterationsNumber, pNewPreconditioner)
    {
        ma1 = a1;
        ma2 = a2;
        mb1 = b1;
        mb2 = b2;
    }

    /// Copy constructor.
     BicgstabBlockPressureSolver(const  BicgstabBlockPressureSolver& Other) : BaseType(Other) {}

    /// Destructor.
    ~BicgstabBlockPressureSolver() override {}

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    BicgstabBlockPressureSolver& operator=(const  BicgstabBlockPressureSolver& Other)
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
        //count pressure dofs
        unsigned int n_pressure_dofs = 0;
        unsigned int tot_active_dofs = 0;
        unsigned int system_size = TSparseSpaceType::Size(rB);
        for (auto it = rdof_set.begin(); it != rdof_set.end(); ++it)
            if (it->EquationId() < system_size)
            {
                ++tot_active_dofs;
                if ( (it->GetVariable().Key() == WATER_PRESSURE) || (it->GetVariable().Key() == AIR_PRESSURE) )
                    ++n_pressure_dofs;
            }
        if (tot_active_dofs != rA.size1() )
            KRATOS_ERROR << "total system size does not coincide with the free dof map";

        // KRATOS_WATCH(tot_active_dofs)
        // KRATOS_WATCH(n_pressure_dofs)

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
        for (auto it = rdof_set.begin(); it != rdof_set.end(); ++it)
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

        if(BaseType::GetPreconditioner()->AdditionalPhysicalDataIsNeeded())
            BaseType::GetPreconditioner()->ProvideAdditionalData(rA, rX, rB, rdof_set, r_model_part);
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

        //GetTimeTable()->Start(Info());

        std::cout << "Scaling begin" << std::endl;
        double start = OpenMPUtils::GetCurrentTime();
        RowScale(rA);
        std::cout << "Row scale completed..." << OpenMPUtils::GetCurrentTime() - start << " s" << std::endl;
        double start1 = OpenMPUtils::GetCurrentTime();
        ColumnScale(rA);
        std::cout << "Column scale completed..." << OpenMPUtils::GetCurrentTime() - start1 << " s" << std::endl;

        for(int i = 0; i < mother_indices.size(); ++i)
            rB(mother_indices[i]) *= ma1;

        for(int i = 0; i < mpressure_indices.size(); ++i)
            rB(mpressure_indices[i]) *= ma2;

        std::cout << "Scaling completed... " << OpenMPUtils::GetCurrentTime() - start << " s" << std::endl;

        #ifdef CHECK_EIGENVALUES
        EigenSolverType eigen_solver;
        VectorType lambdas(10);
//        std::cout << "Checking largest eigenvalues..." << std::endl;
//        eigen_solver.SolveLargest(rA, 5, lambdas);
//        std::cout << "Checking smallest eigenvalues..." << std::endl;
//        eigen_solver.SolveSmallest(rA, 5, lambdas);

        std::cout << "Checking eigenvalues between -1.0 and 1.0 ..." << std::endl;
        int nlambda0 = 50;
        VectorType rLamdas;
        eigen_solver.Solve(rA, nlambda0, rLamdas, -1.0, 1.0);
        #endif


        BaseType::GetPreconditioner()->Initialize(rA, rX, rB);

//        BaseType::GetPreconditioner()->ApplyInverseRight(rX);

//        BaseType::GetPreconditioner()->ApplyLeft(rB);

        int ierr = IterativeSolve(rA, rX, rB);

        if(ierr != 0)
            std::cout << "Warning: the iterative solver encountered some problem, error code = " << ierr << std::endl;

        BaseType::GetPreconditioner()->Finalize(rX);

        for(int i = 0; i < mother_indices.size(); ++i)
            rX(mother_indices[i]) *= mb1;

        for(int i = 0; i < mpressure_indices.size(); ++i)
            rX(mpressure_indices[i]) *= mb2;

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
        buffer << "Biconjugate-gradient stabilized iterative solver with pressure scaling strategy, preconditioner = " << BaseType::GetPreconditioner()->Info();
        return  buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& OStream) const override
    {
        OStream << "Biconjugate-gradient stabilized iterative solver with pressure scaling strategy, preconditioner = ";
        BaseType::GetPreconditioner()->PrintInfo(OStream);
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

    // scaling factors
    DataType ma1, ma2, mb1, mb2;

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

    void RowScale(SparseMatrixType& rA)
    {
        int n = rA.size1();
        std::size_t*   ia = rA.index1_data().begin();
        std::size_t*   ja = rA.index2_data().begin();
        double*         a = rA.value_data().begin();

        #pragma omp parallel for
        for(int i = 0; i < mother_indices.size(); ++i)
        {
            int r = mother_indices[i];
            int nz = ia[r + 1] - ia[r];
            for(int j = 0; j < nz; ++j)
                a[ia[r] + j] *= ma1;
        }

        #pragma omp parallel for
        for(int i = 0; i < mpressure_indices.size(); ++i)
        {
            int r = mpressure_indices[i];
            int nz = ia[r + 1] - ia[r];
            for(int j = 0; j < nz; ++j)
                a[ia[r] + j] *= ma2;
        }
    }

    void ColumnScale(SparseMatrixType& rA)
    {
        int n = rA.size1();
        std::size_t*   ia = rA.index1_data().begin();
        std::size_t*   ja = rA.index2_data().begin();
        double*         a = rA.value_data().begin();

        std::vector<std::vector<int> > all_indices_column_pos(n);
        for(int i = 0; i < n; ++i)
        {
            int nz = ia[i + 1] - ia[i];
            for(int j = 0; j < nz; ++j)
                all_indices_column_pos[ja[ia[i] + j]].push_back(ia[i] + j);
        }

        #pragma omp parallel for
        for(int i = 0; i < mother_indices.size(); ++i)
        {
            int m = all_indices_column_pos[mother_indices[i]].size();
            for(int j = 0; j < m; ++j)
                a[all_indices_column_pos[mother_indices[i]][j]] *= mb1;
        }

        #pragma omp parallel for
        for(int i = 0; i < mpressure_indices.size(); ++i)
        {
            int m = all_indices_column_pos[mpressure_indices[i]].size();
            for(int j = 0; j < m; ++j)
                a[all_indices_column_pos[mpressure_indices[i]][j]] *= mb2;
        }
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

}; // Class  BicgstabBlockPressureSolver

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.

#undef CHECK_EIGENVALUES

#undef BicgstabBlockPressureSolver_define_DBL_MIN
#undef DBL_MIN

#endif //  KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BICGSTAB_SOLVER_H_INCLUDED  defined
