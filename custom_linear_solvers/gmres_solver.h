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


#if !defined(KRATOS_MULTITHREADED_SOLVERS_APPLICATION_GMRES_SOLVER_H_INCLUDED )
#define  KRATOS_MULTITHREADED_SOLVERS_APPLICATION_GMRES_SOLVER_H_INCLUDED



// System includes
#include <cmath>


// External includes


// Project includes
#include "includes/define.h"
#include "utilities/progress.h"
#include "linear_solvers/iterative_solver.h"


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
class GmresSolver : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TModelPartType, TPreconditionerType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of GmresSolver
    KRATOS_CLASS_POINTER_DEFINITION(GmresSolver);

    typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TModelPartType, TPreconditionerType, TReordererType> BaseType;

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
    GmresSolver() {}

    GmresSolver(ValueType NewTolerance, unsigned int NewMaxIterationsNumber, unsigned int NewRestart, typename TPreconditionerType::Pointer pNewPreconditioner)
    : BaseType(NewTolerance, NewMaxIterationsNumber, pNewPreconditioner)
    {
        mRestart = NewRestart;
    }

    /// Copy constructor.
    GmresSolver(const GmresSolver& Other) : BaseType(Other)
    {
        mRestart = Other.mRestart;
    }

    /// Destructor.
    ~GmresSolver() override {}

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    GmresSolver& operator=(const GmresSolver& Other)
    {
        BaseType::operator=(Other);
        mRestart = Other.mRestart;
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

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

        //GetTimeTable()->Start(Info());

        BaseType::GetPreconditioner()->Initialize(rA, rX, rB);

//        BaseType::GetPreconditioner()->ApplyInverseRight(rX);

//        BaseType::GetPreconditioner()->ApplyLeft(rB);

        bool is_solved = IterativeSolve(rA, rX, rB);

        BaseType::GetPreconditioner()->Finalize(rX);

        //GetTimeTable()->Stop(Info());

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

            is_solved &= IterativeSolve(rA, x, b);

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
        buffer << "Generalized minimal residual iterative solver with " << BaseType::GetPreconditioner()->Info();
        return  buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& OStream) const override
    {
        OStream << "Generalized minimal residual iterative solver with ";
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

    unsigned int mRestart;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void Update(VectorType& x, const int k, const DenseMatrixType h, const VectorType s, const std::vector<VectorType> v) {
        VectorType y(s);

        // Backsolve:
        for (int i = k; i >= 0; i--) {
            y(i) /= h(i, i);
            for (int j = i - 1; j >= 0; j--)
                y(j) -= h(j, i) * y(i);
        }

        for (int j = 0; j <= k; j++)
            TSparseSpaceType::UnaliasedAdd(x, y(j), v[j]);
    }

    inline void GeneratePlaneRotation(const double dx, const double dy, double& cs, double& sn) {
        if (dy == 0.0) {
            cs = 1.0;
            sn = 0.0;
        } else if (std::abs(dy) > std::abs(dx)) {
            double temp = dx / dy;
            sn = 1.0 / sqrt(1.0 + temp * temp);
            cs = temp * sn;
        } else {
            double temp = dy / dx;
            cs = 1.0 / sqrt(1.0 + temp * temp);
            sn = temp * cs;
        }
    }

    inline void ApplyPlaneRotation(double& dx, double& dy, const double cs, const double sn) {
        double temp = cs * dx + sn * dy;
        dy = -sn * dx + cs * dy;
        dx = temp;
    }

    bool IterativeSolve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
        ValueType resid, tol = this->GetTolerance();
        unsigned int max_iter = this->GetMaxIterationsNumber();
        int i, j = 1, k, m = mRestart, size = rX.size();
        VectorType s(m + 1), cs(m + 1), sn(m + 1), w(size);
        VectorType r(rB);
        DenseMatrixType H(m+1,m);

        this->GetPreconditioner()->ApplyLeft(r); //r=P^-1*x
        DataType normb = TSparseSpaceType::TwoNorm(r);
        TSparseSpaceType::Mult(rA, rX, r); //r=A*x
        TSparseSpaceType::ScaleAndAdd(1.0, rB, -1.0, r); //r=b-A*x
        this->GetPreconditioner()->ApplyLeft(r); //r=P^-1*(b-A*x)
        DataType beta = TSparseSpaceType::TwoNorm(r);

        if (normb == 0.0)
            normb = 1.0;

        if ((resid = beta / normb) <= tol) {
            tol = resid;
            max_iter = 0;
            this->SetIterationsNumber(max_iter);
            BaseType::mBNorm = normb;
            this->SetResidualNorm(resid);
            return true;
        }

        std::vector<VectorType> v;
        for(i = 0; i < m + 1; ++i)
        {
            VectorType tmpv(size);
            v.push_back(tmpv);
        }

        Kratos::progress_display show_progress( max_iter );
        while (j <= max_iter) {
            TSparseSpaceType::Assign(v[0], 1.0 / beta, r);

            TSparseSpaceType::Set(s, 0.0);
            s(0) = beta;

            for (i = 0; i < m && j <= max_iter; i++, j++) {
                TSparseSpaceType::Mult(rA, v[i], w);
                this->GetPreconditioner()->ApplyLeft(w);
                for (k = 0; k <= i; k++) {
                    H(k, i) = TSparseSpaceType::Dot(w, v[k]);
                    TSparseSpaceType::UnaliasedAdd(w, -H(k, i), v[k]);
                }
                H(i + 1, i) = TSparseSpaceType::TwoNorm(w);
                TSparseSpaceType::Assign(v[i + 1], 1.0 / H(i + 1, i), w);

                for (k = 0; k < i; k++)
                    ApplyPlaneRotation(H(k, i), H(k + 1, i), cs(k), sn(k));

                GeneratePlaneRotation(H(i, i), H(i + 1, i), cs(i), sn(i));
                ApplyPlaneRotation(H(i, i), H(i + 1, i), cs(i), sn(i));
                ApplyPlaneRotation(s(i), s(i + 1), cs(i), sn(i));

                if ((resid = std::abs(s(i + 1)) / normb) < tol) {
                    Update(rX, i, H, s, v);
                    tol = resid;
                    max_iter = j;
                    this->SetIterationsNumber(max_iter);
                    BaseType::mBNorm = normb;
                    this->SetResidualNorm(resid);
                    return true;
                }
            }

            Update(rX, m - 1, H, s, v);

            TSparseSpaceType::Mult(rA, rX, r);
            TSparseSpaceType::ScaleAndAdd(1.0, rB, -1.0, r);
            this->GetPreconditioner()->ApplyLeft(r);
            beta = TSparseSpaceType::TwoNorm(r);
            if ((resid = beta / normb) < tol) {
                tol = resid;
                max_iter = j;
                this->SetIterationsNumber(max_iter);
                BaseType::mBNorm = normb;
                this->SetResidualNorm(resid);
                return true;
            }
            ++show_progress;
        }

        tol = resid;
        this->SetIterationsNumber(max_iter);
        BaseType::mBNorm = normb;
        this->SetResidualNorm(resid);
        return false;
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

}; // Class GmresSolver

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

}  // namespace Kratos.

#endif // KRATOS_MULTITHREADED_SOLVERS_APPLICATION_GMRES_SOLVER_H_INCLUDED  defined
