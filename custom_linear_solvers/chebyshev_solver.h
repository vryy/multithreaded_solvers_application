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
//   Date:                $Date: 29 Sep 2018 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_MULTITHREADED_SOLVERS_APPLICATION_CHEBYSHEV_SOLVER_H_INCLUDED )
#define KRATOS_MULTITHREADED_SOLVERS_APPLICATION_CHEBYSHEV_SOLVER_H_INCLUDED



// System includes
#include <math.h>


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
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
         class TPreconditionerType = Preconditioner<TSparseSpaceType, TDenseSpaceType>,
         class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class ChebyshevSolver : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of ChebyshevSolver
    KRATOS_CLASS_POINTER_DEFINITION(ChebyshevSolver);

    typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ChebyshevSolver() : BaseType(), mMinEigenvalueMode("Gershgorin"), mMaxEigenvalueMode("Gershgorin")
    {}

    ChebyshevSolver(double NewTolerance, unsigned int NewMaxIterationsNumber)
    : BaseType(NewTolerance, NewMaxIterationsNumber), mMinEigenvalueMode("Gershgorin"), mMaxEigenvalueMode("Gershgorin")
    {}

    ChebyshevSolver(double NewTolerance, unsigned int NewMaxIterationsNumber, typename TPreconditionerType::Pointer pNewPreconditioner)
    : BaseType(NewTolerance, NewMaxIterationsNumber, pNewPreconditioner), mMinEigenvalueMode("Gershgorin"), mMaxEigenvalueMode("Gershgorin")
    {}

    ChebyshevSolver(double NewTolerance, unsigned int NewMaxIterationsNumber, typename TPreconditionerType::Pointer pNewPreconditioner,
        const double& lambda_min, const double& lambda_max)
    : BaseType(NewTolerance, NewMaxIterationsNumber, pNewPreconditioner)
    , mLambdaMin(lambda_min), mLambdaMax(lambda_max)
    , mMinEigenvalueMode("Provided"), mMaxEigenvalueMode("Provided")
    {}

    /// Copy constructor.
    ChebyshevSolver(const ChebyshevSolver& Other) 
    : BaseType(Other), mLambdaMin(Other.mLambdaMin), mLambdaMax(Other.mLambdaMax)
    , mMinEigenvalueMode(Other.mMinEigenvalueMode), mMaxEigenvalueMode(Other.mMaxEigenvalueMode)
    {}

    /// Destructor.
    virtual ~ChebyshevSolver() {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    ChebyshevSolver& operator=(const ChebyshevSolver& Other)
    {
        BaseType::operator=(Other);
        mLambdaMin = Other.mLambdaMin;
        mLambdaMax = Other.mLambdaMax;
        mMinEigenvalueMode = Other.mMinEigenvalueMode;
        mMaxEigenvalueMode = Other.mMaxEigenvalueMode;
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Set minimum eigenvalue estimation
    void SetMinEigenvalue(const double& lambda_min)
    {
        mLambdaMin = lambda_min;
        mMinEigenvalueMode = std::string("Provided");
    }

    /// Set maximum eigenvalue estimation
    void SetMaxEigenvalue(const double& lambda_max)
    {
        mLambdaMax = lambda_max;
        mMaxEigenvalueMode = std::string("Provided");
    }

    /// Set estimate minimum eigenvalue based on Gershgorin theorem. It is noted that the eigenvalue estimation can be negative.
    void SetEstimateMinEigenvalueByGershgorin()
    {
         mMinEigenvalueMode = std::string("Gershgorin");
    }

    /// Set estimate maximum eigenvalue based on Gershgorin theorem. It is noted that the eigenvalue estimation can be negative.
    void SetEstimateMaxEigenvalueByGershgorin()
    {
         mMaxEigenvalueMode = std::string("Gershgorin");
    }

    /** Normal solve method.
    Solves the linear system Ax=b and puts the result on SystemVector& rX.
    rX is also th initial guess for iterative methods.
    @param rA. System matrix
    @param rX. Solution vector. it's also the initial
    guess for iterative linear solvers.
    @param rB. Right hand side vector.
    REF: Fig 2.11: Templates book
    */
    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
        if(this->IsNotConsistent(rA, rX, rB))
            return false;

//        GetTimeTable()->Start(Info());

        if (this->GetPreconditioner() != NULL)
            BaseType::GetPreconditioner()->Initialize(rA, rX, rB);

        double tol = this->GetTolerance(), resid, normb, normr;
        unsigned int max_iter = this->GetMaxIterationsNumber();
        std::size_t size = rX.size();
        VectorType r(size), p(size);

        normb = TSparseSpaceType::TwoNorm(rB);
        if (normb == 0.0)
            normb = 1.0;

        // estimate the eigenvalues of needed
        if (mMinEigenvalueMode == std::string("Gershgorin") && mMaxEigenvalueMode == std::string("Gershgorin"))
        {
            MultithreadedSolversMathUtils::EstimateMinMaxEigenvaluesByGershgorin(rA, mLambdaMin, mLambdaMax);
        }
        else
        {
            if (mMinEigenvalueMode == std::string("Gershgorin") || mMaxEigenvalueMode == std::string("Gershgorin"))
            {
                if (mMinEigenvalueMode == std::string("Gershgorin"))
                {
                    double lambda_max;
                    MultithreadedSolversMathUtils::EstimateMinMaxEigenvaluesByGershgorin(rA, mLambdaMin, lambda_max);
                }

                if (mMaxEigenvalueMode == std::string("Gershgorin"))
                {
                    double lambda_min;
                    MultithreadedSolversMathUtils::EstimateMinMaxEigenvaluesByGershgorin(rA, lambda_min, mLambdaMax);
                }
            }
        }

        if (this->GetEchoLevel() > 1)
            std::cout << "ChebyshevSolver: Lambda min: " << mLambdaMin << ", max: " << mLambdaMax << std::endl;

        double c = (mLambdaMax - mLambdaMin) / 2;
        double d = (mLambdaMax + mLambdaMin) / 2;
        double alpha, beta;

        TSparseSpaceType::Mult(rA, rX, r); //r=A*x
        TSparseSpaceType::ScaleAndAdd(1.0, rB, -1.0, r); //r=b-A*x

        for (unsigned int i = 1; i <= max_iter; ++i)
        {
            normr = TSparseSpaceType::TwoNorm(r);

            if (this->GetEchoLevel() > 0)
                std::cout << "ChebyshevSolver iteration #" << i
                          << ", normr = " << normr << ", tol = " << tol
                          << std::endl;

            if ((resid = normr / normb) <= tol) {
                tol = resid;
                max_iter = i;
                this->SetIterationsNumber(max_iter);
                BaseType::mBNorm = normb;
                this->SetResidualNorm(resid);
                return true;
            }

            if (this->GetPreconditioner() != NULL)
                this->GetPreconditioner()->ApplyLeft(r); //r=P^(-1)*r (r=z)

            if (i == 1)
            {
                noalias(p) = r;
                alpha = 2.0/d;
            }
            else
            {
                beta = pow(c*alpha/2, 2);
                alpha = 1.0/(d - beta);
                TSparseSpaceType::ScaleAndAdd(1.0, r, beta, p); //p=r+beta*p
            }
/*            KRATOS_WATCH(alpha)*/
/*            KRATOS_WATCH(beta)*/

            TSparseSpaceType::UnaliasedAdd(rX, alpha, p); // x=x+alpha*p
            TSparseSpaceType::Mult(rA, rX, r); //r=A*x
            TSparseSpaceType::ScaleAndAdd(1.0, rB, -1.0, r); //r=b-A*x
        }

        return true;
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
        //GetTimeTable()->Start(Info());

        // BaseType::GetPreconditioner()->Initialize(rA, rX, rB);

        bool is_solved = true;
        VectorType x(TDenseSpaceType::Size1(rX));
        VectorType b(TDenseSpaceType::Size1(rB));
        for(unsigned int i = 0 ; i < TDenseSpaceType::Size2(rX) ; ++i)
        {
            TDenseSpaceType::GetColumn(i, rX, x);
            TDenseSpaceType::GetColumn(i, rB, b);

            Solve(rA, x, b);

            noalias(column(rX, i)) = x;
            noalias(column(rB, i)) = b;
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
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Chebyshev iterative solver with " << BaseType::GetPreconditioner()->Info();
        return buffer.str();
    }

    /// Print information about this object.
    void  PrintInfo(std::ostream& OStream) const
    {
        OStream << Info();
    }

    /// Print object's data.
    void  PrintData(std::ostream& OStream) const
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

    double mLambdaMin;
    double mLambdaMax;
    std::string mMinEigenvalueMode;
    std::string mMaxEigenvalueMode;

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

}; // Class ChebyshevSolver

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType,
         class TPreconditionerType,
         class TReordererType>
inline std::istream& operator >> (std::istream& IStream,
                                  ChebyshevSolver<TSparseSpaceType, TDenseSpaceType,
                                  TPreconditionerType, TReordererType>& rThis)
{
    return IStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType,
         class TPreconditionerType,
         class TReordererType>
inline std::ostream& operator << (std::ostream& OStream,
                                  const ChebyshevSolver<TSparseSpaceType, TDenseSpaceType,
                                  TPreconditionerType, TReordererType>& rThis)
{
    rThis.PrintInfo(OStream);
    OStream << std::endl;
    rThis.PrintData(OStream);

    return OStream;
}
///@}


}  // namespace Kratos.

#endif //  KRATOS_MULTITHREADED_SOLVERS_APPLICATION_CHEBYSHEV_SOLVER_H_INCLUDED  defined
