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
//   Date:                $Date: 28 Aug 2018 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_MULTITHREADED_SOLVERS_APPLICATION_RICHARDSON_SOLVER_H_INCLUDED )
#define KRATOS_MULTITHREADED_SOLVERS_APPLICATION_RICHARDSON_SOLVER_H_INCLUDED



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
class RichardsonSolver : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of RichardsonSolver
    KRATOS_CLASS_POINTER_DEFINITION(RichardsonSolver);

    typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RichardsonSolver(double Omega) : mOmega(Omega), mPrecondType("Left"), BaseType() {}

    RichardsonSolver(double Omega, double NewTolerance, unsigned int NewMaxIterationsNumber)
    : mOmega(Omega), BaseType(NewTolerance, NewMaxIterationsNumber) , mPrecondType("Left")
    {}

    RichardsonSolver(double Omega, double NewTolerance, unsigned int NewMaxIterationsNumber, typename TPreconditionerType::Pointer pNewPreconditioner)
    : mOmega(Omega), BaseType(NewTolerance, NewMaxIterationsNumber, pNewPreconditioner) , mPrecondType("Left")
    {}

    RichardsonSolver(double Omega, double NewTolerance, unsigned int NewMaxIterationsNumber, std::string PrecondType, typename TPreconditionerType::Pointer pNewPreconditioner)
    : mOmega(Omega), BaseType(NewTolerance, NewMaxIterationsNumber, pNewPreconditioner) , mPrecondType(PrecondType)
    {}

    /// Copy constructor.
    RichardsonSolver(const RichardsonSolver& Other) : mOmega(Other.mOmega), BaseType(Other) {}

    /// Destructor.
    virtual ~RichardsonSolver() {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    RichardsonSolver& operator=(const RichardsonSolver& Other)
    {
        BaseType::operator=(Other);
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
    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
        if(this->IsNotConsistent(rA, rX, rB))
            return false;

//        GetTimeTable()->Start(Info());

        if(mPrecondType == std::string("Left"))
        {
            if (this->GetPreconditioner() != NULL)
                BaseType::GetPreconditioner()->Initialize(rA, rX, rB);
        }
        else
        {
            std::stringstream ss;
            ss << Info() << " does not yet support " << mPrecondType << " preconditioning";
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
        }

        double tol = this->GetTolerance(), resid, normb, normr;
        unsigned int max_iter = this->GetMaxIterationsNumber();
        std::size_t size = rX.size();
        VectorType r(size);

        normb = TSparseSpaceType::TwoNorm(rB);
        if (normb == 0.0)
            normb = 1.0;

        for (unsigned int i = 1; i <= max_iter; ++i)
        {
            TSparseSpaceType::Mult(rA, rX, r); //r=A*x
            TSparseSpaceType::ScaleAndAdd(1.0, rB, -1.0, r); //r=b-A*x

            normr = TSparseSpaceType::TwoNorm(r);

            if (this->GetEchoLevel() > 0)
                std::cout << "RichardsonSolver iteration #" << i
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
                this->GetPreconditioner()->ApplyLeft(r); //r=P^(-1)*r
            TSparseSpaceType::ScaleAndAdd(mOmega, r, 1.0, rX);
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
        buffer << "Richardson iterative solver with " << BaseType::GetPreconditioner()->Info()
               << ", Preconditioner type = " << mPrecondType;
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

    double mOmega;
    std::string mPrecondType;

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

}; // Class RichardsonSolver

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
                                  RichardsonSolver<TSparseSpaceType, TDenseSpaceType,
                                  TPreconditionerType, TReordererType>& rThis)
{
    return IStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType,
         class TPreconditionerType,
         class TReordererType>
inline std::ostream& operator << (std::ostream& OStream,
                                  const RichardsonSolver<TSparseSpaceType, TDenseSpaceType,
                                  TPreconditionerType, TReordererType>& rThis)
{
    rThis.PrintInfo(OStream);
    OStream << std::endl;
    rThis.PrintData(OStream);

    return OStream;
}
///@}


}  // namespace Kratos.

#endif //  KRATOS_MULTITHREADED_SOLVERS_APPLICATION_RICHARDSON_SOLVER_H_INCLUDED  defined
