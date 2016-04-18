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


#if !defined(KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BICGSTAB_SOLVER_H_INCLUDED )
#define KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BICGSTAB_SOLVER_H_INCLUDED



// System includes
#include <math.h>


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "linear_solvers/iterative_solver.h"
#include "custom_utilities/feast_solver.h"
#include "custom_utilities/arpack_solver.h"

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
class BicgstabSolver : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of BicgstabSolver
    KRATOS_CLASS_POINTER_DEFINITION(BicgstabSolver);

    typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BicgstabSolver() {}

    BicgstabSolver(double NewTolerance, unsigned int NewMaxIterationsNumber, typename TPreconditionerType::Pointer pNewPreconditioner) : BaseType(NewTolerance, NewMaxIterationsNumber, pNewPreconditioner) , mPrecondType("Left")
    {}

    BicgstabSolver(double NewTolerance, unsigned int NewMaxIterationsNumber, std::string PrecondType, typename TPreconditionerType::Pointer pNewPreconditioner) : BaseType(NewTolerance, NewMaxIterationsNumber, pNewPreconditioner) , mPrecondType(PrecondType)
    {}

    /// Copy constructor.
    BicgstabSolver(const BicgstabSolver& Other) : BaseType(Other) {}

    /// Destructor.
    virtual ~BicgstabSolver() {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    BicgstabSolver& operator=(const BicgstabSolver& Other)
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

        int ierr;
        if(mPrecondType == std::string("Left"))
        {
            BaseType::GetPreconditioner()->Initialize(rA, rX, rB);
            ierr = IterativeSolveLeft(rA, rX, rB);
            if(ierr != 0)
                std::cout << "Warning: the iterative solver encountered some problem, error code = " << ierr << std::endl;
        }
        else if(mPrecondType == std::string("Right"))
        {
            KRATOS_THROW_ERROR(std::logic_error, Info(), " does not support right preconditioning yet")
        }
        else if(mPrecondType == std::string("Left-Right"))
        {
            BaseType::GetPreconditioner()->Initialize(rA, rX, rB);
            BaseType::GetPreconditioner()->ApplyInverseRight(rX); // x0_tilde = K2 * x0
            BaseType::GetPreconditioner()->ApplyLeft(rB); // b_tilde = K1^(-1) * b;

            ierr = IterativeSolveLeftRight(rA, rX, rB);
            
            if(ierr != 0)
                std::cout << "Warning: the iterative solver is not successful, error code = " << ierr << std::endl;

            BaseType::GetPreconditioner()->ApplyRight(rX); // x = K2^(-1) * x_tilde;
        }

//        GetTimeTable()->Stop(Info());

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
    bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
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
        buffer << "Biconjugate-gradient stabilized iterative solver with " << BaseType::GetPreconditioner()->Info()
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

    std::string mPrecondType;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /// Iterative solve with left preconditioning
    /// REF: https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method
    int IterativeSolveLeft(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
        double resid, tol = this->GetTolerance();
        unsigned int max_iter = this->GetMaxIterationsNumber();
	    int i, j = 1, k, size = rX.size();
	    double rho_1, rho_2, alpha, beta, omega, norms, normr, normb;
	    VectorType p(size), phat(size), s(size), shat(size), t(size), v(size), r(size), rtilde(size);

	    normb = TSparseSpaceType::TwoNorm(rB);
	    TSparseSpaceType::Mult(rA, rX, r); //r=A*x
	    TSparseSpaceType::ScaleAndAdd(1.0, rB, -1.0, r); //r=b-A*x
	    TSparseSpaceType::Assign(rtilde, 1.0, r);

	    if (normb == 0.0)
		    normb = 1.0;

        normr = TSparseSpaceType::TwoNorm(r);
	    if ((resid = normr / normb) <= tol) {
		    tol = resid;
		    max_iter = 0;
		    this->SetIterationsNumber(max_iter);
		    BaseType::mBNorm = normb;
		    this->SetResidualNorm(resid);
		    return 0;
	    }

        boost::progress_display show_progress( max_iter );
	    for (int i = 1; i <= max_iter; ++i) {
	        rho_1 = TSparseSpaceType::Dot(rtilde, r);
//	        KRATOS_WATCH(normr)
	        if(rho_1 == 0.0) {
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
//	        KRATOS_WATCH(norms)
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
	        this->GetPreconditioner()->ApplyLeft(shat);
	        TSparseSpaceType::Mult(rA, shat, t);
	        omega = TSparseSpaceType::Dot(t, s) / TSparseSpaceType::Dot(t, t);
	        TSparseSpaceType::UnaliasedAdd(rX, alpha, phat);
	        TSparseSpaceType::UnaliasedAdd(rX, omega, shat);
	        TSparseSpaceType::ScaleAndAdd(1.0, s, -omega, t, r);
	        
	        rho_2 = rho_1;
	        normr = TSparseSpaceType::TwoNorm(r);
	        if((resid = normr / normb) < tol) {
	            tol = resid;
	            max_iter = i;
	            this->SetIterationsNumber(max_iter);
		        BaseType::mBNorm = normb;
		        this->SetResidualNorm(resid);
	            return 0;
	        }
	        if(omega == 0.0) {
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

    /// Iterative solve with left-right preconditioning
    /// REF: http://www.netlib.org/linalg/html_templates/node54.html
    ///      https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method
    int IterativeSolveLeftRight(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
        double resid, tol = this->GetTolerance();
        unsigned int max_iter = this->GetMaxIterationsNumber();
	    int i, j = 1, k, size = rX.size();
	    double rho_1, rho_2, alpha, beta, omega, norms, normr, normb;
	    VectorType p(size), s(size), t(size), v(size), r(size), rtilde(size);

	    normb = TSparseSpaceType::TwoNorm(rB);

	    TSparseSpaceType::Mult(rA, rX, r); //r=A*x_tilde
	    this->GetPreconditioner()->ApplyLeft(r); //r=K1^(-1)*A*x_tilde
	    TSparseSpaceType::ScaleAndAdd(1.0, rB, -1.0, r); //r=b_tilde-K1^(-1)*A*x_tilde
	    TSparseSpaceType::Assign(rtilde, 1.0, r);

	    if (normb == 0.0)
		    normb = 1.0;

        normr = TSparseSpaceType::TwoNorm(r);
	    if ((resid = normr / normb) <= tol) {
		    tol = resid;
		    max_iter = 0;
		    this->SetIterationsNumber(max_iter);
		    BaseType::mBNorm = normb;
		    this->SetResidualNorm(resid);
		    return 0;
	    }

        boost::progress_display show_progress( max_iter );
	    for (int i = 1; i <= max_iter; ++i) {
	        rho_1 = TSparseSpaceType::Dot(rtilde, r);
//	        KRATOS_WATCH(normr)
	        if(rho_1 == 0.0) {
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

            this->PreconditionedMult(rA, p, v);
	        alpha = rho_1 / TSparseSpaceType::Dot(rtilde, v);
	        TSparseSpaceType::ScaleAndAdd(1.0, r, -alpha, v, s);
	        norms = TSparseSpaceType::TwoNorm(s);
//	        KRATOS_WATCH(norms)
	        if((resid = norms / normb) < tol) {
	            TSparseSpaceType::UnaliasedAdd(rX, alpha, p);
	            tol = resid;
	            max_iter = i;
	            this->SetIterationsNumber(max_iter);
		        BaseType::mBNorm = normb;
		        this->SetResidualNorm(resid);
	            return 0;
	        }

            this->PreconditionedMult(rA, s, t);
	        omega = TSparseSpaceType::Dot(t, s) / TSparseSpaceType::Dot(t, t);
	        TSparseSpaceType::UnaliasedAdd(rX, alpha, p);
	        TSparseSpaceType::UnaliasedAdd(rX, omega, s);
	        TSparseSpaceType::ScaleAndAdd(1.0, s, -omega, t, r);

	        rho_2 = rho_1;
	        normr = TSparseSpaceType::TwoNorm(r);
	        if((resid = normr / normb) < tol) {
	            tol = resid;
	            max_iter = i;
	            this->SetIterationsNumber(max_iter);
		        BaseType::mBNorm = normb;
		        this->SetResidualNorm(resid);
	            return 0;
	        }
	        if(omega == 0.0) {
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

}; // Class BicgstabSolver

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
                                  BicgstabSolver<TSparseSpaceType, TDenseSpaceType,
                                  TPreconditionerType, TReordererType>& rThis)
{
    return IStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType,
         class TPreconditionerType,
         class TReordererType>
inline std::ostream& operator << (std::ostream& OStream,
                                  const BicgstabSolver<TSparseSpaceType, TDenseSpaceType,
                                  TPreconditionerType, TReordererType>& rThis)
{
    rThis.PrintInfo(OStream);
    OStream << std::endl;
    rThis.PrintData(OStream);

    return OStream;
}
///@}


}  // namespace Kratos.

#endif //  KRATOS_MULTITHREADED_SOLVERS_APPLICATION_BICGSTAB_SOLVER_H_INCLUDED  defined 


