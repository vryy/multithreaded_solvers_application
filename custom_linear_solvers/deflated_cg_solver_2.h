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
//   Date:                $Date: 1 July 2015 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(MULTITHREADED_SOLVERS_APP_DEFLATED_CG_SOLVER_2_H_INCLUDED )
#define  MULTITHREADED_SOLVERS_APP_DEFLATED_CG_SOLVER_2_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <vector>
#include <set>

// External includes


// Project includes
#include "includes/define.h"
#include "linear_solvers/iterative_solver.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"
#include "utilities/deflation_utils.h"

#define PRINT_W_DISTRIBUTION

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
class DeflatedCGSolver2 : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DeflatedCGSolver2
    KRATOS_CLASS_POINTER_DEFINITION(DeflatedCGSolver2);

    typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType;

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> LinearSolverType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DeflatedCGSolver2(double NewMaxTolerance, unsigned int NewMaxIterationsNumber, typename TPreconditionerType::Pointer pNewPreconditioner)
    : BaseType(NewMaxTolerance, NewMaxIterationsNumber, pNewPreconditioner)
    {
    }

    /// Copy constructor.

    DeflatedCGSolver2(const DeflatedCGSolver2& Other) : BaseType(Other)
    {
    }


    /// Destructor.
    virtual ~DeflatedCGSolver2()
    {
    }


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.

    DeflatedCGSolver2 & operator=(const DeflatedCGSolver2& Other)
    {
        BaseType::operator=(Other);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    virtual bool AdditionalPhysicalDataIsNeeded()
    {
        return this->GetPreconditioner()->AdditionalPhysicalDataIsNeeded();
    }


    virtual void ProvideAdditionalData(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB,
        typename ModelPart::DofsArrayType& rdof_set,
        ModelPart& r_model_part
    )
    {
        this->GetPreconditioner()->ProvideAdditionalData(rA, rX, rB, rdof_set, r_model_part);
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
        if (this->IsNotConsistent(rA, rX, rB))
            return false;

//        GetTimeTable()->Start(Info());

        BaseType::GetPreconditioner()->Initialize(rA,rX,rB);
        BaseType::GetPreconditioner()->ApplyInverseRight(rX);
        BaseType::GetPreconditioner()->ApplyLeft(rB);

        bool is_solved = IterativeSolve(rA, rX, rB);
        
        if(!is_solved)
            std::cout << "Warning: the iterative solver is not successful" << std::endl;

        BaseType::GetPreconditioner()->Finalize(rX);

//        GetTimeTable()->Stop(Info());

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
        std::cout << "************ DeflatedCGSolver2::Solve(SparseMatrixType&, DenseMatrixType&, DenseMatrixType&) not defined! ************" << std::endl;
        return false;
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

    /// Turn back information as a string.

    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Deflated Conjugate gradient linear solver (v2) with " << BaseType::GetPreconditioner()->Info();
        return buffer.str();
    }

    /// Print information about this object.

    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.

    virtual void PrintData(std::ostream& rOStream) const
    {
        BaseType::PrintData(rOStream);
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

    virtual void ConstructW(SparseMatrixType& rA, std::vector<int>& w, SparseMatrixType& Adeflated) const
    {
        int max_reduced_size = 100;
        std::cout << "Construct the deflation space using standard aggregation assuming 1 dof per node, max_reduced_size = " << max_reduced_size << std::endl;
        DeflationUtils::ConstructW(max_reduced_size, rA, w, Adeflated);
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


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /// Perform the multiplication rY = W (W^T K W)^(-1) W^T * rX
    /// Deflation matrix W is represented by indices vector w
    void ApplyE(std::size_t full_size,
                std::size_t reduced_size,
                std::vector<int>& w,
                VectorType& rX,
                VectorType& rY,
                LUSkylineFactorization<TSparseSpaceType, TDenseSpaceType>& rFact)
    {
        VectorType t(reduced_size), d(reduced_size);
        DeflationUtils::ApplyWtranspose(w, rX, t);
        rFact.backForwardSolve(reduced_size, t, d);
        DeflationUtils::ApplyW(w, d, rY);
    }

    /// Perform the multiplication rY = P * rX, where P is the projection matrix P = I - KW(W^T K W)^(-1)W^T
    /// Deflation matrix W is represented by indices vector w
    void ApplyP(std::size_t full_size,
                std::size_t reduced_size,
                std::vector<int>& w,
                SparseMatrixType& rA,
                VectorType& rX,
                VectorType& rY,
                LUSkylineFactorization<TSparseSpaceType, TDenseSpaceType>& rFact)
    {
        VectorType v(full_size);
        ApplyE(full_size, reduced_size, w, rX, v, rFact);
        TSparseSpaceType::Mult(rA, v, rY);
        TSparseSpaceType::ScaleAndAdd(1.00, rX, -1.00, rY);
    }

    /// Perform the multiplication rY = P^T * rX, where P is the projection matrix P = I - KW(W^T K W)^(-1)W^T
    /// Deflation matrix W is represented by indices vector w
    void ApplyPtranspose(std::size_t full_size,
                         std::size_t reduced_size,
                         std::vector<int>& w,
                         SparseMatrixType& rA,
                         VectorType& rX,
                         VectorType& rY,
                         LUSkylineFactorization<TSparseSpaceType, TDenseSpaceType>& rFact)
    {
        VectorType v(full_size);
        TSparseSpaceType::TransposeMult(rA, rX, v);
        ApplyE(full_size, reduced_size, w, v, rY, rFact);
        TSparseSpaceType::ScaleAndAdd(1.00, rX, -1.00, rY);
    }

    /// REFERENCE: Jönsthövel, Comparison of the deflated preconditioned conjugate gradient method and algebraic multigrid for composite materials
    bool IterativeSolve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
        const int full_size = TSparseSpaceType::Size(rX);

        std::vector<int> w; // use to build the deflation subspace, where w[i] = 1 means that the row i is accounted in deflation space and vice versa
        SparseMatrixType Adeflated;

        // construct the reduction structure
        this->ConstructW(rA, w, Adeflated);

        //fill reduced matrix Adeflated
        DeflationUtils::FillDeflatedMatrix(rA, w, Adeflated);

        std::size_t reduced_size = Adeflated.size1();
        KRATOS_WATCH(reduced_size)

        #ifdef PRINT_W_DISTRIBUTION
//        std::cout << "w:";
//        for(unsigned int i = 0; i < w.size(); ++i)
//            std::cout << " " << w[i];
//        std::cout << std::endl;
        std::vector<unsigned int> w_dist(reduced_size);
        for(unsigned int i = 0; i < w.size(); ++i)
            if(w[i] >= reduced_size || w[i] < 0)
                KRATOS_THROW_ERROR(std::logic_error, "w contains strange entry at location", i)
            else
                ++w_dist[w[i]];
        std::cout << "w distribution:" << std::endl;
        for(unsigned int i = 0; i < reduced_size; ++i)
            std::cout << "\tw[i] = " << i << ": " << w_dist[i] << " elements" << std::endl;
        std::cout << std::endl;
        #endif

        // To save some time, we do the factorization once, and do the solve several times.
        // When this is available through the LinearSolver interface, replace this.
        LUSkylineFactorization<TSparseSpaceType, TDenseSpaceType> Factorization;
        Factorization.copyFromCSRMatrix(Adeflated);
        Factorization.factorize();
        std::cout << "********** Factorization of deflated matrix done!" << std::endl;

        VectorType r(full_size), rh(full_size), y(full_size), p(full_size), q(full_size), wh(full_size), uh(full_size);

        // r = b - A x
        TSparseSpaceType::Mult(rA, rX, r);
        TSparseSpaceType::ScaleAndAdd(1.00, rB, -1.00, r);

        // rh = P * r
        ApplyP(full_size, reduced_size, w, rA, r, rh, Factorization);
        double normrh0 = TSparseSpaceType::TwoNorm(r);

        // solve A y = rh
        this->PreconditionedMult(rA, rh, y);

        // assign p = y;
        TSparseSpaceType::Assign(p, 1.0, y);

        // temporarily set uh to zero. I don't know how uh_0 is calculated
        TSparseSpaceType::SetToZero(uh);

        boost::progress_display show_progress( this->GetMaxIterationsNumber() );
        BaseType::mIterationsNumber = 0;
        double eps;
        do
        {
            TSparseSpaceType::Mult(rA, p, q);

            ApplyP(full_size, reduced_size, w, rA, q, wh, Factorization);

            double alpha = TSparseSpaceType::Dot(rh, y) / TSparseSpaceType::Dot(wh, p);

            TSparseSpaceType::UnaliasedAdd(uh, alpha, p);

            double temp = TSparseSpaceType::Dot(rh, y);

            TSparseSpaceType::UnaliasedAdd(rh, -alpha, wh);

            this->PreconditionedMult(rA, rh, y);

            double beta = TSparseSpaceType::Dot(rh, y) / temp;

            TSparseSpaceType::ScaleAndAdd(1.0, y, beta, p);

            eps = TSparseSpaceType::TwoNorm(rh) / normrh0;
            
            std::cout << "iter " << BaseType::mIterationsNumber << ": eps = " << eps << std::endl;
            
            ++show_progress;
        }
        while(++BaseType::mIterationsNumber < this->GetMaxIterationsNumber() && eps > this->GetTolerance());

        ApplyE(full_size, reduced_size, w, rB, y, Factorization);
        ApplyPtranspose(full_size, reduced_size, w, rA, uh, rX, Factorization);
        TSparseSpaceType::UnaliasedAdd(rX, 1.0, y);

        return true;
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

}; // Class DeflatedCGSolver2

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
inline std::istream & operator >>(std::istream& IStream,
                                  DeflatedCGSolver2<TSparseSpaceType, TDenseSpaceType,
                                  TPreconditionerType, TReordererType>& rThis)
{
    return IStream;
}

/// output stream function

template<class TSparseSpaceType, class TDenseSpaceType,
         class TPreconditionerType,
         class TReordererType>
inline std::ostream & operator <<(std::ostream& OStream,
                                  const DeflatedCGSolver2<TSparseSpaceType, TDenseSpaceType,
                                  TPreconditionerType, TReordererType>& rThis)
{
    rThis.PrintInfo(OStream);
    OStream << std::endl;
    rThis.PrintData(OStream);

    return OStream;
}
///@}


} // namespace Kratos.

#undef PRINT_W_DISTRIBUTION

#endif // MULTITHREADED_SOLVERS_APP_DEFLATED_SUBDOMAIN_CG_SOLVER_H_INCLUDED  defined 


