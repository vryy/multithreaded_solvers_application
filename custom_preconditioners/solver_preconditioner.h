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
//   Date:                $Date: 2013-14-1 9:34:00 $
//   Revision:            $Revision: 1.1 $
//
//

#if !defined(KRATOS_SOLVER_PRECONDITIONER_H_INCLUDED)
#define KRATOS_SOLVER_PRECONDITIONER_H_INCLUDED




// System includes



// External includes
#include <boost/smart_ptr.hpp>
#include <boost/numeric/ublas/vector.hpp>


// Project includes
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "linear_solvers/preconditioner.h"
#include "linear_solvers/linear_solver.h"



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

/// SolverPreconditioner class.
/**   */
template<class TSparseSpaceType, class TDenseSpaceType>
class SolverPreconditioner : public Preconditioner<TSparseSpaceType, TDenseSpaceType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SolverPreconditioner
    KRATOS_CLASS_POINTER_DEFINITION (SolverPreconditioner);

    typedef Preconditioner<TSparseSpaceType, TDenseSpaceType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
    
    typedef typename TSparseSpaceType::MatrixPointerType SparseMatrixPointerType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType> LinearSolverType;
    
    typedef typename LinearSolverType::Pointer LinearSolverPointerType;

    typedef std::size_t  SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SolverPreconditioner(LinearSolverPointerType pSolver)
    {
        mpSolver = pSolver;
    }


    /// Copy constructor.
    SolverPreconditioner(const SolverPreconditioner& Other)
    {
        mpSolver = Other.mpSolver;
    }


    /// Destructor.
    virtual ~SolverPreconditioner()
    {
    }



    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    SolverPreconditioner& operator=(const SolverPreconditioner& Other)
    {
        mpSolver = Other.mpSolver;
        return *this;
    }




    ///@}
    ///@name Operations
    ///@{

    virtual void Initialize(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
        mpSolver->Initialize(rA, rX, rB);
        mA = rA;
    }
    
    
    virtual bool AdditionalPhysicalDataIsNeeded()
    {
        return false;
    }


    virtual void Mult(SparseMatrixType& rA, VectorType& rX, VectorType& rY)
    {
        TSparseSpaceType::Mult(rA, rX, rY);
        ApplyLeft(rY);
    }


    /** calculate preconditioned_X = A^{-1} * X;
        @param rX  Unknows of preconditioner suystem
    */
    virtual VectorType& ApplyLeft(VectorType& rX)
    {
        SizeType size = TSparseSpaceType::Size(rX);
        
        VectorType pX(size, 0.00);
        
        SparseMatrixType tmpA = mA; // to avoid mA is modified
        
        mpSolver->Solve(tmpA, pX, rX);
        
//        KRATOS_WATCH(pX);
        
        TSparseSpaceType::Copy(pX, rX);
        
        return rX;
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
        buffer << "SolverPreconditioner using ";
        buffer << mpSolver->Info();
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

    LinearSolverPointerType mpSolver;

    SparseMatrixType mA;

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
    ///@name Unaccessible methods
    ///@{


    ///@}

}; // Class SolverPreconditioner

///@}


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType>
inline std::istream& operator >> (std::istream& IStream, SolverPreconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
{
    return IStream;
}


/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType>
inline std::ostream& operator << (std::ostream& OStream, const SolverPreconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
{
    rThis.PrintInfo(OStream);
    OStream << std::endl;
    rThis.PrintData(OStream);
    return OStream;
}
///@}


}  // namespace Kratos.


#endif // KRATOS_SOLVER_PRECONDITIONER_H_INCLUDED  defined 

