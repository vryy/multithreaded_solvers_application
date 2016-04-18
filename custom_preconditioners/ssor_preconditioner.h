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
//   Date:                $Date: 2 July 2015 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_SSOR_PRECONDITIONER_H_INCLUDED )
#define  KRATOS_SSOR_PRECONDITIONER_H_INCLUDED



// System includes


// External includes


// Project includes
#include "includes/define.h"


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

/// SSORPreconditioner class.
/** SSORPreconditioner for linesr system solvers.
REFERENCE: http://www.netlib.org/linalg/html_templates/node58.html
 */
template<class TSparseSpaceType, class TDenseSpaceType>
class SSORPreconditioner : public Preconditioner<TSparseSpaceType, TDenseSpaceType>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(SSORPreconditioner);

    typedef  Preconditioner<TSparseSpaceType, TDenseSpaceType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SSORPreconditioner(double Omega) : mOmega(Omega)
    {
        if(mOmega >= 2.0 || mOmega <= 0.0)
            KRATOS_THROW_ERROR(std::logic_error, "The Omega for SSOR must be in [0.0 2.0]", "")
    }

    /// Copy constructor.
    SSORPreconditioner(const SSORPreconditioner& Other)
    : BaseType(Other), mD(Other.mD) {}

    /// Destructor.
    virtual ~SSORPreconditioner() {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    SSORPreconditioner& operator=(const SSORPreconditioner& Other)
    {
        BaseType::operator=(Other);
        mD = Other.mD;
        return *this;
    }


    ///@}
    ///@name Operations
    ///@{


    /** Initialize
    Initialize preconditioner for linear system rA*rX=rB
    @param rA  system matrix
    @param rX Unknows vector
    @param rB Right side linear system of equations
    */
    void Initialize(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
        // clear data
        TSparseSpaceType::ClearData(mD);
        TSparseSpaceType::ClearData(mL);
        TSparseSpaceType::ClearData(mU);
    
        // fill the diagonal
        unsigned int n = TSparseSpaceType::Size1(rA);
        mD.resize(n);
        #pragma omp parallel for
        for(unsigned int i = 0 ; i < n; ++i)
        {
            double diag_Aii = rA(i, i);
            if(diag_Aii != 0.0)
                mD[i] = diag_Aii * (2.0 - mOmega) / mOmega;
            else
                KRATOS_THROW_ERROR(std::logic_error,"zero found in the diagonal. Diagonal preconditioner can not be used","");
        }

        // pointer of the matrix data
        const std::size_t* index1 = rA.index1_data().begin();
        const std::size_t* index2 = rA.index2_data().begin();
        const double*      values = rA.value_data().begin();

        // fill the upper and lower matrix
        TSparseSpaceType::Resize(mU, n, n);
        TSparseSpaceType::Resize(mL, n, n);
        TSparseSpaceType::SetToZero(mU);
        TSparseSpaceType::SetToZero(mL);
        for(unsigned int i = 0; i < n; ++i)
            for(unsigned int j = index1[i]; j < index1[i + 1]; ++j)
            {
                if(index2[j] == i)
                {
                    mL.push_back(i, i, values[j] / mOmega);
                    mU.push_back(i, i, values[j] / mOmega);
                }
                else if(index2[j] > i)
                    mU.push_back(i, index2[j], values[j]);
                else if(index2[j] < i)
                    mL.push_back(i, index2[j], values[j]);
            }
        mU.complete_index1_data();
        mL.complete_index1_data();
    }

    void Initialize(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
    {
        BaseType::Initialize(rA, rX, rB);
    }

    VectorType& ApplyLeft(VectorType& rX)
    {
        std::size_t j, n = TSparseSpaceType::Size(rX);

        // apply inverse L by back solve
        const std::size_t* index1 = mL.index1_data().begin();
        const std::size_t* index2 = mL.index2_data().begin();
        const double*      values = mL.value_data().begin();
        for(std::size_t i = 0; i < n; ++i)
        {
            double sum = rX(i);
            for(j = index1[i] + 1; j < index1[i + 1]; ++j)
                sum -= values[j] * rX(index2[j]);
            rX(i) = sum / mL(i, i);
        }

        // apply the diagonal
        #pragma omp parallel for
        for(std::size_t i = 0; i < n; ++i)
            rX(i) *= mD(i);

        // apply inverse U by back solve
        index1 = mU.index1_data().begin();
        index2 = mU.index2_data().begin();
        values = mU.value_data().begin();
        for(int i = n - 1; i >= 0; --i)
        {
            double sum = rX(i);
            for(j = index1[i] + 1; j < index1[i + 1]; ++j)
                sum -= values[j] * rX(index2[j]);
            rX(i) = sum / mU(i, i);
        }

        return rX;
    }

    VectorType& ApplyRight(VectorType& rX)
    {
        return rX;
    }

    VectorType& ApplyTransposeLeft(VectorType& rX)
    {
        return rX;
    }

    VectorType& ApplyTransposeRight(VectorType& rX)
    {
        return rX;
    }

    VectorType& ApplyInverseRight(VectorType& rX)
    {
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
        std::stringstream ss;
        ss << "SSOR preconditioner with omega = " << mOmega;
        return ss.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& OStream) const
    {
        OStream << Info();
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

    VectorType mD;

    SparseMatrixType mU;

    SparseMatrixType mL;

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

}; // Class SSORPreconditioner

///@}

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType>
inline std::istream& operator >> (std::istream& IStream,
                                  SSORPreconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
{
    return IStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType>
inline std::ostream& operator << (std::ostream& OStream,
                                  const SSORPreconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
{
    rThis.PrintInfo(OStream);
    OStream << std::endl;
    rThis.PrintData(OStream);

    return OStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_SSOR_PRECONDITIONER_H_INCLUDED  defined 

