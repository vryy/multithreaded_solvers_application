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


#if !defined(KRATOS_ILU_LR_PRECONDITIONER_H_INCLUDED )
#define  KRATOS_ILU_LR_PRECONDITIONER_H_INCLUDED




// System includes



// External includes
#include "boost/smart_ptr.hpp"



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

/// ILU_LR_Preconditioner class.
/**   */
template<class TSparseSpaceType, class TDenseSpaceType, class TModelPartType>
class ILU_LR_Preconditioner : public Preconditioner<TSparseSpaceType, TDenseSpaceType, TModelPartType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of ILU_LR_Preconditioner
    typedef boost::shared_ptr<ILU_LR_Preconditioner> Pointer;

    typedef Preconditioner<TSparseSpaceType, TDenseSpaceType, TModelPartType> BaseType;

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
    ILU_LR_Preconditioner()
    {
        L = NULL;
        iL = NULL;
        jL = NULL;
        U = NULL;
        iU = NULL;
        jU = NULL;
    }

    /// Copy constructor.
    ILU_LR_Preconditioner(const ILU_LR_Preconditioner& Other) {}

    /// Destructor.
    ~ILU_LR_Preconditioner() override
    {
        if ( L!=NULL) delete[]  L;
        if (iL!=NULL) delete[] iL;
        if (jL!=NULL) delete[] jL;
        if ( U!=NULL) delete[]  U;
        if (iU!=NULL) delete[] iU;
        if (jU!=NULL) delete[] jU;

        L = NULL;
        iL = NULL;
        jL = NULL;
        U = NULL;
        iU = NULL;
        jU = NULL;
    }

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    ILU_LR_Preconditioner& operator=(const ILU_LR_Preconditioner& Other)
    {
        mILUSize = Other.mILUSize;
        unsigned int size = Other.iL[mILUSize];
        L = new double[size];
        U = new double[size];
        iL = new int[mILUSize+1];
        jL = new int[size];
        iU = new int[mILUSize+1];
        jU = new int[size];


        std::copy(Other.L, Other.L+size, L);
        std::copy(Other.U, Other.U+size, U);
        std::copy(Other.iL, Other.iL+mILUSize+1, iL);
        std::copy(Other.jL, Other.jL+size, jL);
        std::copy(Other.iU, Other.iU+mILUSize+1, iU);
        std::copy(Other.jU, Other.jU+size, jU);

        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /** multiply rX by L^-1
        @param rX  Unknows of preconditioner suystem
    */
    VectorType& ApplyLeft(VectorType& rX) override
    {
        const int size = TSparseSpaceType::Size(rX);
        double sum;
        int i, indexj;
        for (i=0; i<size; i++)
        {
            sum=rX[i];
            for (indexj=iL[i]; indexj<iL[i+1]; indexj++)
            {
                sum=sum-L[indexj]*rX[jL[indexj]];
            }
            rX[i]=sum;
        }
        return rX;
    }

    /** multiply rX by U^-1
        @param rX  Unknows of preconditioner suystem
    */
    VectorType& ApplyRight(VectorType& rX) override
    {
        const int size = TSparseSpaceType::Size(rX);
        double sum;
        int i, indexj;
        for (i=size-1; i>=0; i--)
        {
            sum=rX[i];
            for (indexj=iU[i]+1; indexj<iU[i+1]; indexj++)
            {
                sum=sum-U[indexj]*rX[jU[indexj]];
            }
            rX[i]=sum/U[iU[i]];
        }
        return rX;
    }

    /** multiply rX by U
        @param rX  Unknows of preconditioner suystem
        REMARKS: To be debugged
    */
    VectorType& ApplyInverseRight(VectorType& rX) override
    {
        const int size = TSparseSpaceType::Size(rX);
        VectorType temp(size);
        double sum;
        int i, indexj;
        for (i=0; i<size; i++)
        {
            temp[i] = 0.0;
            for (indexj=iU[i]+1; indexj<iU[i+1]; indexj++)
            {
                temp[i]+=U[indexj]*rX[jU[indexj]];
            }
        }
        for (i=0; i<size; i++) rX[i] = temp[i];
        return rX;
    }

    /** Multiply rX by U^-T
        @param rX  Unknows of preconditioner suystem
    */
    VectorType& ApplyTransposeLeft(VectorType& rX) override
    {
        const int size = TSparseSpaceType::Size(rX);
        int i, indexj;
        double tempi;
        for (i=0; i<size; i++)
        {
            rX[i]=rX[i]/U[iU[i]];
            tempi=rX[i];
            for (indexj=iU[i]+1; indexj<iU[i+1]; indexj++)
            {
                rX[jU[indexj]]=rX[jU[indexj]]-tempi*U[indexj];
            }
        }
        return rX;
    }

    /** Multiply rX by L^-T
        @param rX  Unknows of preconditioner suystem
    */
    VectorType& ApplyTransposeRight(VectorType& rX) override
    {
        const int size = TSparseSpaceType::Size(rX);
        int i, indexj;
        double rxi;
        for (i=size-1; i>=0; i--)
        {
            rxi=rX[i];
            for (indexj=iL[i]; indexj<iL[i+1]; indexj++)
            {
                rX[jL[indexj]]=rX[jL[indexj]]-rxi*L[indexj];
            }
        }
        return rX;
    }

    ///@}
    ///@name Access
    ///@{

    const int* GetPointer_iL()
    {
        return iL;
    }

    const int* GetPointer_jL()
    {
        return jL;
    }

    const double* GetPointer_L()
    {
        return L;
    }

    const int* GetPointer_iU()
    {
        return iU;
    }

    const int* GetPointer_jU()
    {
        return jU;
    }

    const double* GetPointer_U()
    {
        return U;
    }

    const int GetNZL()
    {
        //return mILUSize;
        return iL[mILUSize];
    }

    const int GetNZU()
    {
        //return mILUSize;
        return iU[mILUSize];
    }

    const int GetSize()
    {
        //return mILUSize;
        return mILUSize;
    }

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Return information about this object.
    std::string Info() const override
    {
        return "ILU_LR_Preconditioner";
    }

    /// Print information about this object.
    void  PrintInfo(std::ostream& OStream) const override
    {
        OStream << "ILU_LR_Preconditioner";
    }

    void PrintData(std::ostream& OStream) const override
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

    unsigned int mILUSize;
    int *iL, *jL, *iU, *jU;
    double *L, *U;

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

}; // Class ILU_LR_Preconditioner

///@}


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

}  // namespace Kratos.

#endif // KRATOS_ILU_PRECONDITIONER_H_INCLUDED  defined
