/*
 * =======================================================================*
 * kkkk   kkkk  kkkkkkkkkk   kkkkk    kkkkkkkkkk kkkkkkkkkk kkkkkkkkkK    *
 * kkkk  kkkk   kkkk   kkkk  kkkkkk   kkkkkkkkkk kkkkkkkkkk kkkkkkkkkK    *
 * kkkkkkkkk    kkkk   kkkk  kkkkkkk     kkkk    kkk    kkk  kkkk         *
 * kkkkkkkkk    kkkkkkkkkkk  kkkk kkk    kkkk    kkk    kkk    kkkk       *
 * kkkk  kkkk   kkkk  kkkk   kkkk kkkk   kkkk    kkk    kkk      kkkk     *
 * kkkk   kkkk  kkkk   kkkk  kkkk  kkkk  kkkk    kkkkkkkkkk  kkkkkkkkkk   *
 * kkkk    kkkk kkkk    kkkk kkkk   kkkk kkkk    kkkkkkkkkk  kkkkkkkkkk   *
 *                                                                        *
 * krATos: a fREe opEN sOURce CoDE for mULti-pHysIC aDaptIVe SoLVErS,     *
 * aN extEnsIBLe OBjeCt oRiEnTEd SOlutION fOR fInITe ELemEnt fORmULatIONs *
 * Copyleft by 2003 ciMNe                                                 *
 * Copyleft by 2003 originary authors Copyleft by 2003 your name          *
 * This library is free software; you can redistribute it and/or modify   *
 * it under the terms of the GNU Lesser General Public License as         *
 * published by the Free Software Foundation; either version 2.1 of       *
 * the License, or any later version.                                     *
 *                                                                        *
 * This library is distributed in the hope that it will be useful, but    *
 * WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                   *
 * See the GNU Lesser General Public License for more details.            *
 *                                                                        *
 * You should have received a copy of the GNU Lesser General Public       *
 * License along with this library; if not, write to International Centre *
 * for Numerical Methods in Engineering (CIMNE),                          *
 * Edifici C1 - Campus Nord UPC, Gran Capità s/n, 08034 Barcelona.        *
 *                                                                        *
 * You can also contact us to the following email address:                *
 * kratos@cimne.upc.es                                                    *
 * or fax number: +34 93 401 65 17                                        *
 *                                                                        *
 * Created at Institute for Structural Mechanics                          *
 * Ruhr-University Bochum, Germany                                        *
 * Last modified by:    $Author: hbui $                                   *
 * Date:                $Date: Sep 3, 2014 $                              *
 * Revision:            $Revision: 1.0 $                                  *
 *========================================================================*
 * International Center of Numerical Methods in Engineering - CIMNE       *
 * Barcelona - Spain                                                      *
 *========================================================================*
 */

#if !defined(KRATOS_MULTITHREADED_SOLVERS_APPLICATION_ARPACK_SOLVER_H_INCLUDED )
#define  KRATOS_MULTITHREADED_SOLVERS_APPLICATION_ARPACK_SOLVER_H_INCLUDED

// External includes
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <boost/numeric/bindings/traits/sparse_traits.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_sparse.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "linear_solvers/preconditioner.h"
#include "utilities/openmp_utils.h"

namespace Kratos {

template<class TSparseSpaceType, class TDenseSpaceType>
class ArpackSolver
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ArpackSolver);

    typedef Preconditioner<TSparseSpaceType, TDenseSpaceType> PreconditionerType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    template<class TThisSparseSpaceType>
    class MatrixWithProduct
    {
    public:
        typedef typename TThisSparseSpaceType::MatrixType MatrixType;
        typedef typename TThisSparseSpaceType::VectorType VectorType;
        typedef typename TThisSparseSpaceType::DataType   DataType;
        MatrixWithProduct(MatrixType* pA) : mpA(pA) {}
        void MultMv(DataType* v, DataType* w)
        {
            int n = TThisSparseSpaceType::Size1(*mpA);
            VectorType V(n);
            for(unsigned int i = 0; i < n; ++i) V(i) = v[i];
            VectorType W(n);
            TThisSparseSpaceType::Mult(*mpA, V, W);
            for(unsigned int i = 0; i < n; ++i) w[i] = W(i);
        }
    private:
        MatrixType* mpA;
    }; // class MatrixWithProduct

    template<class TThisSparseSpaceType, class TThisDenseSpaceType>
    class InverseMatrixWithProduct
    {
    public:
        typedef typename TThisSparseSpaceType::MatrixType MatrixType;
        typedef typename TThisSparseSpaceType::VectorType VectorType;
        typedef typename TThisSparseSpaceType::DataType   DataType;
        InverseMatrixWithProduct(MatrixType* pA, typename PreconditionerType::Pointer pPrecond) : mpA(pA), mPrecond(pPrecond)
        {
            VectorType x, b;
            mPrecond->Initialize(*mpA, x, b);
        }
        void MultMv(DataType* v, DataType* w)
        {
            int n = TThisSparseSpaceType::Size1(*mpA);
            VectorType V(n);
            for(unsigned int i = 0; i < n; ++i) V(i) = v[i];
            mPrecond->ApplyLeft(V);
            for(unsigned int i = 0; i < n; ++i) w[i] = V(i);
        }
    private:
        MatrixType* mpA;
        typename PreconditionerType::Pointer mPrecond;
    }; // class MatrixWithProduct

    /**
     * Default constructor
     */
    ArpackSolver()
    {
    }

    /**
     * Constructor with solver to find smallest eigenvalue
     */
    ArpackSolver(typename PreconditionerType::Pointer pPrecond)
    {
        mPrecond = pPrecond;
    }

    /**
     * Destructor
     */
    virtual ~ArpackSolver() {}

    /**
     * Find the nlambda largest eigenvalues of A
     * @param rA. System matrix
     * @param rX. Solution vector.
     * @param rB. Right hand side vector.
     */
    bool SolveLargest(SparseMatrixType& rA, int nlambda, VectorType& rLambda);

    /**
     * Find the nlambda smallest eigenvalues of A
     * @param rA. System matrix
     * @param rX. Solution vector.
     * @param rB. Right hand side vector.
     */
    bool SolveSmallest(SparseMatrixType& rA, int nlambda, VectorType& rLambda);

    /**
     * Find the nlambda largest and nlambda smallest eigenvalues of A
     * @param rA. System matrix
     * @param rX. Solution vector.
     * @param rB. Right hand side vector.
     */
    bool Solve(SparseMatrixType& rA, int nlambda, VectorType& rLambda);

    /// Return information about this object.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "ARPACK solver";
        return buffer.str();
    }

    /**
     * Print information about this object.
     */
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ARPACK solver finished.";
    }

    /**
     * Print object's data.
     */
    void PrintData(std::ostream& rOStream) const
    {
    }

private:

    typename PreconditionerType::Pointer mPrecond;

    /**
     * Prints eigenvalues of nonsymmetric eigen-problems on standard "std::cout" stream.
     */
    template<class TProblem>
    void Solution(TProblem& Prob)
    {
        int i, n, nconv, mode;

        n     = Prob.GetN();
        nconv = Prob.ConvergedEigenvalues();
        mode  = Prob.GetMode();

        std::cout << "Real nonsymmetric eigenvalue problem: A*x - lambda*x" << std::endl;
        switch (mode) {
            case 1:
                std::cout << "Regular mode" << std::endl;
                break;
            case 3:
                std::cout << "Shift and invert mode" << std::endl;
                break;
        }

        std::cout << "Dimension of the system            : " << n              << std::endl;
        std::cout << "Number of 'requested' eigenvalues  : " << Prob.GetNev()  << std::endl;
        std::cout << "Number of 'converged' eigenvalues  : " << nconv          << std::endl;
        std::cout << "Number of Arnoldi vectors generated: " << Prob.GetNcv()  << std::endl;
        std::cout << "Number of iterations taken         : " << Prob.GetIter() << std::endl;

        if (Prob.EigenvaluesFound()) {
            // Printing eigenvalues.
            std::cout << "Eigenvalues:" << std::endl;
            for (i = 0; i < nconv; ++i) {
                std::cout << "  lambda[" << (i+1) << "]: " << Prob.EigenvalueReal(i);
                if (Prob.EigenvalueImag(i) >= 0.0) {
                    std::cout << " + " << Prob.EigenvalueImag(i) << " I" << std::endl;
                }
                else {
                    std::cout << " - " << fabs(Prob.EigenvalueImag(i)) << " I" << std::endl;
                }
            }
            std::cout << std::endl;
        }
    } // Solution

    /**
     * Assignment operator.
     */
    ArpackSolver& operator=(const ArpackSolver& Other);

    /**
     * Copy constructor.
     */

};
// Class ArpackSolver

} // namespace Kratos.

/* Remarks:
    which = 'LM' : Eigenvalues with largest magnitude (eigs, eigsh), that is, largest eigenvalues in the euclidean norm of complex numbers.
    which = 'SM' : Eigenvalues with smallest magnitude (eigs, eigsh), that is, smallest eigenvalues in the euclidean norm of complex numbers.
    which = 'LR' : Eigenvalues with largest real part (eigs)
    which = 'SR' : Eigenvalues with smallest real part (eigs)
    which = 'LI' : Eigenvalues with largest imaginary part (eigs)
    which = 'SI' : Eigenvalues with smallest imaginary part (eigs)
    which = 'LA' : Eigenvalues with largest algebraic value (eigsh), that is, largest eigenvalues inclusive of any negative sign.
    which = 'SA' : Eigenvalues with smallest algebraic value (eigsh), that is, smallest eigenvalues inclusive of any negative sign.
    which = 'BE' : Eigenvalues from both ends of the spectrum (eigsh)
 */


#endif // KRATOS_MULTITHREADED_SOLVERS_APPLICATION_ARPACK_SOLVER_H_INCLUDED  defined
