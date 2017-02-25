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

#include "spaces/ublas_space.h"
#include "arpack_solver.h"
#include "arsnsym.h"

namespace Kratos {

    template<class TSparseSpaceType, class TDenseSpaceType>
    bool ArpackSolver<TSparseSpaceType, TDenseSpaceType>::SolveLargest(ArpackSolver::SparseMatrixType& rA,
            int nlambda, ArpackSolver::VectorType& rLambda)
    {
       double start_solver = OpenMPUtils::GetCurrentTime();

       if(TSparseSpaceType::Size(rLambda) != nlambda)
           TSparseSpaceType::Resize(rLambda, nlambda);

       int m = TSparseSpaceType::Size1(rA);
       int n = TSparseSpaceType::Size2(rA);

       typedef MatrixWithProduct<TSparseSpaceType> MatrixWithProductType;
       MatrixWithProductType Awrap(&rA);
       ARNonSymStdEig<double, MatrixWithProductType> dprob(n, nlambda, &Awrap, &MatrixWithProductType::MultMv);
       dprob.FindEigenvectors();
       Solution(dprob);

       std::cout << "#### EIGENSOLVER TIME: " << OpenMPUtils::GetCurrentTime() - start_solver << " ####" << std::endl;
       return true;
    }

    template<class TSparseSpaceType, class TDenseSpaceType>
    bool ArpackSolver<TSparseSpaceType, TDenseSpaceType>::SolveSmallest(ArpackSolver::SparseMatrixType& rA,
            int nlambda, ArpackSolver::VectorType& rLambda)
    {
       double start_solver = OpenMPUtils::GetCurrentTime();

       if(TSparseSpaceType::Size(rLambda) != nlambda)
           TSparseSpaceType::Resize(rLambda, nlambda);

       int m = TSparseSpaceType::Size1(rA);
       int n = TSparseSpaceType::Size2(rA);

//        typedef MatrixWithProduct<TSparseSpaceType> MatrixWithProductType;
//        MatrixWithProductType Awrap(&rA);
//        ARNonSymStdEig<double, MatrixWithProductType> dprob(n, nlambda, &Awrap, &MatrixWithProduct<TSparseSpaceType>::MultMv, "SM");

       // this code uses inverse of the matrix to find smallest eigenvalues. It will use a preconditioner (typically a solver wrapped by preconditioner) provided from outside for solving
       typedef InverseMatrixWithProduct<TSparseSpaceType, TDenseSpaceType> InverseMatrixWithProductType;
       InverseMatrixWithProductType Awrap(&rA, mPrecond);
       ARNonSymStdEig<double, InverseMatrixWithProductType> dprob(n, nlambda, &Awrap, &InverseMatrixWithProductType::MultMv);
       dprob.FindEigenvectors();
       Solution(dprob);

       std::cout << "#### EIGENSOLVER TIME: " << OpenMPUtils::GetCurrentTime() - start_solver << " ####" << std::endl;
       return true;
    }

    template<class TSparseSpaceType, class TDenseSpaceType>
    bool ArpackSolver<TSparseSpaceType, TDenseSpaceType>::Solve(ArpackSolver::SparseMatrixType& rA,
            int nlambda, ArpackSolver::VectorType& rLambda)
    {
       double start_solver = OpenMPUtils::GetCurrentTime();

       if(TSparseSpaceType::Size(rLambda) != 2*nlambda)
           TSparseSpaceType::Resize(rLambda, 2*nlambda);

       int m = TSparseSpaceType::Size1(rA);
       int n = TSparseSpaceType::Size2(rA);

       typedef MatrixWithProduct<TSparseSpaceType> MatrixWithProductType;
       MatrixWithProductType Awrap(&rA);
       ARNonSymStdEig<double, MatrixWithProductType> dprob(n, 2*nlambda, &Awrap, &MatrixWithProduct<TSparseSpaceType>::MultMv, "BE");
       dprob.FindEigenvectors();
       Solution(dprob);

       std::cout << "#### EIGENSOLVER TIME: " << OpenMPUtils::GetCurrentTime() - start_solver << " ####" << std::endl;
       return true;
    }

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    template class ArpackSolver<SparseSpaceType, LocalSpaceType>;

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
