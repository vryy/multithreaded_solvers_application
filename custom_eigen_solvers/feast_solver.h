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

#if !defined(KRATOS_MULTITHREADED_SOLVERS_APPLICATION_FEAST_SOLVER_H_INCLUDED )
#define  KRATOS_MULTITHREADED_SOLVERS_APPLICATION_FEAST_SOLVER_H_INCLUDED

// System includes
#include <cstdio>
#include <cstdlib>
#include <cmath>

// External includes
#include <boost/numeric/bindings/traits/sparse_traits.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_sparse.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/foreach.hpp>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "utilities/openmp_utils.h"

namespace Kratos {

/**
FEAST is an eigenvalue solver that's good to find the eigenvalues in the specified interval. However, how to know that interval is unknown, maybe from practice/experience
*/
template<class TSparseSpaceType, class TDenseSpaceType>
class FeastSolver
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(FeastSolver);

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    /**
     * @param niter number of iterative refinements allowed
     */
    FeastSolver()
    {
    }

    /**
     * Destructor
     */
    virtual ~FeastSolver() {}

    /**
     * Find the nlambda eigenvalues of A in the interval
     * @param rA        System matrix
     * @param nlambda   initial guess of number of eigenvalues in [a, b]
     * @param rLambda   vector contains eigenvalues
     * @param emin      lower bound to find eigenvalues
     * @param emax      upper bound to find eigenvalues
     */
    bool Solve(SparseMatrixType& rA, int nlambda, VectorType& rLambda, double emin, double emax);

    /**
     * Find the nlambda eigenvalues of generalized eigenvalues problem AX = lambda*BX in the interval
     * @param rA        System matrix
     * @param nlambda   initial guess of number of eigenvalues in [a, b]
     * @param rLambda   vector contains eigenvalues
     * @param emin      lower bound to find eigenvalues
     * @param emax      upper bound to find eigenvalues
     */
    bool SolveGeneralized(SparseMatrixType& rA, SparseMatrixType& rB,
            int nlambda, VectorType& rLambda, boost::python::list EigenVectors,
            double emin, double emax);

    /// Return information about this object.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "FEAST solver";
        return buffer.str();
    }

    /**
     * Print information about this object.
     */
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "FEAST solver finished.";
    }

    /**
     * Print object's data.
     */
    void PrintData(std::ostream& rOStream) const
    {
    }

private:

    /**
     * Assignment operator.
     */
    FeastSolver& operator=(const FeastSolver& Other);

}; // Class FeastSolver

} // namespace Kratos.

#endif // KRATOS_MULTITHREADED_SOLVERS_APPLICATION_FEAST_SOLVER_H_INCLUDED  defined
