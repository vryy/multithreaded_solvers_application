//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         multithreaded_solvers_application/LICENSE.txt
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Hoang-Giang Bui
//  Date:            27 Nov 2018
//


#if !defined(KRATOS_SPECTRA_EIGENVALUES_SOLVER_H_INCLUDED )
#define  KRATOS_SPECTRA_EIGENVALUES_SOLVER_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes
#include <Spectra/GenEigsSolver.h>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/SymEigsShiftSolver.h>

#include "external_includes/boost/numeric/bindings/traits/sparse_traits.hpp"
#include "external_includes/boost/numeric/bindings/traits/ublas_matrix.hpp"
#include "external_includes/boost/numeric/bindings/traits/ublas_sparse.hpp"

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"


namespace Kratos
{
///@addtogroup MultithreadedSolversApplication
///@{

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

struct CompressedMatrixSpectraOp
{
    CompressedMatrixSpectraOp(CompressedMatrix& rA) : mA(rA)
    {}

    int rows() { return static_cast<int>(mA.size1()); }
    int cols() { return static_cast<int>(mA.size2()); }

    // y_out = M * x_in
    void perform_op(double *x_in, double *y_out)
    {
        typedef boost::numeric::bindings::traits::sparse_matrix_traits<CompressedMatrix> matraits;
        double* a = matraits::value_storage(mA);

        for (int r = 0; r < rows(); ++r)
        {
            int nz = mA.index1_data()[r+1] - mA.index1_data()[r];
            y_out[r] = 0.0;
            for (int j = 0; j < nz; ++j)
            {
                int ind = mA.index1_data()[r] + j;
                int c = mA.index2_data()[ind];
                y_out[r] += a[ind] * x_in[c];
            }
        }
    }

    CompressedMatrix& mA;
};

struct CompressedMatrixInversedSpectraOp
{
    typedef UblasSpace<double, CompressedMatrix, Vector> SpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SpaceType, LocalSpaceType, ModelPart> LinearSolverType;

    CompressedMatrixInversedSpectraOp(LinearSolverType::Pointer pLinearSolver, CompressedMatrix& rA)
    : mpLinearSystemSolver(pLinearSolver), mA(rA), msigma(0.0)
    {}

    int rows() { return static_cast<int>(mA.size1()); }
    int cols() { return static_cast<int>(mA.size2()); }

    void set_shift(double sigma) { msigma = sigma; }

    // y_out = M * x_in
    void perform_op(double *x_in, double *y_out)
    {
        Vector X(rows()), Y(rows());

        for (int i = 0; i < rows(); ++i)
            Y(i) = x_in[i];

        if (msigma != 0.0)
        {
            CompressedMatrix temp = mA - msigma*IdentityMatrix(rows());
            mpLinearSystemSolver->Solve(temp, X, Y);
        }
        else
            mpLinearSystemSolver->Solve(mA, X, Y);

        for (int i = 0; i < rows(); ++i)
            y_out[i] = X(i);
    }

    CompressedMatrix& mA;
    LinearSolverType::Pointer mpLinearSystemSolver;
    double msigma;
};

/// Short class definition.
/** Eigenvalues solver base on Spectra
*/
class SpectraEigenvaluesSolver
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SpectraEigenvaluesSolver
    KRATOS_CLASS_POINTER_DEFINITION(SpectraEigenvaluesSolver);

    typedef UblasSpace<double, CompressedMatrix, Vector> SpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SpaceType,  LocalSpaceType, ModelPart> LinearSolverType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SpectraEigenvaluesSolver() {}

    /// Destructor.
    virtual ~SpectraEigenvaluesSolver() {}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /// Find the largest eigenvalues of the unsymmetric matrix
    static int SolveLargestUnsym(CompressedMatrix& rA, const int& ne,
        std::vector<double>& eigenvalues_real, std::vector<double>& eigenvalues_imag, const int mem_factor = 2)
    {
        CompressedMatrixSpectraOp op(rA);
        Spectra::GenEigsSolver<double, Spectra::LARGEST_MAGN, CompressedMatrixSpectraOp> eigs(&op, ne, mem_factor*ne);

        // Initialize and compute
        eigs.init();
        int nconv = eigs.compute();

        // Retrieve results
        Eigen::VectorXcd evalues;
        if(eigs.info() == Spectra::SUCCESSFUL)
        {
            evalues = eigs.eigenvalues();
            eigenvalues_real.resize(evalues.size());
            eigenvalues_imag.resize(evalues.size());
            for (std::size_t i = 0; i < evalues.size(); ++i)
            {
                eigenvalues_real[i] = evalues(i).real();
                eigenvalues_imag[i] = evalues(i).imag();
            }
        }
        else
        {
            KRATOS_ERROR << "Spectra failed to compute largest eigenvalues, error code = " << eigs.info();
        }

        return eigs.info();
    }

    /// Find the largest eigenvalues of the symmetric matrix
    static int SolveLargestSym(CompressedMatrix& rA, const int& ne,
        std::vector<double>& eigenvalues, const int mem_factor = 2)
    {
        CompressedMatrixSpectraOp op(rA);
        Spectra::SymEigsSolver<double, Spectra::LARGEST_MAGN, CompressedMatrixSpectraOp> eigs(&op, ne, mem_factor*ne);

        // Initialize and compute
        eigs.init();
        int nconv = eigs.compute();

        // Retrieve results
        Eigen::VectorXd evalues;
        if(eigs.info() == Spectra::SUCCESSFUL)
        {
            evalues = eigs.eigenvalues();
            eigenvalues.resize(evalues.size());
            for (std::size_t i = 0; i < evalues.size(); ++i)
            {
                eigenvalues[i] = evalues(i);
            }
        }
        else
        {
            KRATOS_ERROR << "Spectra failed to compute largest eigenvalues, error code = " << eigs.info();
        }

        return eigs.info();
    }

    /// Find the smallest eigenvalues of the symmetric matrix
    /// The matrix shall be SPD so that the eigenvalues around 0.0 are determined
    static int SolveSmallestSPD(CompressedMatrix& rA, LinearSolverType::Pointer pLinearSolver,
        const int& ne, std::vector<double>& eigenvalues, const int mem_factor = 2)
    {
        CompressedMatrixInversedSpectraOp op(pLinearSolver, rA);
        Spectra::SymEigsShiftSolver<double, Spectra::LARGEST_MAGN, CompressedMatrixInversedSpectraOp> eigs(&op, ne, mem_factor*ne, 0.0);

        // Initialize and compute
        eigs.init();
        int nconv = eigs.compute();

        // Retrieve results
        Eigen::VectorXd evalues;
        if(eigs.info() == Spectra::SUCCESSFUL)
        {
            evalues = eigs.eigenvalues();
            eigenvalues.resize(evalues.size());
            for (std::size_t i = 0; i < evalues.size(); ++i)
            {
                eigenvalues[i] = evalues(i);
            }
        }
        else
        {
            KRATOS_ERROR << "Spectra failed to compute smallest eigenvalues, error code = " << eigs.info();
        }

        return eigs.info();
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
        return "Eigenvalues Solver based on SPECTRA";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
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

    /// Assignment operator.
    SpectraEigenvaluesSolver& operator=(SpectraEigenvaluesSolver const& rOther);

    /// Copy constructor.
    SpectraEigenvaluesSolver(SpectraEigenvaluesSolver const& rOther);


    ///@}

}; // Class SpectraEigenvaluesSolver

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                SpectraEigenvaluesSolver& rThis)
{}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const SpectraEigenvaluesSolver& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SPECTRA_EIGENVALUES_SOLVER_H_INCLUDED  defined
