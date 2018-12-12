//
//   Project Name:        Kratos
//   Last modified by:    $Author: hbui $
//   Date:                $Date: Apr 19, 2012 $
//   Revision:            $Revision: 1.0 $
//
//


// System includes

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/preconditioner.h"

#ifdef MULTITHREADED_SOLVERS_APP_USE_SPECTRA
#include "custom_eigen_solvers/spectra_eigenvalues_solver.h"
#endif

#ifdef MULTITHREADED_SOLVERS_APP_USE_ARPACK
#include "custom_eigen_solvers/arpack_solver.h"
#endif

#ifdef MULTITHREADED_SOLVERS_APP_USE_FEAST
#include "custom_eigen_solvers/feast_solver.h"
#endif

#include "multithreaded_solvers_application/multithreaded_solvers_application.h"

namespace Kratos
{

namespace Python
{

    #ifdef MULTITHREADED_SOLVERS_APP_USE_SPECTRA
    boost::python::list SpectraEigenvaluesSolver_SolveLargestUnsym(SpectraEigenvaluesSolver& rDummy, CompressedMatrix& rA, const int& ne)
    {
        boost::python::list values;

        std::vector<double> eigenvalues_real;
        std::vector<double> eigenvalues_imag;
        rDummy.SolveLargestUnsym(rA, ne, eigenvalues_real, eigenvalues_imag);

        for (std::size_t i = 0; i < eigenvalues_real.size(); ++i)
        {
            boost::python::list eval;
            eval.append(eigenvalues_real[i]);
            eval.append(eigenvalues_imag[i]);
            values.append(eval);
        }

        return values;
    }

    boost::python::list SpectraEigenvaluesSolver_SolveLargestSym(SpectraEigenvaluesSolver& rDummy, CompressedMatrix& rA, const int& ne)
    {
        boost::python::list values;

        std::vector<double> eigenvalues;
        rDummy.SolveLargestSym(rA, ne, eigenvalues);

        for (std::size_t i = 0; i < eigenvalues.size(); ++i)
        {
            values.append(eigenvalues[i]);
        }

        return values;
    }

    boost::python::list SpectraEigenvaluesSolver_SolveSmallestSPD(SpectraEigenvaluesSolver& rDummy, CompressedMatrix& rA,
        SpectraEigenvaluesSolver::LinearSolverType::Pointer pLinearSolver, const int& ne)
    {
        boost::python::list values;

        std::vector<double> eigenvalues;
        rDummy.SolveSmallestSPD(rA, pLinearSolver, ne, eigenvalues);

        for (std::size_t i = 0; i < eigenvalues.size(); ++i)
        {
            values.append(eigenvalues[i]);
        }

        return values;
    }
    #endif

    void MultithreadedSolversApplication_AddLinearSolversToPython()
    {
        typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
        typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

        typedef LinearSolver<SparseSpaceType,  LocalSpaceType> LinearSolverType;
        typedef DirectSolver<SparseSpaceType,  LocalSpaceType> DirectSolverType;
        typedef IterativeSolver<SparseSpaceType, LocalSpaceType> IterativeSolverType;
        typedef Preconditioner<SparseSpaceType,  LocalSpaceType> PreconditionerType;

        using namespace boost::python;

        //***************************************************************************
        //eigenvalue solvers
        //***************************************************************************

        #ifdef MULTITHREADED_SOLVERS_APP_USE_SPECTRA
        class_<SpectraEigenvaluesSolver, SpectraEigenvaluesSolver::Pointer, boost::noncopyable>
        ("SpectraEigenvaluesSolver", init<>())
        .def("SolveLargestUnsym", &SpectraEigenvaluesSolver_SolveLargestUnsym)
        .def("SolveLargestSym", &SpectraEigenvaluesSolver_SolveLargestSym)
        .def("SolveSmallestSPD", &SpectraEigenvaluesSolver_SolveSmallestSPD)
        ;
        #endif

        #ifdef MULTITHREADED_SOLVERS_APP_USE_ARPACK
        typedef ArpackSolver<SparseSpaceType, LocalSpaceType> ArpackSolverType;
        class_<ArpackSolverType, ArpackSolverType::Pointer, boost::noncopyable>
        ("ArpackSolver", init<>())
        .def("SolveLargest", &ArpackSolverType::SolveLargest)
        .def("SolveSmallest", &ArpackSolverType::SolveSmallest)
        .def("Solve", &ArpackSolverType::Solve)
        ;
        #endif

        #ifdef MULTITHREADED_SOLVERS_APP_USE_FEAST
        typedef FeastSolver<SparseSpaceType, LocalSpaceType> FeastSolverType;
        class_<FeastSolverType, FeastSolverType::Pointer, boost::noncopyable>
        ("FeastSolver", init<>())
        .def("Solve", &FeastSolverType::Solve)
        .def("SolveGeneralized", &FeastSolverType::SolveGeneralized)
        ;
        #endif
    }

}  // namespace Python.

} // Namespace Kratos

