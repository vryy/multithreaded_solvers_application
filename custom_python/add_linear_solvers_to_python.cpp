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

#ifdef _OPENMP
#include "spaces/parallel_ublas_space.h"
#endif
#include "linear_solvers/linear_solver.h"
#include "linear_solvers/direct_solver.h"
#include "linear_solvers/iterative_solver.h"
#include "custom_linear_solvers/gmres_solver.h"
#include "custom_linear_solvers/bicgstab_solver.h"
#include "custom_linear_solvers/bicgstab_block_pressure_solver.h"
#include "custom_linear_solvers/block_pressure_staggered_solver.h"
#include "custom_linear_solvers/block_2phase_schur_solver.h"
#include "custom_linear_solvers/block_2phase_index_based_schur_solver.h"
#include "custom_linear_solvers/block_pressure_schur_solver.h"
#include "custom_linear_solvers/drained_solver.h"
#ifdef MULTITHREADED_SOLVERS_APP_USE_MKL
#include "custom_linear_solvers/bicgstab_scaling_solver.h"
#include "custom_linear_solvers/scaling_solver2.h"
#endif
#include "custom_linear_solvers/deflated_cg_solver_2.h"
#include "custom_linear_solvers/deflated_subdomain_nodal_based_cg_solver.h"
#include "custom_linear_solvers/variable_solver.h"
#include "custom_linear_solvers/diagonal_fit_solver.h"
#include "custom_linear_solvers/richardson_solver.h"
#include "custom_linear_solvers/chebyshev_solver.h"

#ifdef MULTITHREADED_SOLVERS_APPLICATION_USE_SUPERLU_MT
#include "custom_linear_solvers/superlu_mt_solver.h"
#endif

#ifdef MULTITHREADED_SOLVERS_APPLICATION_USE_UMFPACK
#include "custom_linear_solvers/umfpack_solver.h"
#endif

#ifdef MULTITHREADED_SOLVERS_APPLICATION_USE_PARDISO
#include "custom_linear_solvers/pardiso_solver.h"
#endif

#ifdef MULTITHREADED_SOLVERS_APP_USE_SPECTRA
#include "custom_eigen_solvers/spectra_eigenvalues_solver.h"
#endif

#include "custom_python/add_linear_solvers_to_python.h"

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

        typedef LinearSolver<SparseSpaceType, LocalSpaceType, ModelPart> LinearSolverType;
        typedef DirectSolver<SparseSpaceType, LocalSpaceType, ModelPart> DirectSolverType;
        typedef IterativeSolver<SparseSpaceType, LocalSpaceType, ModelPart> IterativeSolverType;
        typedef Preconditioner<SparseSpaceType, LocalSpaceType, ModelPart> PreconditionerType;

        using namespace boost::python;

        //***************************************************************************
        //linear solvers
        //***************************************************************************

        #ifdef MULTITHREADED_SOLVERS_APPLICATION_USE_UMFPACK
        typedef UmfPackSolver<SparseSpaceType, LocalSpaceType> UmfPackSolverType;
        class_<UmfPackSolverType, UmfPackSolverType::Pointer, bases<DirectSolverType>, boost::noncopyable >
        ( "UmfPackSolver", init<>() )
        .def(self_ns::str(self))
        ;
        #endif

        #ifdef MULTITHREADED_SOLVERS_APPLICATION_USE_SUPERLU_MT
        typedef SuperLUMTSolver<SparseSpaceType, LocalSpaceType> SuperLUMTSolverType;
        class_<SuperLUMTSolverType, SuperLUMTSolverType::Pointer, bases<DirectSolverType>, boost::noncopyable >
        ( "SuperLUMTSolver", init<>() )
        ;
        #endif

        #ifdef MULTITHREADED_SOLVERS_APPLICATION_USE_PARDISO
        typedef PardisoSolver<SparseSpaceType, LocalSpaceType> PardisoSolverType;
        class_<PardisoSolverType, PardisoSolverType::Pointer, bases<DirectSolverType>, boost::noncopyable >
        ( "PardisoSolver", init<>() )
        .def(init<unsigned int>() )
        .def(self_ns::str(self))
        ;
        #endif

        typedef GmresSolver<SparseSpaceType, LocalSpaceType, ModelPart> GmresSolverType;
        class_<GmresSolverType, GmresSolverType::Pointer, bases<IterativeSolverType> >
        ("GmresSolver", init<double, unsigned int, unsigned int, PreconditionerType::Pointer>())
        .def(self_ns::str(self))
        ;

        typedef BicgstabSolver<SparseSpaceType, LocalSpaceType, ModelPart> BicgstabSolverType;
        class_<BicgstabSolverType, BicgstabSolverType::Pointer, bases<IterativeSolverType> >
        ("BicgstabSolver", init<double, unsigned int, PreconditionerType::Pointer>())
        .def(init<double, unsigned int, std::string, PreconditionerType::Pointer>())
        .def(self_ns::str(self))
        ;

        typedef BicgstabBlockPressureSolver<SparseSpaceType, LocalSpaceType, ModelPart> BicgstabBlockPressureSolverType;
        class_<BicgstabBlockPressureSolverType, BicgstabBlockPressureSolverType::Pointer, bases<LinearSolverType> >
        ("BicgstabBlockPressureSolver", init<double, unsigned int, PreconditionerType::Pointer, double, double, double, double>())
        .def(self_ns::str(self))
        .def("AdditionalPhysicalDataIsNeeded", &BicgstabBlockPressureSolverType::AdditionalPhysicalDataIsNeeded)
        .def("ProvideAdditionalData", &BicgstabBlockPressureSolverType::ProvideAdditionalData)
        ;

        #ifdef MULTITHREADED_SOLVERS_APP_USE_MKL
        // TODO: to be removed (superseded by ScalingSolver2)
        typedef BicgstabScalingSolver<SparseSpaceType, LocalSpaceType, ModelPart> BicgstabScalingSolverType;
        class_<BicgstabScalingSolverType, BicgstabScalingSolverType::Pointer, bases<LinearSolverType> >
        ("BicgstabScalingSolver", init<double, unsigned int, PreconditionerType::Pointer>())
        .def(init<double, unsigned int, PreconditionerType::Pointer, int>())
        .def(self_ns::str(self))
        .def("AdditionalPhysicalDataIsNeeded", &BicgstabScalingSolverType::AdditionalPhysicalDataIsNeeded)
        .def("ProvideAdditionalData", &BicgstabScalingSolverType::ProvideAdditionalData)
        .def("DisableCheckConditionNumber", &BicgstabScalingSolverType::DisableCheckConditionNumber)
        .def("EnableCheckConditionNumber", &BicgstabScalingSolverType::EnableCheckConditionNumber)
        ;

        typedef ScalingSolver2<SparseSpaceType, LocalSpaceType, ModelPart> ScalingSolver2Type;
        class_<ScalingSolver2Type, ScalingSolver2Type::Pointer, bases<LinearSolverType> >
        ("ScalingSolver2", init<LinearSolverType::Pointer>())
        .def(init<LinearSolverType::Pointer, int>())
        .def(self_ns::str(self))
        .def("AdditionalPhysicalDataIsNeeded", &ScalingSolver2Type::AdditionalPhysicalDataIsNeeded)
        .def("ProvideAdditionalData", &ScalingSolver2Type::ProvideAdditionalData)
        .def("DisableCheckConditionNumber", &ScalingSolver2Type::DisableCheckConditionNumber)
        .def("EnableCheckConditionNumber", &ScalingSolver2Type::EnableCheckConditionNumber)
        ;
        #endif

        typedef BlockPressureStaggeredSolver<SparseSpaceType, LocalSpaceType, ModelPart> BlockPressureStaggeredSolverType;
        class_<BlockPressureStaggeredSolverType, BlockPressureStaggeredSolverType::Pointer, bases<LinearSolverType> >
        ("BlockPressureStaggeredSolver", init<LinearSolverType::Pointer, LinearSolverType::Pointer>())
        .def(self_ns::str(self))
        .def("AdditionalPhysicalDataIsNeeded", &BlockPressureStaggeredSolverType::AdditionalPhysicalDataIsNeeded)
        .def("ProvideAdditionalData", &BlockPressureStaggeredSolverType::ProvideAdditionalData)
        ;

        typedef Block2PhaseSchurSolver<SparseSpaceType, LocalSpaceType, ModelPart> Block2PhaseSchurSolverType;
        class_<Block2PhaseSchurSolverType, Block2PhaseSchurSolverType::Pointer, bases<LinearSolverType> >
        ("Block2PhaseSchurSolver", init<LinearSolverType::Pointer>())
        .def("AdditionalPhysicalDataIsNeeded", &Block2PhaseSchurSolverType::AdditionalPhysicalDataIsNeeded)
        .def("ProvideAdditionalData", &Block2PhaseSchurSolverType::ProvideAdditionalData)
        .def(self_ns::str(self))
        ;

        typedef BlockPressureSchurSolver<SparseSpaceType, LocalSpaceType, ModelPart> BlockPressureSchurSolverType;
        class_<BlockPressureSchurSolverType, BlockPressureSchurSolverType::Pointer, bases<Block2PhaseSchurSolverType> >
        ("BlockPressureSchurSolver", init<LinearSolverType::Pointer>())
        .def(self_ns::str(self))
        .def("AdditionalPhysicalDataIsNeeded", &BlockPressureSchurSolverType::AdditionalPhysicalDataIsNeeded)
        .def("ProvideAdditionalData", &BlockPressureSchurSolverType::ProvideAdditionalData)
        ;

        typedef Block2PhaseIndexBasedSchurSolver<SparseSpaceType, LocalSpaceType, ModelPart> Block2PhaseIndexBasedSchurSolverType;
        class_<Block2PhaseIndexBasedSchurSolverType, Block2PhaseIndexBasedSchurSolverType::Pointer, bases<Block2PhaseSchurSolverType> >
        ("Block2PhaseIndexBasedSchurSolver", init<LinearSolverType::Pointer>())
        .def(init<LinearSolverType::Pointer, const unsigned int&, const unsigned int&>())
        .def("AdditionalPhysicalDataIsNeeded", &Block2PhaseIndexBasedSchurSolverType::AdditionalPhysicalDataIsNeeded)
        .def("ProvideAdditionalData", &Block2PhaseIndexBasedSchurSolverType::ProvideAdditionalData)
        .def(self_ns::str(self))
        ;

        typedef DrainedSolver<SparseSpaceType, LocalSpaceType, ModelPart> DrainedSolverType;
        class_<DrainedSolverType, DrainedSolverType::Pointer, bases<LinearSolverType> >
        ("DrainedSolver", init<LinearSolverType::Pointer>())
        .def("AdditionalPhysicalDataIsNeeded", &DrainedSolverType::AdditionalPhysicalDataIsNeeded)
        .def("ProvideAdditionalData", &DrainedSolverType::ProvideAdditionalData)
        .def(self_ns::str(self))
        ;

        typedef VariableSolver<SparseSpaceType, LocalSpaceType, ModelPart> VariableSolverType;
        class_<VariableSolverType, VariableSolverType::Pointer, bases<LinearSolverType> >
        ("VariableSolver", init<>())
        .def("AddSolver", &VariableSolverType::AddSolver)
        .def("AdditionalPhysicalDataIsNeeded", &VariableSolverType::AdditionalPhysicalDataIsNeeded)
        .def("ProvideAdditionalData", &VariableSolverType::ProvideAdditionalData)
        .def(self_ns::str(self))
        ;

        typedef DiagonalFitSolver<SparseSpaceType, LocalSpaceType, ModelPart> DiagonalFitSolverType;
        class_<DiagonalFitSolverType, DiagonalFitSolverType::Pointer, bases<LinearSolverType> >
        ("DiagonalFitSolver", init<LinearSolverType::Pointer>())
        .def(init<LinearSolverType::Pointer, const double>())
        .def("AdditionalPhysicalDataIsNeeded", &DiagonalFitSolverType::AdditionalPhysicalDataIsNeeded)
        .def("ProvideAdditionalData", &DiagonalFitSolverType::ProvideAdditionalData)
        .def(self_ns::str(self))
        ;

        typedef DeflatedCGSolver2<SparseSpaceType, LocalSpaceType, ModelPart> DeflatedCGSolver2Type;
        class_<DeflatedCGSolver2Type, DeflatedCGSolver2Type::Pointer, bases<LinearSolverType> >
        ("DeflatedCGSolver2", init<double, unsigned int, PreconditionerType::Pointer>())
        .def(self_ns::str(self))
        .def("AdditionalPhysicalDataIsNeeded", &DeflatedCGSolver2Type::AdditionalPhysicalDataIsNeeded)
        .def("ProvideAdditionalData", &DeflatedCGSolver2Type::ProvideAdditionalData)
        ;

        typedef DeflatedSubdomainNodalBasedCGSolver<SparseSpaceType, LocalSpaceType, ModelPart> DeflatedSubdomainNodalBasedCGSolverType;
        class_<DeflatedSubdomainNodalBasedCGSolverType, DeflatedSubdomainNodalBasedCGSolverType::Pointer, bases<LinearSolverType> >
        ("DeflatedSubdomainNodalBasedCGSolver", init<double, unsigned int, PreconditionerType::Pointer, boost::python::list&>())
        .def(self_ns::str(self))
        .def("AdditionalPhysicalDataIsNeeded", &DeflatedSubdomainNodalBasedCGSolverType::AdditionalPhysicalDataIsNeeded)
        .def("ProvideAdditionalData", &DeflatedSubdomainNodalBasedCGSolverType::ProvideAdditionalData)
        ;

        typedef RichardsonSolver<SparseSpaceType, LocalSpaceType, ModelPart> RichardsonSolverType;
        class_<RichardsonSolverType, RichardsonSolverType::Pointer, bases<IterativeSolverType> >
        ("RichardsonSolver", init<double, double, unsigned int, PreconditionerType::Pointer>())
        .def(init<double>())
        .def(init<double, double, unsigned int>())
        .def(init<double, double, unsigned int, std::string, PreconditionerType::Pointer>())
        .def(self_ns::str(self))
        ;

        typedef ChebyshevSolver<SparseSpaceType, LocalSpaceType, ModelPart> ChebyshevSolverType;
        class_<ChebyshevSolverType, ChebyshevSolverType::Pointer, bases<IterativeSolverType> >
        ("ChebyshevSolver", init<>())
        .def(init<double, unsigned int>())
        .def(init<double, unsigned int, PreconditionerType::Pointer>())
        .def(init<double, unsigned int, PreconditionerType::Pointer, const double&, const double&>())
        .def(self_ns::str(self))
        .def("SetMinEigenvalue", &ChebyshevSolverType::SetMinEigenvalue)
        .def("SetMaxEigenvalue", &ChebyshevSolverType::SetMaxEigenvalue)
        .def("SetEstimateMinEigenvalueByGershgorin", &ChebyshevSolverType::SetEstimateMinEigenvalueByGershgorin)
        .def("SetEstimateMaxEigenvalueByGershgorin", &ChebyshevSolverType::SetEstimateMaxEigenvalueByGershgorin)
        ;

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
    }

}  // namespace Python.

} // Namespace Kratos
