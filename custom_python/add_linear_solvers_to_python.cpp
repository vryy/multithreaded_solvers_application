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
#include "custom_linear_solvers/superlu_mt_solver.h"
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
#include "custom_preconditioners/solver_preconditioner.h"
#include "custom_preconditioners/ilut_preconditioner.h"
#include "custom_preconditioners/iluk_preconditioner.h"
#include "custom_preconditioners/ilu_lr_preconditioner.h"
#include "custom_preconditioners/ilu0_lr_preconditioner.h"
#include "custom_preconditioners/ilut_lr_preconditioner.h"
#include "custom_preconditioners/block_2phase_schur_preconditioner.h"
#include "custom_preconditioners/block_pressure_schur_preconditioner.h"
#include "custom_preconditioners/block_pressure_index_based_schur_preconditioner.h"
#include "custom_preconditioners/block_jacobi_pressure_preconditioner.h"
#include "custom_preconditioners/block_jacobi_preconditioner.h"
#include "custom_preconditioners/block_jacobi_nodal_based_preconditioner.h"
#include "custom_preconditioners/block_jacobi_nodal_based_pressure_preconditioner.h"
#include "custom_preconditioners/ssor_preconditioner.h"
#include "custom_preconditioners/ssor_lr_preconditioner.h"

#ifdef MULTITHREADED_SOLVERS_APPLICATION_USE_UMFPACK
#include "custom_linear_solvers/umfpack_solver.h"
#endif

#ifdef MULTITHREADED_SOLVERS_APPLICATION_USE_PARDISO
#include "custom_linear_solvers/pardiso_solver.h"
#endif

#ifdef MULTITHREADED_SOLVERS_APP_USE_SPECTRA
#include "custom_eigen_solvers/spectra_eigenvalues_solver.h"
#endif

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
        //linear solvers
        //***************************************************************************

        #ifdef MULTITHREADED_SOLVERS_APPLICATION_USE_UMFPACK
	    typedef UmfPackSolver<SparseSpaceType, LocalSpaceType> UmfPackSolverType;
	    class_<UmfPackSolverType, UmfPackSolverType::Pointer, bases<DirectSolverType>, boost::noncopyable >
	    ( "UmfPackSolver", init<>() )
        .def(self_ns::str(self))
        ;
        #endif

        typedef SuperLUMTSolver<SparseSpaceType, LocalSpaceType> SuperLUMTSolverType;
        class_<SuperLUMTSolverType, SuperLUMTSolverType::Pointer, bases<DirectSolverType>, boost::noncopyable >
        ( "SuperLUMTSolver", init<>() )
        ;

        #ifdef MULTITHREADED_SOLVERS_APPLICATION_USE_PARDISO
        typedef PardisoSolver<SparseSpaceType, LocalSpaceType> PardisoSolverType;
	    class_<PardisoSolverType, PardisoSolverType::Pointer, bases<DirectSolverType>, boost::noncopyable >
	    ( "PardisoSolver", init<>() )
        .def(init<unsigned int>() )
//        .def(self_ns::str(self))
        ;
        #endif

        typedef GmresSolver<SparseSpaceType, LocalSpaceType> GmresSolverType;
        class_<GmresSolverType, GmresSolverType::Pointer, bases<IterativeSolverType> >
        ("GmresSolver", init<double, unsigned int, unsigned int, PreconditionerType::Pointer>())
        .def(self_ns::str(self))
        ;

        typedef BicgstabSolver<SparseSpaceType, LocalSpaceType> BicgstabSolverType;
        class_<BicgstabSolverType, BicgstabSolverType::Pointer, bases<IterativeSolverType> >
        ("BicgstabSolver", init<double, unsigned int, PreconditionerType::Pointer>())
        .def(init<double, unsigned int, std::string, PreconditionerType::Pointer>())
        .def(self_ns::str(self))
        ;

        typedef BicgstabBlockPressureSolver<SparseSpaceType, LocalSpaceType> BicgstabBlockPressureSolverType;
        class_<BicgstabBlockPressureSolverType, BicgstabBlockPressureSolverType::Pointer, bases<LinearSolverType> >
        ("BicgstabBlockPressureSolver", init<double, unsigned int, PreconditionerType::Pointer, double, double, double, double>())
        .def(self_ns::str(self))
        .def("AdditionalPhysicalDataIsNeeded", &BicgstabBlockPressureSolverType::AdditionalPhysicalDataIsNeeded)
        .def("ProvideAdditionalData", &BicgstabBlockPressureSolverType::ProvideAdditionalData)
        ;

        #ifdef MULTITHREADED_SOLVERS_APP_USE_MKL
        // TODO: to be removed (superseded by ScalingSolver2)
        typedef BicgstabScalingSolver<SparseSpaceType, LocalSpaceType> BicgstabScalingSolverType;
        class_<BicgstabScalingSolverType, BicgstabScalingSolverType::Pointer, bases<LinearSolverType> >
        ("BicgstabScalingSolver", init<double, unsigned int, PreconditionerType::Pointer>())
        .def(init<double, unsigned int, PreconditionerType::Pointer, int>())
        .def(self_ns::str(self))
        .def("AdditionalPhysicalDataIsNeeded", &BicgstabScalingSolverType::AdditionalPhysicalDataIsNeeded)
        .def("ProvideAdditionalData", &BicgstabScalingSolverType::ProvideAdditionalData)
        .def("DisableCheckConditionNumber", &BicgstabScalingSolverType::DisableCheckConditionNumber)
        .def("EnableCheckConditionNumber", &BicgstabScalingSolverType::EnableCheckConditionNumber)
        ;

        typedef ScalingSolver2<SparseSpaceType, LocalSpaceType> ScalingSolver2Type;
        class_<ScalingSolver2Type, ScalingSolver2Type::Pointer, bases<LinearSolverType> >
        ("ScalingSolver2", init<LinearSolverType::Pointer>())
        .def(init<LinearSolverType::Pointer, int>())
        .def("AdditionalPhysicalDataIsNeeded", &ScalingSolver2Type::AdditionalPhysicalDataIsNeeded)
        .def("ProvideAdditionalData", &ScalingSolver2Type::ProvideAdditionalData)
        .def("DisableCheckConditionNumber", &ScalingSolver2Type::DisableCheckConditionNumber)
        .def("EnableCheckConditionNumber", &ScalingSolver2Type::EnableCheckConditionNumber)
        ;
        #endif

        typedef BlockPressureStaggeredSolver<SparseSpaceType, LocalSpaceType> BlockPressureStaggeredSolverType;
        class_<BlockPressureStaggeredSolverType, BlockPressureStaggeredSolverType::Pointer, bases<LinearSolverType> >
        ("BlockPressureStaggeredSolver", init<LinearSolverType::Pointer, LinearSolverType::Pointer>())
        .def(self_ns::str(self))
        .def("AdditionalPhysicalDataIsNeeded", &BlockPressureStaggeredSolverType::AdditionalPhysicalDataIsNeeded)
        .def("ProvideAdditionalData", &BlockPressureStaggeredSolverType::ProvideAdditionalData)
        ;

        typedef Block2PhaseSchurSolver<SparseSpaceType, LocalSpaceType> Block2PhaseSchurSolverType;
        class_<Block2PhaseSchurSolverType, Block2PhaseSchurSolverType::Pointer, bases<LinearSolverType> >
        ("Block2PhaseSchurSolver", init<LinearSolverType::Pointer>())
        .def("AdditionalPhysicalDataIsNeeded", &Block2PhaseSchurSolverType::AdditionalPhysicalDataIsNeeded)
        .def("ProvideAdditionalData", &Block2PhaseSchurSolverType::ProvideAdditionalData)
        .def(self_ns::str(self))
        ;

        typedef BlockPressureSchurSolver<SparseSpaceType, LocalSpaceType> BlockPressureSchurSolverType;
        class_<BlockPressureSchurSolverType, BlockPressureSchurSolverType::Pointer, bases<Block2PhaseSchurSolverType> >
        ("BlockPressureSchurSolver", init<LinearSolverType::Pointer>())
        .def(self_ns::str(self))
        .def("AdditionalPhysicalDataIsNeeded", &BlockPressureSchurSolverType::AdditionalPhysicalDataIsNeeded)
        .def("ProvideAdditionalData", &BlockPressureSchurSolverType::ProvideAdditionalData)
        ;

        typedef Block2PhaseIndexBasedSchurSolver<SparseSpaceType, LocalSpaceType> Block2PhaseIndexBasedSchurSolverType;
        class_<Block2PhaseIndexBasedSchurSolverType, Block2PhaseIndexBasedSchurSolverType::Pointer, bases<Block2PhaseSchurSolverType> >
        ("Block2PhaseIndexBasedSchurSolver", init<LinearSolverType::Pointer>())
        .def(init<LinearSolverType::Pointer, const unsigned int&, const unsigned int&>())
        .def("AdditionalPhysicalDataIsNeeded", &Block2PhaseIndexBasedSchurSolverType::AdditionalPhysicalDataIsNeeded)
        .def("ProvideAdditionalData", &Block2PhaseIndexBasedSchurSolverType::ProvideAdditionalData)
        .def(self_ns::str(self))
        ;

        typedef DrainedSolver<SparseSpaceType, LocalSpaceType> DrainedSolverType;
        class_<DrainedSolverType, DrainedSolverType::Pointer, bases<LinearSolverType> >
        ("DrainedSolver", init<LinearSolverType::Pointer>())
        .def("AdditionalPhysicalDataIsNeeded", &DrainedSolverType::AdditionalPhysicalDataIsNeeded)
        .def("ProvideAdditionalData", &DrainedSolverType::ProvideAdditionalData)
        .def(self_ns::str(self))
        ;

        typedef VariableSolver<SparseSpaceType, LocalSpaceType> VariableSolverType;
        class_<VariableSolverType, VariableSolverType::Pointer, bases<LinearSolverType> >
        ("VariableSolver", init<>())
        .def("AddSolver", &VariableSolverType::AddSolver)
        .def("AdditionalPhysicalDataIsNeeded", &VariableSolverType::AdditionalPhysicalDataIsNeeded)
        .def("ProvideAdditionalData", &VariableSolverType::ProvideAdditionalData)
        ;

        typedef DiagonalFitSolver<SparseSpaceType, LocalSpaceType> DiagonalFitSolverType;
        class_<DiagonalFitSolverType, DiagonalFitSolverType::Pointer, bases<LinearSolverType> >
        ("DiagonalFitSolver", init<LinearSolverType::Pointer>())
        .def(init<LinearSolverType::Pointer, const double&>())
        .def("AdditionalPhysicalDataIsNeeded", &DiagonalFitSolverType::AdditionalPhysicalDataIsNeeded)
        .def("ProvideAdditionalData", &DiagonalFitSolverType::ProvideAdditionalData)
        ;

        typedef DeflatedCGSolver2<SparseSpaceType, LocalSpaceType> DeflatedCGSolver2Type;
        class_<DeflatedCGSolver2Type, DeflatedCGSolver2Type::Pointer, bases<LinearSolverType> >
        ("DeflatedCGSolver2", init<double, unsigned int, PreconditionerType::Pointer>())
        .def(self_ns::str(self))
        .def("AdditionalPhysicalDataIsNeeded", &DeflatedCGSolver2Type::AdditionalPhysicalDataIsNeeded)
        .def("ProvideAdditionalData", &DeflatedCGSolver2Type::ProvideAdditionalData)
        ;

        typedef DeflatedSubdomainNodalBasedCGSolver<SparseSpaceType, LocalSpaceType> DeflatedSubdomainNodalBasedCGSolverType;
        class_<DeflatedSubdomainNodalBasedCGSolverType, DeflatedSubdomainNodalBasedCGSolverType::Pointer, bases<LinearSolverType> >
        ("DeflatedSubdomainNodalBasedCGSolver", init<double, unsigned int, PreconditionerType::Pointer, boost::python::list&>())
        .def(self_ns::str(self))
        .def("AdditionalPhysicalDataIsNeeded", &DeflatedSubdomainNodalBasedCGSolverType::AdditionalPhysicalDataIsNeeded)
        .def("ProvideAdditionalData", &DeflatedSubdomainNodalBasedCGSolverType::ProvideAdditionalData)
        ;

        typedef RichardsonSolver<SparseSpaceType, LocalSpaceType> RichardsonSolverType;
        class_<RichardsonSolverType, RichardsonSolverType::Pointer, bases<IterativeSolverType> >
        ("RichardsonSolver", init<double, double, unsigned int, PreconditionerType::Pointer>())
        .def(init<double>())
        .def(init<double, double, unsigned int>())
        .def(init<double, double, unsigned int, std::string, PreconditionerType::Pointer>())
        .def(self_ns::str(self))
        ;

        typedef ChebyshevSolver<SparseSpaceType, LocalSpaceType> ChebyshevSolverType;
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
        //preconditioners
        //***************************************************************************

        typedef SolverPreconditioner<SparseSpaceType, LocalSpaceType> SolverPreconditionerType;
        class_<SolverPreconditionerType, SolverPreconditionerType::Pointer, bases<PreconditionerType> >
        ("SolverPreconditioner", init<LinearSolverType::Pointer>())
        .def(self_ns::str(self))
        ;

        typedef SSORPreconditioner<SparseSpaceType, LocalSpaceType> SSORPreconditionerType;
        class_<SSORPreconditionerType, SSORPreconditionerType::Pointer, bases<PreconditionerType> >
        ("SSORPreconditioner", init<double>())
        .def(self_ns::str(self))
        ;

        typedef SSOR_LR_Preconditioner<SparseSpaceType, LocalSpaceType> SSOR_LR_PreconditionerType;
        class_<SSOR_LR_PreconditionerType, SSOR_LR_PreconditionerType::Pointer, bases<PreconditionerType> >
        ("SSOR_LR_Preconditioner", init<double>())
        .def(self_ns::str(self))
        ;

        typedef ILUtPreconditioner<SparseSpaceType, LocalSpaceType> ILUtPreconditionerType;
        class_<ILUtPreconditionerType, ILUtPreconditionerType::Pointer, bases<PreconditionerType> >
        ("ILUtPreconditioner", init<double, double>())
        .def(self_ns::str(self))
        ;

        typedef ILUkPreconditioner<SparseSpaceType, LocalSpaceType> ILUkPreconditionerType;
        class_<ILUkPreconditionerType, ILUkPreconditionerType::Pointer, bases<PreconditionerType> >
        ("ILUkPreconditioner", init<int>())
        .def(self_ns::str(self))
        ;

//        typedef ILU_LR_Preconditioner<SparseSpaceType, LocalSpaceType> ILU_LR_PreconditionerType;
//        class_<ILU_LR_PreconditionerType, ILU_LR_PreconditionerType::Pointer, bases<PreconditionerType> >
//        ("ILU_LR_Preconditioner", init<>())
//        .def(self_ns::str(self))
//        ;

        typedef ILU0_LR_Preconditioner<SparseSpaceType, LocalSpaceType> ILU0_LR_PreconditionerType;
        class_<ILU0_LR_PreconditionerType, ILU0_LR_PreconditionerType::Pointer, bases<PreconditionerType> >
        ("ILU0_LR_Preconditioner", init<>())
        .def(self_ns::str(self))
        ;

        typedef ILUt_LR_Preconditioner<SparseSpaceType, LocalSpaceType> ILUt_LR_PreconditionerType;
        class_<ILUt_LR_PreconditionerType, ILUt_LR_PreconditionerType::Pointer, bases<PreconditionerType> >
        ("ILUt_LR_Preconditioner", init<double, double>())
        .def(self_ns::str(self))
        ;

        typedef Block2PhaseSchurPreconditioner<SparseSpaceType, LocalSpaceType> Block2PhaseSchurPreconditionerType;
        class_<Block2PhaseSchurPreconditionerType, Block2PhaseSchurPreconditionerType::Pointer, bases<PreconditionerType> >
        ("Block2PhaseSchurPreconditionerType", init<PreconditionerType::Pointer, PreconditionerType::Pointer, const std::string&>())
        .def(init<PreconditionerType::Pointer, PreconditionerType::Pointer, const std::string&, LinearSolverType::Pointer>())
        .def("SetSchurMatrix", &Block2PhaseSchurPreconditionerType::SetSchurMatrix)
        .def(self_ns::str(self))
        ;

        typedef BlockPressureSchurPreconditioner<SparseSpaceType, LocalSpaceType> BlockPressureSchurPreconditionerType;
        class_<BlockPressureSchurPreconditionerType, BlockPressureSchurPreconditionerType::Pointer, bases<Block2PhaseSchurPreconditionerType> >
        ("BlockPressureSchurPreconditioner", init<PreconditionerType::Pointer, PreconditionerType::Pointer, const std::string&>())
        .def(init<PreconditionerType::Pointer, PreconditionerType::Pointer, const std::string&, const std::string&>())
        .def(init<PreconditionerType::Pointer, PreconditionerType::Pointer, const std::string&, LinearSolverType::Pointer>())
        .def(init<PreconditionerType::Pointer, PreconditionerType::Pointer, const std::string&, LinearSolverType::Pointer, const std::string&>())
        .def(self_ns::str(self))
        ;

        typedef BlockPressureIndexBasedSchurPreconditioner<SparseSpaceType, LocalSpaceType> BlockPressureIndexBasedSchurPreconditionerType;
        class_<BlockPressureIndexBasedSchurPreconditionerType, BlockPressureIndexBasedSchurPreconditionerType::Pointer, bases<Block2PhaseSchurPreconditionerType> >
        ("BlockPressureIndexBasedSchurPreconditioner", init<PreconditionerType::Pointer, PreconditionerType::Pointer, const std::string&, const std::size_t&, const std::size_t&>())
        .def(init<PreconditionerType::Pointer, PreconditionerType::Pointer, const std::string&, LinearSolverType::Pointer, const std::size_t&, const std::size_t&>())
        .def(init<PreconditionerType::Pointer, PreconditionerType::Pointer, const std::string&, const std::string&, const std::size_t&, const std::size_t&>())
        .def(init<PreconditionerType::Pointer, PreconditionerType::Pointer, const std::string&, LinearSolverType::Pointer, const std::string&, const std::size_t&, const std::size_t&>())
        .def(self_ns::str(self))
        ;

        typedef BlockJacobiPressurePreconditioner<SparseSpaceType, LocalSpaceType> BlockJacobiPressurePreconditionerType;
        class_<BlockJacobiPressurePreconditionerType, BlockJacobiPressurePreconditionerType::Pointer, bases<PreconditionerType> >
        ("BlockJacobiPressurePreconditioner", init<PreconditionerType::Pointer, PreconditionerType::Pointer>())
        .def(self_ns::str(self))
        ;

        typedef BlockJacobiPreconditioner<SparseSpaceType, LocalSpaceType> BlockJacobiPreconditionerType;
        class_<BlockJacobiPreconditionerType, BlockJacobiPreconditionerType::Pointer, bases<PreconditionerType> >
        ("BlockJacobiPreconditioner", init<>())
        .def("AddPreconditioner", &BlockJacobiPreconditionerType::AddPreconditioner)
        .def("SetPreconditioner", &BlockJacobiPreconditionerType::SetPreconditioner)
        .def(self_ns::str(self))
        ;

        typedef BlockJacobiNodalBasedPreconditioner<SparseSpaceType, LocalSpaceType> BlockJacobiNodalBasedPreconditionerType;
        class_<BlockJacobiNodalBasedPreconditionerType, BlockJacobiNodalBasedPreconditionerType::Pointer, bases<PreconditionerType> >
        ("BlockJacobiNodalBasedPreconditioner", init<boost::python::list&>())
        .def("AddPreconditioner", &BlockJacobiNodalBasedPreconditionerType::AddPreconditioner)
        .def("SetPreconditioner", &BlockJacobiNodalBasedPreconditionerType::SetPreconditioner)
        .def(self_ns::str(self))
        ;

        typedef BlockJacobiNodalBasedPressurePreconditioner<SparseSpaceType, LocalSpaceType> BlockJacobiNodalBasedPressurePreconditionerType;
        class_<BlockJacobiNodalBasedPressurePreconditionerType, BlockJacobiNodalBasedPressurePreconditionerType::Pointer, bases<PreconditionerType> >
        ("BlockJacobiNodalBasedPressurePreconditioner", init<boost::python::list&>())
        .def("AddPreconditioner", &BlockJacobiNodalBasedPressurePreconditionerType::AddPreconditioner)
        .def("SetPreconditioner", &BlockJacobiNodalBasedPressurePreconditionerType::SetPreconditioner)
        .def(self_ns::str(self))
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

