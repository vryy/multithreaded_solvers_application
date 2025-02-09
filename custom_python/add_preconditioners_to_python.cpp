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
#include "custom_preconditioners/solver_preconditioner.h"
#ifdef MULTITHREADED_SOLVERS_APP_USE_FORTRAN
#include "custom_preconditioners/ilut_preconditioner.h"
#include "custom_preconditioners/ilut_lr_preconditioner.h"
#endif
#include "custom_preconditioners/iluk_preconditioner.h"
#include "custom_preconditioners/ilu_lr_preconditioner.h"
#include "custom_preconditioners/ilu0_lr_preconditioner.h"
#include "custom_preconditioners/block_2phase_schur_preconditioner.h"
#include "custom_preconditioners/block_pressure_schur_preconditioner.h"
#include "custom_preconditioners/block_pressure_index_based_schur_preconditioner.h"
#include "custom_preconditioners/block_jacobi_pressure_preconditioner.h"
#include "custom_preconditioners/block_jacobi_preconditioner.h"
#include "custom_preconditioners/block_jacobi_nodal_based_preconditioner.h"
#include "custom_preconditioners/block_jacobi_nodal_based_pressure_preconditioner.h"
#include "custom_preconditioners/ssor_preconditioner.h"
#include "custom_preconditioners/ssor_lr_preconditioner.h"

#include "custom_python/add_preconditioners_to_python.h"

namespace Kratos
{

namespace Python
{

    void MultithreadedSolversApplication_AddPreconditionersToPython()
    {
        typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
        typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

        typedef LinearSolver<SparseSpaceType,  LocalSpaceType> LinearSolverType;
        typedef Preconditioner<SparseSpaceType,  LocalSpaceType> PreconditionerType;

        using namespace boost::python;

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

        #ifdef MULTITHREADED_SOLVERS_APP_USE_FORTRAN
        typedef ILUtPreconditioner<SparseSpaceType, LocalSpaceType> ILUtPreconditionerType;
        class_<ILUtPreconditionerType, ILUtPreconditionerType::Pointer, bases<PreconditionerType> >
        ("ILUtPreconditioner", init<double, double>())
        .def(self_ns::str(self))
        ;
        #endif

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

        #ifdef MULTITHREADED_SOLVERS_APP_USE_FORTRAN
        typedef ILUt_LR_Preconditioner<SparseSpaceType, LocalSpaceType> ILUt_LR_PreconditionerType;
        class_<ILUt_LR_PreconditionerType, ILUt_LR_PreconditionerType::Pointer, bases<PreconditionerType> >
        ("ILUt_LR_Preconditioner", init<double, double>())
        .def(self_ns::str(self))
        ;
        #endif

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
    }

}  // namespace Python.

} // Namespace Kratos
