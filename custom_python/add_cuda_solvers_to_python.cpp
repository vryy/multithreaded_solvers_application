//
//   Project Name:        Kratos
//   Last modified by:    $Author: hbui $
//   Date:                $Date: Feb 9, 2025 $
//   Revision:            $Revision: 1.0 $
//
//


// System includes

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/direct_solver.h"

#ifdef MULTITHREADED_SOLVERS_APP_USE_CUDSS
#include "custom_linear_solvers/cudss_solver.h"
#endif

#include "custom_python/add_cuda_solvers_to_python.h"

namespace Kratos
{

namespace Python
{

    void MultithreadedSolversApplication_AddCudaSolversToPython()
    {
        typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
        typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

        typedef DirectSolver<SparseSpaceType,  LocalSpaceType> DirectSolverType;

        using namespace boost::python;

        //***************************************************************************
        //linear solvers
        //***************************************************************************

        #ifdef MULTITHREADED_SOLVERS_APP_USE_CUDSS
        typedef cuDSSSolver<SparseSpaceType, LocalSpaceType> cuDSSSolverType;
        class_<cuDSSSolverType, cuDSSSolverType::Pointer, bases<DirectSolverType>, boost::noncopyable>
        ("cuDSSSolver", init<>())
        .def("AdditionalPhysicalDataIsNeeded", &cuDSSSolverType::AdditionalPhysicalDataIsNeeded)
        ;
        #endif
    }

}  // namespace Python.

} // Namespace Kratos
