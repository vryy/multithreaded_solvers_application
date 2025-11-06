/*
see multithreaded_solvers_application/LICENSE.txt
*/

//
//   Project Name:        Kratos
//   Last modified by:    $Author: hbui $
//   Date:                $Date: Jan 25, 2013 $
//   Revision:            $Revision: 1.1 $
//
//


// System includes

#if defined(KRATOS_PYTHON)

// Project includes
#include "includes/define_python.h"
#include "multithreaded_solvers_application.h"
#include "custom_python/add_linear_solvers_to_python.h"
#include "custom_python/add_cuda_solvers_to_python.h"
#include "custom_python/add_eigen_solvers_to_python.h"
#include "custom_python/add_utilities_to_python.h"
#include "custom_python/add_processes_to_python.h"

namespace Kratos
{

namespace Python
{

    using namespace boost::python;

    BOOST_PYTHON_MODULE(KratosMultithreadedSolversApplication)
    {

        class_<KratosMultithreadedSolversApplication,
            KratosMultithreadedSolversApplication::Pointer,
            bases<KratosApplication>, boost::noncopyable >("KratosMultithreadedSolversApplication")
            .def("Has", &KratosMultithreadedSolversApplication::Has)
            .staticmethod("Has")
            ;

        MultithreadedSolversApplication_AddLinearSolversToPython();
        MultithreadedSolversApplication_AddCudaSolversToPython();
        MultithreadedSolversApplication_AddEigenSolversToPython();
        MultithreadedSolversApplication_AddUtilitiesToPython();
        MultithreadedSolversApplication_AddProcessesToPython();

  }

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
