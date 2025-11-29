//
//   Project Name:        Kratos
//   Last modified by:    $Author: hbui $
//   Date:                $Date: 11 Jun 2015 $
//   Revision:            $Revision: 1.0 $
//
//


// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "custom_processes/system_rcm_reorderer_process.h"
#include "custom_processes/system_boost_rcm_reorderer_process.h"
#ifdef MULTITHREADED_SOLVERS_APPLICATION_USE_METIS
#include "custom_processes/system_metis_reorderer_process.h"
#endif
#ifdef MULTITHREADED_SOLVERS_APPLICATION_USE_AMD
#include "custom_processes/system_amd_reorderer_process.h"
#endif
#include "add_processes_to_python.h"

namespace Kratos
{

namespace Python
{
void MultithreadedSolversApplication_AddProcessesToPython()
{
    using namespace boost::python;

    class_<SystemRCMReordererProcess, bases<Process> >("SystemRCMReordererProcess", init<ModelPart&>())
    ;

    class_<SystemBoostRCMReordererProcess, bases<Process> >("SystemBoostRCMReordererProcess", init<ModelPart&>())
    ;

    #ifdef MULTITHREADED_SOLVERS_APPLICATION_USE_METIS
    class_<SystemMetisReordererProcess, bases<Process> >("SystemMetisReordererProcess", init<ModelPart&>())
    .def("test_metis", &SystemMetisReordererProcess::test_metis)
    ;
    #endif

    #ifdef MULTITHREADED_SOLVERS_APPLICATION_USE_AMD
    class_<SystemAMDReordererProcess, bases<Process> >("SystemAMDReordererProcess", init<ModelPart&>())
    .def("test_amd", &SystemAMDReordererProcess::test_amd)
    ;
    #endif

}


}  // namespace Python.

} // Namespace Kratos

