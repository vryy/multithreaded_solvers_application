/*
==============================================================================
KratosULFApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Pawel Ryzhakov
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/


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

