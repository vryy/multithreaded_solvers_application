/*
==============================================================================
KratosStructuralApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


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
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 22 Aug 2014 $
//   Revision:            $Revision: 1.0 $
//
//


// System includes

// External includes
#include <boost/python.hpp>
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "spaces/ublas_space.h"
#include "custom_utilities/mesh_rcm.h"

#ifdef MULTITHREADED_SOLVERS_APP_USE_ARPACK
#include "custom_utilities/arpack_solver.h"
#endif

#ifdef MULTITHREADED_SOLVERS_APP_USE_FEAST
#include "custom_utilities/feast_solver.h"
#endif

#include "multithreaded_solvers_application.h"



namespace Kratos
{

namespace Python
{

using namespace boost::python;

void  MultithreadedSolversApplication_AddUtilitiesToPython()
{

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    class_< MeshRCM, MeshRCM::Pointer, boost::noncopyable >
    ( "MeshRCM", init<>() )
    .def("Renumber", &MeshRCM::Renumber)
//    .def(self_ns::str(self))
    ;

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
}  // namespace Kratos.
