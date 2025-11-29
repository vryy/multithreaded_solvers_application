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

}

}  // namespace Python.
}  // namespace Kratos.
