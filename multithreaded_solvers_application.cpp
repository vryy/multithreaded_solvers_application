//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date:  Jan 25, 2013$
//   Revision:            $Revision: 1.1 $
//
//
//Change log:


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "multithreaded_solvers_application.h"


namespace Kratos
{

KRATOS_CREATE_VARIABLE(int, SYSTEM_SIZE )
KRATOS_CREATE_VARIABLE(boost::numeric::ublas::vector<int>, SYSTEM_PERMUTATION_VECTOR )

void KratosMultithreadedSolversApplication::Register()
{
    std::cout << "Initializing KratosMultithreadedSolversApplication..." << std::endl;

    KRATOS_REGISTER_VARIABLE(SYSTEM_SIZE)
    KRATOS_REGISTER_VARIABLE(SYSTEM_PERMUTATION_VECTOR)
}

}  // namespace Kratos.
