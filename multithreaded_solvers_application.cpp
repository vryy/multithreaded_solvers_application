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
#include "multithreaded_solvers_application.h"

#ifdef MULTITHREADED_SOLVERS_APPLICATION_USE_PARDISO
extern "C" void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
#endif

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

bool KratosMultithreadedSolversApplication::Has(const std::string& SolverName)
{
    if (SolverName == "UmfPackSolver")
    {
#ifdef MULTITHREADED_SOLVERS_APPLICATION_USE_UMFPACK
        return true;
#else
        return false;
#endif
    }
    else if (SolverName == "SuperLUMTSolver")
    {
#ifdef MULTITHREADED_SOLVERS_APPLICATION_USE_SUPERLU_MT
        return true;
#else
        return false;
#endif
    }
    else if (SolverName == "PardisoSolver")
    {
#ifdef MULTITHREADED_SOLVERS_APPLICATION_USE_PARDISO
        int mtype = 1;
        void *pt[64];
        int iparm[64];
        double dparm[64];

        int solver = 0, error = 0;
        pardisoinit (pt, &mtype, &solver, iparm, dparm, &error);
        if (error == 0)
            return true;
#endif
        return false;
    }

    return false;
}

}  // namespace Kratos.
