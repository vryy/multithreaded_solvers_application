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

 	void KratosMultithreadedSolversApplication::Register()
 	{
 		// calling base class register to register Kratos components
 		KratosApplication::Register();
 		std::cout << "Initializing KratosMultithreadedSolversApplication... " << std::endl;

 	}

}  // namespace Kratos.


