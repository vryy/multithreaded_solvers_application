//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Jan 25, 2013 $
//   Revision:            $Revision: 1.1 $
//
//
//Change log:

#if !defined(KRATOS_MULTITHREADED_SOLVERS_APPLICATION_H_INCLUDED )
#define  KRATOS_MULTITHREADED_SOLVERS_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"

namespace Kratos
{

///@name Kratos Globals
///@{

// Variables definition
KRATOS_DEFINE_APPLICATION_VARIABLE(MULTITHREADED_SOLVERS_APPLICATION, int, SYSTEM_SIZE )
KRATOS_DEFINE_APPLICATION_VARIABLE(MULTITHREADED_SOLVERS_APPLICATION, boost::numeric::ublas::vector<int>, SYSTEM_PERMUTATION_VECTOR )

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(MULTITHREADED_SOLVERS_APPLICATION) KratosMultithreadedSolversApplication : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{


    /// Pointer definition of KratosMultiphaseApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosMultithreadedSolversApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosMultithreadedSolversApplication() {}

    /// Destructor.
    ~KratosMultithreadedSolversApplication() override {}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void Register() override;

    static bool Has(const std::string& SolverName);

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "KratosMultithreadedSolversApplication";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    ///// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        KRATOS_WATCH("in Kratos multithreaded solvers's application");
        KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() );
        rOStream << "Variables:" << std::endl;
        KratosComponents<VariableData>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Elements:" << std::endl;
        KratosComponents<Element>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Conditions:" << std::endl;
        KratosComponents<Condition>().PrintData(rOStream);
    }

    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    KratosMultithreadedSolversApplication& operator=(KratosMultithreadedSolversApplication const& rOther);

    /// Copy constructor.
    KratosMultithreadedSolversApplication(KratosMultithreadedSolversApplication const& rOther);

    ///@}

}; // Class KratosMultithreadedSolversApplication

///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.

#endif // KRATOS_MULTITHREADED_SOLVERS_APPLICATION_H_INCLUDED  defined
