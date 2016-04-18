/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

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
//   Original Date:       $Date: 11 Jun 2015 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_MULTITHREADED_SOLVERS_APP_SYSTEM_AMD_REORDERER_PROCESS_H_INCLUDED)
#define KRATOS_MULTITHREADED_SOLVERS_APP_SYSTEM_AMD_REORDERER_PROCESS_H_INCLUDED

// System includes
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <vector>

// External includes
#include "amd.h"

// Project includes
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "processes/process.h"
#include "utilities/openmp_utils.h"

// extern "C"
// {
//    int amd_order(int, int*, int*, int*, double*, double*);
// };

namespace Kratos
{

/**
  * Calculate the permutation of the model_part based on element and condition connectivity.
  */
class SystemAMDReordererProcess : public Process {

public:
    SystemAMDReordererProcess(ModelPart& r_model_part) : mr_model_part(r_model_part) {}
    virtual ~SystemAMDReordererProcess() {}
    
    KRATOS_CLASS_POINTER_DEFINITION(SystemAMDReordererProcess);
    
    typedef typename ModelPart::ElementsContainerType ElementsContainerType;
    typedef typename ModelPart::ConditionsContainerType ConditionsContainerType;
    
    virtual void Execute();

    void test_amd()
    {
        int n = 5;
        int Ap[] = {0, 2, 6, 10, 12, 14};
        int Ai[] = {0, 1, 0, 1, 2, 4, 1, 2, 3, 4, 2, 3, 1, 4};
        int P[5];
        
        amd_order(n, Ap, Ai, P, (double*)NULL, (double*)NULL);
        std::cout << "P:";
        for(int k = 0; k < n; ++k)
            std::cout << " " << P[k];
        std::cout << std::endl;
    }

private:
    ModelPart& mr_model_part;
};


//****************************************************************************
/// Calculate the index permutation using METIS_NodeND
void SystemAMDReordererProcess::Execute()
{
    ElementsContainerType& rElements = mr_model_part.Elements();
    ConditionsContainerType& rConditions = mr_model_part.Conditions();
    ProcessInfo& rProcessInfo = mr_model_part.GetProcessInfo();
    int n; // system_size

    if(rProcessInfo.Has(SYSTEM_SIZE) == false)
        KRATOS_THROW_ERROR(std::logic_error, "SYSTEM_SIZE is not set", "")
    else
        n = rProcessInfo[SYSTEM_SIZE];

    double start = OpenMPUtils::GetCurrentTime();

    // build the adjacency list
    std::map<int, std::set<int> > AdjMap;
    std::size_t i, j;
    Element::EquationIdVectorType ids;
    for(ElementsContainerType::ptr_iterator it = rElements.ptr_begin(); it != rElements.ptr_end(); ++it)
    {
        (*it)->EquationIdVector(ids, rProcessInfo);
        for(i = 0; i < ids.size(); ++i)
            for(j = 0; j < ids.size(); ++j)
                if((ids[i] < n) && (ids[j] < n))
                    AdjMap[ids[i]].insert(ids[j]);
    }
    for(ConditionsContainerType::ptr_iterator it = rConditions.ptr_begin(); it != rConditions.ptr_end(); ++it)
    {
        (*it)->EquationIdVector(ids, rProcessInfo);
        for(i = 0; i < ids.size(); ++i)
            for(j = 0; j < ids.size(); ++j)
                if((ids[i] < n) && (ids[j] < n))
                    AdjMap[ids[i]].insert(ids[j]);
    }
    std::cout << typeid(*this).name() << " builds adjacency map completed" << std::endl;

    if(AdjMap.size() != n)
        KRATOS_THROW_ERROR(std::logic_error, "Adjacent Map must have the same size as the number of vertices. There are something wrong in building the adjacency map", "")

    int nz = 0;
    for(std::map<int, std::set<int> >::iterator it = AdjMap.begin(); it != AdjMap.end(); ++it)
        nz += it->second.size();

    // build the system rowptr and colind
    int* Ap = new int[n+1];
    int* Ai = new int[nz];
    Ap[0] = 0;
    int cnt1 = 0, cnt2 = 0;
    for(std::map<int, std::set<int> >::iterator it = AdjMap.begin(); it != AdjMap.end(); ++it)
    {
        Ap[cnt1 + 1] = Ap[cnt1] + it->second.size();
        ++cnt1;
        for(std::set<int>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
            Ai[cnt2++] = (*it2);
    }

    int* P = new int[n];
    int result = amd_order(n, Ap, Ai, P, (double*)NULL, (double*)NULL);

    int* Pi = new int[n];
    for(int k = 0; k < n; ++k)
        Pi[P[k]] = k;

    // export the values
    vector<int> perm_vector;
    perm_vector.resize(n);
//    std::copy(P, P + n, perm_vector.begin());
    std::copy(Pi, Pi + n, perm_vector.begin());

    delete Ap;
    delete Ai;
    delete P;
    delete Pi;

    rProcessInfo[SYSTEM_PERMUTATION_VECTOR] = perm_vector;
    std::cout << "Compute system permutation vector by using AMD completed" << std::endl;
}

}

#endif

