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

#if !defined(KRATOS_MULTITHREADED_SOLVERS_APP_SYSTEM_METIS_REORDERER_PROCESS_H_INCLUDED)
#define KRATOS_MULTITHREADED_SOLVERS_APP_SYSTEM_METIS_REORDERER_PROCESS_H_INCLUDED

// System includes
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <vector>

// External includes
#include "parmetis.h"

// Project includes
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "processes/process.h"
#include "utilities/openmp_utils.h"

#define PRINT_METIS_OPTIONS

extern "C"
{
    int METIS_NodeND(idx_t*, idx_t*, idx_t*, idx_t*, idx_t*, idx_t*, idx_t*);
};

namespace Kratos
{

/**
  * Calculate the permutation of the model_part based on element and condition connectivity.
  */
class SystemMetisReordererProcess : public Process {

public:
    SystemMetisReordererProcess(ModelPart& r_model_part) : mr_model_part(r_model_part) {}
    virtual ~SystemMetisReordererProcess() {}
    
    KRATOS_CLASS_POINTER_DEFINITION(SystemMetisReordererProcess);
    
    typedef typename ModelPart::ElementsContainerType ElementsContainerType;
    typedef typename ModelPart::ConditionsContainerType ConditionsContainerType;
    
    virtual void Execute();

    void test_metis();

private:
    ModelPart& mr_model_part;
};


//****************************************************************************
/// Calculate the index permutation using METIS_NodeND
void SystemMetisReordererProcess::Execute()
{
    ElementsContainerType& rElements = mr_model_part.Elements();
    ConditionsContainerType& rConditions = mr_model_part.Conditions();
    ProcessInfo& rProcessInfo = mr_model_part.GetProcessInfo();
    idx_t nvtxs; // system_size

    if(rProcessInfo.Has(SYSTEM_SIZE) == false)
        KRATOS_THROW_ERROR(std::logic_error, "SYSTEM_SIZE is not set", "")
    else
        nvtxs = rProcessInfo[SYSTEM_SIZE];

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
                if((ids[i] != ids[j]) && (ids[i] < nvtxs) && (ids[j] < nvtxs))
                {
                    AdjMap[ids[i]].insert(ids[j]);
                    AdjMap[ids[j]].insert(ids[i]);
                }
    }
    for(ConditionsContainerType::ptr_iterator it = rConditions.ptr_begin(); it != rConditions.ptr_end(); ++it)
    {
        (*it)->EquationIdVector(ids, rProcessInfo);
        for(i = 0; i < ids.size(); ++i)
            for(j = 0; j < ids.size(); ++j)
                if((ids[i] != ids[j]) && (ids[i] < nvtxs) && (ids[j] < nvtxs))
                {
                    AdjMap[ids[i]].insert(ids[j]);
                    AdjMap[ids[j]].insert(ids[i]);
                }
    }
    std::cout << "Builds adjacency map completed" << std::endl;

    if(AdjMap.size() != nvtxs)
        KRATOS_THROW_ERROR(std::logic_error, "Adjacent Map must have the same size as the number of vertices. There are something wrong in building the adjacency map", "")

    int nedges = 0;
    for(std::map<int, std::set<int> >::iterator it = AdjMap.begin(); it != AdjMap.end(); ++it)
        nedges += it->second.size();

    KRATOS_WATCH(nvtxs)
    KRATOS_WATCH(nedges)

    idx_t* xadj = new idx_t[nvtxs + 1];
    idx_t* adjncy = new idx_t[nedges];
    int cnt1 = 0, cnt2 = 0;
    xadj[0] = 0;
    for(std::map<int, std::set<int> >::iterator it = AdjMap.begin(); it != AdjMap.end(); ++it)
    {
        xadj[cnt1 + 1] = xadj[cnt1] + it->second.size();
        ++cnt1;
        for(std::set<int>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
            adjncy[cnt2++] = (*it2);
    }
    std::cout << "Converts to adjacency list completed" << std::endl;

    // METIS options
    idx_t* options = new idx_t[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options); // set default options for Metis

    options[METIS_OPTION_NUMBERING] = 0; // C-style numbering
    options[METIS_OPTION_DBGLVL] = METIS_DBG_INFO | METIS_DBG_TIME;
    options[METIS_OPTION_COMPRESS] = 0; // does not try to combine identical vertices
    options[METIS_OPTION_RTYPE] = METIS_RTYPE_SEP1SIDED;

    #ifdef PRINT_METIS_OPTIONS
    KRATOS_WATCH(options[METIS_OPTION_CTYPE])
    KRATOS_WATCH(options[METIS_OPTION_RTYPE])
    KRATOS_WATCH(options[METIS_OPTION_NO2HOP])
    KRATOS_WATCH(options[METIS_OPTION_NSEPS])
    KRATOS_WATCH(options[METIS_OPTION_NITER])
    KRATOS_WATCH(options[METIS_OPTION_UFACTOR])
    KRATOS_WATCH(options[METIS_OPTION_COMPRESS])
    KRATOS_WATCH(options[METIS_OPTION_CCORDER])
    KRATOS_WATCH(options[METIS_OPTION_SEED])
    KRATOS_WATCH(options[METIS_OPTION_PFACTOR])
    KRATOS_WATCH(options[METIS_OPTION_NUMBERING])
    KRATOS_WATCH(options[METIS_OPTION_DBGLVL])
    #endif

    // calling METIS
    idx_t* vwgt = NULL;
    idx_t* perm = new idx_t[nvtxs];
    idx_t* iperm = new idx_t[nvtxs];
    int error = METIS_NodeND(&nvtxs, xadj, adjncy, vwgt, options, perm, iperm);
    std::cout << "METIS_NodeND completed" << std::endl;

    // check for error
    if(error == METIS_ERROR_INPUT)
        KRATOS_THROW_ERROR(std::logic_error, "Error of METIS_NodeND: Input error", "")
    else if(error == METIS_ERROR_MEMORY)
        KRATOS_THROW_ERROR(std::logic_error, "Error of METIS_NodeND: cannot allocate required memory", "")
    else if(error == METIS_ERROR)
        KRATOS_THROW_ERROR(std::logic_error, "Error of METIS_NodeND", "")

//    std::cout << "perm:";
//    for(int i = 0; i < nvtxs; ++i)
//        std::cout << " " << perm[i];
//    std::cout << std::endl;
// 
//    std::cout << "iperm:";
//    for(int i = 0; i < nvtxs; ++i)
//        std::cout << " " << iperm[i];
//    std::cout << std::endl;

    // export the values
    vector<int> perm_vector;
    perm_vector.resize(nvtxs);
    std::copy(iperm, iperm + nvtxs, perm_vector.begin());
//    std::copy(perm, perm + nvtxs, perm_vector.begin());

    delete perm;
    delete iperm;
    delete xadj;
    delete adjncy;
    delete options;

    rProcessInfo[SYSTEM_PERMUTATION_VECTOR] = perm_vector;
    std::cout << "Compute system permutation vector by using METIS completed" << std::endl;
}

void SystemMetisReordererProcess::test_metis()
{
    std::cout << __FUNCTION__ << " starts" << std::endl;
    idx_t nvtxs = 15;
    idx_t xadj[] = {0, 2, 5, 8, 11, 13, 16, 20, 24, 28, 31, 33, 36, 39, 42, 44};
    idx_t adjncy[] = {1, 5, 0, 2, 6, 1, 3, 7, 2, 4, 8, 3, 9, 0, 6, 10, 1, 5, 7, 11, 2, 6, 8, 12, 3, 7, 9, 13, 4, 8, 14, 5, 11, 6, 10, 12, 7, 11, 13, 8, 12, 14, 9, 13};
    idx_t perm[nvtxs];
    idx_t iperm[nvtxs];

    METIS_NodeND(&nvtxs, xadj, adjncy, NULL, NULL, perm, iperm);

    std::cout << "perm:";
    for(int i = 0; i < nvtxs; ++i)
        std::cout << " " << perm[i];
    std::cout << std::endl;

    std::cout << "iperm:";
    for(int i = 0; i < nvtxs; ++i)
        std::cout << " " << iperm[i];
    std::cout << std::endl;
}

}

#undef PRINT_METIS_OPTIONS

#endif

