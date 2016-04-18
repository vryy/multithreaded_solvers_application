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
//   Original Date:       $Date: 26 Aug 2014 $
//   Last Date:           $Date: 11 Jun 2015 $
//   Revision:            $Revision: 1.3 $
//
//

#if !defined(KRATOS_MULTITHREADED_SOLVERS_APP_SYSTEM_BOOST_RCM_REORDERER_PROCESS_H_INCLUDED)
#define KRATOS_MULTITHREADED_SOLVERS_APP_SYSTEM_BOOST_RCM_REORDERER_PROCESS_H_INCLUDED

// System includes
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <vector>

// External includes
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>

// Project includes
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "processes/process.h"
#include "utilities/openmp_utils.h"

namespace Kratos
{

/**
  * Calculate the permutation of the model_part based on element and condition connectivity.
  */
class SystemBoostRCMReordererProcess : public Process {

public:
    SystemBoostRCMReordererProcess(ModelPart& r_model_part) : mr_model_part(r_model_part) {}
    virtual ~SystemBoostRCMReordererProcess() {}
    
    KRATOS_CLASS_POINTER_DEFINITION(SystemBoostRCMReordererProcess);
    
    typedef typename ModelPart::ElementsContainerType ElementsContainerType;
    typedef typename ModelPart::ConditionsContainerType ConditionsContainerType;
    
    virtual void Execute();

private:
    ModelPart& mr_model_part;
};

//****************************************************************************
//std::vector<int>& SystemBoostRCMReordererProcess::CalculateIndexPermutation3(
//    ElementsContainerType& rElements,
//    ConditionsContainerType& rConditions,
//    ProcessInfo& rProcessInfo,
//    int system_size)
//{
//    double start = OpenMPUtils::GetCurrentTime();
//    
//    typedef int IndexType;
//    typedef std::size_t SizeType;
//    typedef std::multimap<SizeType, IndexType, std::greater<SizeType> > level_set_type;

//    std::vector<int> r_index_permutation;

//    r_index_permutation.resize(system_size);

//    std::vector<bool> is_marked(size, false);
//    level_set_type level_set;
//    IndexVectorType connectivity;
//    level_set_type next_level_set;

//    IndexType InitialIndex = 0;
//    r_index_permutation[0] = InitialIndex;
//    is_marked[InitialIndex] = true;

//    unsigned int next = 1;

//    level_set.insert(level_set_type::value_type(TSparseSpaceType::GraphDegree(InitialIndex, rA), InitialIndex));
//    while(next < system_size)
//    {
//        for(level_set_type::iterator i = level_set.begin() ; i != level_set.end() ; ++i)
//        {
//            TSparseSpaceType::GraphNeighbors(i->second, rA, connectivity);

//            for(IndexVectorType::iterator j = connectivity.begin() ; j != connectivity.end() ; ++j)
//            if(is_marked[*j] == false)
//            {
//                r_index_permutation[next++] = *j;
//                is_marked[*j] = true;

//                SizeType degree = TSparseSpaceType::GraphDegree(*j, rA);
//                if(degree != 0)
//                    next_level_set.insert(level_set_type::value_type(degree,*j));
//            }
//        }

//        level_set.swap(next_level_set);
//        next_level_set.clear();
//        if((level_set.size() == 0) && (next < system_size)) // No connected graph
//        {
//            SizeType k = 0;
//            while(is_marked[k]) ++k;
//            level_set.insert(level_set_type::value_type(TSparseSpaceType::GraphDegree(k, rA), k));
//            r_index_permutation[next++] = k;
//            is_marked[k] = true;
//        }
//    }
//        // for testing reverse cuthill mckee
//        //IndexVectorType temp(size);
//        //for(int i = 0 ; i < size ; i++) temp[size - i - 1] = r_index_permutation[i];
//        //r_index_permutation = temp;


//        //std::cout << "r_index_permutation : ";
//        //for(int i = 0 ; i < r_index_permutation.size() ; i++)
//        // std::cout << r_index_permutation[i] << ",";
//        //std::cout << std::endl;

//    return r_index_permutation;
//}

//****************************************************************************
/// Calculate the index permutation using Boost
void SystemBoostRCMReordererProcess::Execute()
{
    ElementsContainerType& rElements = mr_model_part.Elements();
    ConditionsContainerType& rConditions = mr_model_part.Conditions();
    ProcessInfo& rProcessInfo = mr_model_part.GetProcessInfo();
    int system_size;

    if(rProcessInfo.Has(SYSTEM_SIZE) == false)
        KRATOS_THROW_ERROR(std::logic_error, "SYSTEM_SIZE is not set", "")
    else
        system_size = rProcessInfo[SYSTEM_SIZE];

    double start = OpenMPUtils::GetCurrentTime();
    
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                                  boost::property<boost::vertex_color_t, boost::default_color_type,
                                  boost::property<boost::vertex_degree_t, int> > > Graph;
    typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef boost::graph_traits<Graph>::vertices_size_type size_type;

    // add all pairs of adjacency matrix based on element connectivity
    Graph G(10);
    std::size_t i, j;
    Element::EquationIdVectorType ids;
    for(ElementsContainerType::ptr_iterator it = rElements.ptr_begin(); it != rElements.ptr_end(); ++it)
    {
        (*it)->EquationIdVector(ids, rProcessInfo);
        for(i = 0; i < ids.size(); ++i)
            for(j = 0; j < ids.size(); ++j)
                if((i != j) && (ids[i] < system_size) && (ids[j] < system_size) && (ids[i] < ids[j]))
                    boost::add_edge(ids[i], ids[j], G);
    }
    
    // add all pairs of adjacency matrix based on condition connectivity
    for(ConditionsContainerType::ptr_iterator it = rConditions.ptr_begin(); it != rConditions.ptr_end(); ++it)
    {
        (*it)->EquationIdVector(ids, rProcessInfo);
        for(i = 0; i < ids.size(); ++i)
            for(j = 0; j < ids.size(); ++j)
                if((i != j) && (ids[i] < system_size) && (ids[j] < system_size) && (ids[i] < ids[j]))
                    boost::add_edge(ids[i], ids[j], G);
    }
    std::cout << "Reverse Cuthill-Mckee Reordering: build pairs completed, time = " << OpenMPUtils::GetCurrentTime() - start << " s" << std::endl;
    start = OpenMPUtils::GetCurrentTime();
    
    boost::graph_traits<Graph>::vertex_iterator ui, ui_end;

    boost::property_map<Graph, boost::vertex_degree_t>::type deg = boost::get(boost::vertex_degree, G);
    for (boost::tie(ui, ui_end) = vertices(G); ui != ui_end; ++ui)
        deg[*ui] = boost::degree(*ui, G);

    boost::property_map<Graph, boost::vertex_index_t>::type
    index_map = boost::get(boost::vertex_index, G);

    std::cout << "original bandwidth: " << boost::bandwidth(G) << std::endl;
    
    std::vector<Vertex> inv_perm(boost::num_vertices(G));
    std::vector<size_type> perm(boost::num_vertices(G));
  
    boost::cuthill_mckee_ordering(G, inv_perm.rbegin(), boost::get(boost::vertex_color, G), boost::make_degree_map(G));
    
    for (size_type c = 0; c != inv_perm.size(); ++c)
        perm[index_map[inv_perm[c]]] = c;
    
    std::cout << "new bandwidth: " 
              << boost::bandwidth(G, boost::make_iterator_property_map(&perm[0], index_map, perm[0]))
              << std::endl;

    vector<int> perm_vector;
    perm_vector.resize(inv_perm.size());
    std::copy(perm.begin(), perm.end(), perm_vector.begin());

    // bound check
    if(perm_vector.size() != system_size)
        KRATOS_THROW_ERROR(std::logic_error, "Error Reverse Cuthill-Mckee Reordering, perm_vector.size() != system_size", "")

    std::cout << "Reverse Cuthill-Mckee Reordering: RCM completed, time = " << OpenMPUtils::GetCurrentTime() - start << " s" << std::endl;
    rProcessInfo[SYSTEM_PERMUTATION_VECTOR] = perm_vector;
}

}

#endif

