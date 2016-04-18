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

#if !defined(KRATOS_MULTITHREADED_SOLVERS_APP_SYSTEM_RCM_REORDERER_PROCESS_H_INCLUDED)
#define KRATOS_MULTITHREADED_SOLVERS_APP_SYSTEM_RCM_REORDERER_PROCESS_H_INCLUDED

// System includes
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <vector>

// External includes

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
class SystemRCMReordererProcess : public Process
{

public:
    SystemRCMReordererProcess(ModelPart& r_model_part) : mr_model_part(r_model_part) {}
    virtual ~SystemRCMReordererProcess() {}
    
    KRATOS_CLASS_POINTER_DEFINITION(SystemRCMReordererProcess);
    
    typedef typename ModelPart::ElementsContainerType ElementsContainerType;
    typedef typename ModelPart::ConditionsContainerType ConditionsContainerType;

    virtual void Execute();

private:
    ModelPart& mr_model_part;

    int adj_bandwidth ( int node_num, int adj_num, int adj_row[], int adj[] );
    int adj_perm_bandwidth ( int node_num, int adj_num, int adj_row[], int adj[], int perm[], int perm_inv[] );
    void degree ( int root, int adj_num, int adj_row[], int adj[], int mask[], int deg[], int *iccsze, int ls[], int node_num );
    int *genrcm ( int node_num, int adj_num, int adj_row[], int adj[] );
    int i4_max ( int i1, int i2 );
    int i4_min ( int i1, int i2 );
    void i4_swap ( int *i, int *j );
    int i4col_compare ( int m, int n, int a[], int i, int j );
    void i4col_sort_a ( int m, int n, int a[] );
    void i4col_sort2_a ( int m, int n, int a[] );
    int i4col_sorted_unique_count ( int m, int n, int a[] );
    void i4col_swap ( int m, int n, int a[], int icol1, int icol2 );
    void i4vec_reverse ( int n, int a[] );
    void level_set ( int root, int adj_num, int adj_row[], int adj[], int mask[], int *level_num, int level_row[], int level[], int node_num );
    bool perm_check ( int n, int p[], int base );
    int *perm_inverse3 ( int n, int perm[] );
    void r8col_permute ( int m, int n, int p[], int base, double a[] );
    void rcm ( int root, int adj_num, int adj_row[], int adj[], int mask[], int perm[], int *iccsze, int node_num );
    void root_find ( int *root, int adj_num, int adj_row[], int adj[], int mask[], int *level_num, int level_row[], int level[], int node_num );
    void sort_heap_external ( int n, int *indx, int *i, int *j, int isgn );
    void timestamp ( );
};

//****************************************************************************
//std::vector<int>& SystemRCMReordererProcess::CalculateIndexPermutation(
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
/// compute the index permutation based on internal algorithm
void SystemRCMReordererProcess::Execute()
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
    // add all pairs of adjacency matrix based on element connectivity
    std::vector<int> pair;
    int pair_num = 0;
    Element::EquationIdVectorType ids;
    for(ElementsContainerType::ptr_iterator it = rElements.ptr_begin(); it != rElements.ptr_end(); ++it)
    {
        (*it)->EquationIdVector(ids, rProcessInfo);
        for(std::size_t i = 0; i < ids.size(); ++i)
        {
            for(std::size_t j = 0; j < ids.size(); ++j)
            {
                if((i != j) && (ids[i] < system_size) && (ids[j] < system_size))
                {
                    pair.push_back(ids[i]);
                    pair.push_back(ids[j]);
                    ++pair_num;
                }
            }
        }
    }
    
    // add all pairs of adjacency matrix based on condition connectivity
    for(ConditionsContainerType::ptr_iterator it = rConditions.ptr_begin(); it != rConditions.ptr_end(); ++it)
    {
        (*it)->EquationIdVector(ids, rProcessInfo);
        for(std::size_t i = 0; i < ids.size(); ++i)
        {
            for(std::size_t j = 0; j < ids.size(); ++j)
            {
                if((i != j) && (ids[i] < system_size) && (ids[j] < system_size))
                {
                    pair.push_back(ids[i]);
                    pair.push_back(ids[j]);
                    ++pair_num;
                }
            }
        }
    }
    std::cout << "Reverse Cuthill-Mckee Reordering: build pairs completed, time = " << OpenMPUtils::GetCurrentTime() - start << " s" << std::endl;
    start = OpenMPUtils::GetCurrentTime();
    
    /////////////////part 1////////////////////////////
    // build the adjacency matrix
    int *adj;
    int adj_num;
    int *adj_row = new int[system_size + 1];
    int i, j, k, node;
    
    //
    //  Force the nodes of each pair to be listed in ascending order.
    //
    i4col_sort2_a ( 2, pair_num, &pair[0] );
    
    //
    //  Rearrange the columns in ascending order.
    //
    i4col_sort_a ( 2, pair_num, &pair[0] );

    //
    //  Get the number of unique columns.
    //
    int pair_unique_num = i4col_sorted_unique_count ( 2, pair_num, &pair[0] );

    //
    //  The number of adjacencies is TWICE this value, plus the number of nodes.
    //
    adj_num = 2 * pair_unique_num;
    
    //
    //  Now set up the ADJ_ROW counts.
    //
    for ( node = 0; node < system_size; ++node )
    {
        adj_row[node] = 0;
    }

    for ( k = 0; k < pair_num; ++k )
    {
        if ( 0 < k )
        {
            if ( pair[0+(k-1)*2] == pair[0+k*2] && pair[1+(k-1)*2] == pair[1+k*2] )
            {
                continue;
            }
        }
        i = pair[0+k*2];
        j = pair[1+k*2];

        adj_row[i] = adj_row[i] + 1;
        adj_row[j] = adj_row[j] + 1;
    }

    //
    //  We used ADJ_ROW to count the number of entries in each row.
    //  Convert it to pointers into the ADJ array.
    //
    for ( node = system_size - 1; 0 <= node; --node )
    {
        adj_row[node+1] = adj_row[node];
    }

    adj_row[0] = 0;
    for ( node = 1; node <= system_size; ++node )
    {
        adj_row[node] = adj_row[node-1] + adj_row[node];
    }
    
    /////////////////part 2////////////////////////////
    //
    //  Mark all entries of ADJ so we will know later if we missed one.
    //
    adj = new int[adj_num];

    for ( i = 0; i < adj_num; ++i )
    {
        adj[i] = -1;
    }
    
    //
    //  Copy the ADJ_ROW array and use it to keep track of the next
    //  free entry for each row.
    //
    int* adj_row_copy = new int[system_size];

    for ( node = 0; node < system_size; ++node )
    {
        adj_row_copy[node] = adj_row[node];
    }
    
    //
    //  Now set up the ADJ_ROW counts.
    //
    for ( k = 0; k < pair_num; ++k )
    {
        if ( 0 < k )
        {
            if ( pair[0+(k-1)*2] == pair[0+k*2] && pair[1+(k-1)*2] == pair[1+k*2] )
            {
                continue;
            }
        }
        i = pair[0+k*2];
        j = pair[1+k*2];

        adj[adj_row_copy[i]] = j;
        adj_row_copy[i] = adj_row_copy[i] + 1;
        adj[adj_row_copy[j]] = i;
        adj_row_copy[j] = adj_row_copy[j] + 1;
    }
    delete [] adj_row_copy;

    std::cout << "Reverse Cuthill-Mckee Reordering: build adjacency matrix completed, time = " << OpenMPUtils::GetCurrentTime() - start << " s" << std::endl;
    start = OpenMPUtils::GetCurrentTime();
    
    /////////////////part 3////////////////////////////
    //
    //  Compute the bandwidth.
    //
    double bandwidth = adj_bandwidth ( system_size, adj_num, adj_row, adj );
    std::cout << "  ADJ bandwidth = " << bandwidth << "\n";

    //
    //  GENRCM computes the RCM permutation.
    //
    int* perm = genrcm ( system_size, adj_num, adj_row, adj );
    //
    //  Compute the inverse permutation.
    //
    int* perm_inv = perm_inverse3 ( system_size, perm );

    //
    //  Compute the bandwidth of the permuted array.
    //
    bandwidth = adj_perm_bandwidth ( system_size, adj_num, adj_row, adj, perm, perm_inv );
    std::cout << "  ADJ bandwidth after RCM permutation = " << bandwidth << "\n";

    //
    // Export the data
    //
    vector<int> perm_vector;
    perm_vector.resize(system_size);
    std::copy(perm_inv, perm_inv + system_size, perm_vector.begin());

    //
    //  Free memory.
    //
    delete [] adj;
    delete [] adj_row;
    delete [] perm;
    delete [] perm_inv;
    
    //
    //  Terminate.
    //
    std::cout << "Reverse Cuthill-Mckee Reordering completed, time = " << OpenMPUtils::GetCurrentTime() - start << " s" << std::endl;
//    timestamp ( );
    rProcessInfo[SYSTEM_PERMUTATION_VECTOR] = perm_vector;
}

//****************************************************************************80

int SystemRCMReordererProcess::adj_bandwidth ( int node_num, int adj_num, int adj_row[], int adj[] )

//****************************************************************************80
//
//  Purpose:
//
//    ADJ_BANDWIDTH computes the bandwidth of an adjacency matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int ADJ_NUM, the number of adjacency entries.
//
//    Input, int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    Input, int ADJ[ADJ_NUM], the adjacency structure.
//    For each row, it contains the column indices of the nonzero entries.
//
//    Output, int ADJ_BANDWIDTH, the bandwidth of the adjacency
//    matrix.
//
{
  int band_hi;
  int band_lo;
  int col;
  int i;
  int j;
  int value;

  band_lo = 0;
  band_hi = 0;

  for ( i = 0; i < node_num; i++ )
  {
    for ( j = adj_row[i]; j <= adj_row[i+1]-1; j++ )
    {
      col = adj[j];
      band_lo = i4_max ( band_lo, i - col );
      band_hi = i4_max ( band_hi, col - i );
    }
  }

  value = band_lo + 1 + band_hi;

  return value;
}
//****************************************************************************80

int SystemRCMReordererProcess::adj_perm_bandwidth ( int node_num, int adj_num, int adj_row[], int adj[], 
  int perm[], int perm_inv[] )

//****************************************************************************80
//
//  Purpose:
//
//    ADJ_PERM_BANDWIDTH computes the bandwidth of a permuted adjacency matrix.
//
//  Discussion:
//
//    The matrix is defined by the adjacency information and a permutation.  
//
//    The routine also computes the bandwidth and the size of the envelope.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int ADJ_NUM, the number of adjacency entries.
//
//    Input, int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    Input, int ADJ[ADJ_NUM], the adjacency structure.
//    For each row, it contains the column indices of the nonzero entries.
//
//    Input, int PERM[NODE_NUM], PERM_INV(NODE_NUM), the permutation
//    and inverse permutation.
//
//    Output, int ADJ_PERM_BANDWIDTH, the bandwidth of the permuted 
//    adjacency matrix.
//
{
  int band_hi;
  int band_lo;
  int bandwidth;
  int col;
  int i;
  int j;

  band_lo = 0;
  band_hi = 0;

  for ( i = 0; i < node_num; i++ )
  {
    for ( j = adj_row[perm[i]]; j <= adj_row[perm[i]+1]-1; j++ )
    {
      col = perm_inv[adj[j]];
      band_lo = i4_max ( band_lo, i - col );
      band_hi = i4_max ( band_hi, col - i );
    }
  }

  bandwidth = band_lo + 1 + band_hi;

  return bandwidth;
}
//****************************************************************************80

void SystemRCMReordererProcess::degree ( int root, int adj_num, int adj_row[], int adj[], int mask[], 
  int deg[], int *iccsze, int ls[], int node_num )

//****************************************************************************80
//
//  Purpose:
//
//    DEGREE computes the degrees of the nodes in the connected component.
//
//  Discussion:
//
//    The connected component is specified by MASK and ROOT.
//    Nodes for which MASK is zero are ignored.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    Original FORTRAN77 version by Alan George, Joseph Liu.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//  Parameters:
//
//    Input, int ROOT, the node that defines the connected component.
//
//    Input, int ADJ_NUM, the number of adjacency entries.
//
//    Input, int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    Input, int ADJ[ADJ_NUM], the adjacency structure.
//    For each row, it contains the column indices of the nonzero entries.
//
//    Input, int MASK[NODE_NUM], is nonzero for those nodes which are
//    to be considered.
//
//    Output, int DEG[NODE_NUM], contains, for each  node in the connected
//    component, its degree.
//
//    Output, int *ICCSIZE, the number of nodes in the connected component.
//
//    Output, int LS[NODE_NUM], stores in entries 1 through ICCSIZE the nodes
//    in the connected component, starting with ROOT, and proceeding 
//    by levels.
//
//    Input, int NODE_NUM, the number of nodes.
//
{
  int i;
  int ideg;
  int j;
  int jstop;
  int jstrt;
  int lbegin;
  int lvlend;
  int lvsize;
  int nbr;
  int node;
//
//  The sign of ADJ_ROW(I) is used to indicate if node I has been considered.
//
  ls[0] = root;
  adj_row[root-1] = - adj_row[root-1];
  lvlend = 0;
  *iccsze = 1;
//
//  LBEGIN is the pointer to the beginning of the current level, and
//  LVLEND points to the end of this level.
//
  for ( ; ; )
  {
    lbegin = lvlend + 1;
    lvlend = *iccsze;
//
//  Find the degrees of nodes in the current level,
//  and at the same time, generate the next level.
//
    for ( i = lbegin; i <= lvlend; i++ )
    {
      node = ls[i-1];
      jstrt = - adj_row[node-1];
      jstop = abs ( adj_row[node] ) - 1;
      ideg = 0;

      for ( j = jstrt; j <= jstop; j++ )
      {
        nbr = adj[j-1];

        if ( mask[nbr-1] != 0 )
        {
          ideg = ideg + 1;

          if ( 0 <= adj_row[nbr-1] )
          {
            adj_row[nbr-1] = -adj_row[nbr-1];
            *iccsze = *iccsze + 1;
            ls[*iccsze-1] = nbr;
          }
        }
      }
      deg[node-1] = ideg;
    }
//
//  Compute the current level width.
//
    lvsize = *iccsze - lvlend;
//
//  If the current level width is nonzero, generate another level.
//
    if ( lvsize == 0 )
    {
      break;
    }
  }
//
//  Reset ADJ_ROW to its correct sign and return.
//
  for ( i = 0; i < *iccsze; i++ )
  {
    node = ls[i] - 1;
    adj_row[node] = -adj_row[node];
  }

  return;
}
//****************************************************************************80

int* SystemRCMReordererProcess::genrcm ( int node_num, int adj_num, int adj_row[], int adj[] )

//****************************************************************************80
//
//  Purpose:
//
//    GENRCM finds the reverse Cuthill-Mckee ordering for a general graph.
//
//  Discussion:
//
//    For each connected component in the graph, the routine obtains
//    an ordering by calling RCM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 May 2011
//
//  Author:
//
//    Original FORTRAN77 version by Alan George, Joseph Liu.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int ADJ_NUM, the number of adjacency entries.
//
//    Input, int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    Input, int  ADJ[ADJ_NUM], the adjacency structure.
//    For each row, it contains the column indices of the nonzero entries.
//
//    Output, int GENRCM[NODE_NUM], the RCM ordering.
//
//  Local Parameters:
//
//    Local, int  LEVEL_ROW[NODE_NUM+1], the index vector for a level
//    structure.  The level structure is stored in the currently unused
//    spaces in the permutation vector PERM.
//
//    Local, int MASK[NODE_NUM], marks variables that have been numbered.
//
{
  int i;
  int iccsze;
  int level_num;
  int *level_row;
  int *mask;
  int num;
  int *perm;
  int root;
//
//  Assuming the input dat is 0 based, add 1 to ADJ_ROW and ADJ,
//  because GENRCM uses 1-based indexing!
//
  for ( i = 0; i < node_num + 1; i++ )
  {
    adj_row[i] = adj_row[i] + 1;
  }
  for ( i = 0; i < adj_num; i++ )
  {
    adj[i] = adj[i] + 1;
  }

  perm = new int[node_num];
  level_row = new int[node_num+1];
  mask = new int[node_num];

  for ( i = 0; i < node_num; i++ )
  {
    mask[i] = 1;
  }

  num = 1;

  for ( i = 0; i < node_num; i++ )
  {
//
//  For each masked connected component...
//
    if ( mask[i] != 0 )
    {
      root = i + 1;
//
//  Find a pseudo-peripheral node ROOT.  The level structure found by
//  ROOT_FIND is stored starting at PERM(NUM).
//
      root_find ( &root, adj_num, adj_row, adj, mask, &level_num,
        level_row, perm+num-1, node_num );
//
//  RCM orders the component using ROOT as the starting node.
//
      rcm ( root, adj_num, adj_row, adj, mask, perm+num-1, &iccsze,
        node_num );

      num = num + iccsze;
    }
//
//  We can stop once every node is in one of the connected components.
//
    if ( node_num < num )
    {
      break;
    }
  }

  delete [] level_row;
  delete [] mask;
//
//  PERM is computed as a 1-based vector.
//  Rewrite it as a 0-based vector.
//
  for ( i = 0; i < node_num; i++ )
  {
    perm[i] = perm[i] - 1;
  }
//
//  Subtract 1 from ADJ_ROW and ADJ because GENRCM used 1-based indexing!
//
  for ( i = 0; i < node_num + 1; i++ )
  {
    adj_row[i] = adj_row[i] - 1;
  }
  for ( i = 0; i < adj_num; i++ )
  {
    adj[i] = adj[i] - 1;
  }
  return perm;
}
//****************************************************************************80

int SystemRCMReordererProcess::i4_max ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, are two integers to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

int SystemRCMReordererProcess::i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the minimum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

void SystemRCMReordererProcess::i4_swap ( int *i, int *j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_SWAP switches two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 January 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int *I, *J.  On output, the values of I and
//    J have been interchanged.
//
{
  int k;

  k = *i;
  *i = *j;
  *j = k;
 
  return;
}
//****************************************************************************80

int SystemRCMReordererProcess::i4col_compare ( int m, int n, int a[], int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4COL_COMPARE compares columns I and J of an I4COL.
//
//  Discussion:
//
//    An I4COL is an M by N array of integer values, regarded
//    as an array of N columns of length M.
//
//  Example:
//
//    Input:
//
//      M = 3, N = 4, I = 2, J = 4
//
//      A = (
//        1  2  3  4
//        5  6  7  8
//        9 10 11 12 )
//
//    Output:
//
//      I4COL_COMPARE = -1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, int A[M*N], an array of N columns of vectors of length M.
//
//    Input, int I, J, the columns to be compared.
//    I and J must be between 1 and N.
//
//    Output, int I4COL_COMPARE, the results of the comparison:
//    -1, column I < column J,
//     0, column I = column J,
//    +1, column J < column I.
//
{
  int k;
//
//  Check.
//
  if ( i < 1 )
  {
    std::cerr << "\n";
    std::cerr << "I4COL_COMPARE - Fatal error!\n";
    std::cerr << "  Column index I = " << i << " is less than 1.\n";
    exit ( 1 );
  }

  if ( n < i )
  {
    std::cerr << "\n";
    std::cerr << "I4COL_COMPARE - Fatal error!\n";
    std::cerr << "  N = " << n << " is less than column index I = " << i << ".\n";
    exit ( 1 );
  }

  if ( j < 1 )
  {
    std::cerr << "\n";
    std::cerr << "I4COL_COMPARE - Fatal error!\n";
    std::cerr << "  Column index J = " << j << " is less than 1.\n";
    exit ( 1 );
  }

  if ( n < j )
  {
    std::cerr << "\n";
    std::cerr << "I4COL_COMPARE - Fatal error!\n";
    std::cerr << "  N = " << n << " is less than column index J = " << j << ".\n";
    exit ( 1 );
  }

  if ( i == j )
  {
    return 0;
  }

  k = 1;

  while ( k <= m )
  {
    if ( a[k-1+(i-1)*m] < a[k-1+(j-1)*m] )
    {
      return (-1);
    }
    else if ( a[k-1+(j-1)*m] < a[k-1+(i-1)*m] )
    {
      return 1;
    }
    k = k + 1;
  }

  return 0;
}
//****************************************************************************80

void SystemRCMReordererProcess::i4col_sort_a ( int m, int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4COL_SORT_A ascending sorts the columns of an I4COL.
//
//  Discussion:
//
//    An I4COL is an M by N array of integer values, regarded
//    as an array of N columns of length M.
//
//    In lexicographic order, the statement "X < Y", applied to two
//    vectors X and Y of length M, means that there is some index I, with
//    1 <= I <= M, with the property that
//
//      X(J) = Y(J) for J < I,
//    and
//      X(I) < Y(I).
//
//    In other words, X is less than Y if, at the first index where they
//    differ, the X value is less than the Y value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of A.
//
//    Input, int N, the number of columns of A.
//
//    Input/output, int A[M*N].
//    On input, the array of N columns of M vectors;
//    On output, the columns of A have been sorted in ascending
//    lexicographic order.
//
{
  int i;
  int indx;
  int isgn;
  int j;
//
//  Initialize.
//
  i = 0;
  indx = 0;
  isgn = 0;
  j = 0;
//
//  Call the external heap sorter.
//
  for ( ; ; )
  {
    sort_heap_external ( n, &indx, &i, &j, isgn );
//
//  Interchange the I and J objects.
//
    if ( 0 < indx )
    {
      i4col_swap ( m, n, a, i, j );
    }
//
//  Compare the I and J objects.
//
    else if ( indx < 0 )
    {
      isgn = i4col_compare ( m, n, a, i, j );
    }
    else if ( indx == 0 )
    {
      break;
    }

  }

  return;
}
//****************************************************************************80

void SystemRCMReordererProcess::i4col_sort2_a ( int m, int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4COL_SORT2_A ascending sorts the elements of each column of an I4COL.
//
//  Discussion:
//
//    An I4COL is an M by N array of integer values, regarded
//    as an array of N columns of length M.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of A.
//
//    Input, int N, the number of columns of A, and the length
//    of a vector of data.
//
//    Input/output, int A[M*N].
//    On input, the array of N columns of M vectors.
//    On output, the elements of each column of A have been sorted in ascending
//    order.
//
{
  int col;
  int i;
  int indx;
  int isgn;
  int j;
  int row;
  int temp;

  if ( m <= 1 )
  {
    return;
  }

  if ( n <= 0 )
  {
    return;
  }
//
//  Initialize.
//
  for ( col = 0; col < n; col++ )
  {
    i = 0;
    indx = 0;
    isgn = 0;
    j = 0;
//
//  Call the external heap sorter.
//
    for ( ; ; )
    {
      sort_heap_external ( m, &indx, &i, &j, isgn );
//
//  Interchange the I and J objects.
//
      if ( 0 < indx )
      {
        temp       = a[i-1+col*m];
        a[i-1+col*m] = a[j-1+col*m];
        a[j-1+col*m] = temp;
      }
//
//  Compare the I and J objects.
//
      else if ( indx < 0 )
      {
        if ( a[j-1+col*m] < a[i-1+col*m] )
        {
          isgn = +1;
        }
        else
        {
          isgn = -1;
        }
      }
      else if ( indx == 0 )
      {
        break;
      }
    }
  }

  return;
}
//****************************************************************************80

int SystemRCMReordererProcess::i4col_sorted_unique_count ( int m, int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4COL_SORTED_UNIQUE_COUNT counts unique elements in an I4COL.
//
//  Discussion:
//
//    An I4COL is an M by N array of integer values, regarded
//    as an array of N columns of length M.
//
//    The columns of the array may be ascending or descending sorted.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, int A[M*N], a sorted array, containing
//    N columns of data.
//
//    Output, int I4COL_SORTED_UNIQUE_COUNT, the number of unique columns.
//
{
  int i;
  int j1;
  int j2;
  int unique_num;

  if ( n <= 0 )
  {
    unique_num = 0;
    return unique_num;
  }

  unique_num = 1;
  j1 = 0;

  for ( j2 = 1; j2 < n; j2++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( a[i+j1*m] != a[i+j2*m] )
      {
        unique_num = unique_num + 1;
        j1 = j2;
        break;
      }
    }
  }

  return unique_num;
}
//****************************************************************************80

void SystemRCMReordererProcess::i4col_swap ( int m, int n, int a[], int icol1, int icol2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4COL_SWAP swaps two columns of an I4COL.
//
//  Discussion:
//
//    An I4COL is an M by N array of integer values, regarded
//    as an array of N columns of length M.
//
//    The two dimensional information is stored as a one dimensional
//    array, by columns.
//
//    The row indices are 1 based, NOT 0 based!  However, a preprocessor
//    variable, called OFFSET, can be reset from 1 to 0 if you wish to
//    use 0-based indices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 April 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, int A[M*N], an array of data.
//
//    Input, int ICOL1, ICOL2, the two columns to swap.
//    These indices should be between 1 and N.
//
{
# define OFFSET 1

  int i;
  int t;
//
//  Check.
//
  if ( icol1 - OFFSET < 0 || n-1 < icol1 - OFFSET )
  {
    std::cerr << "\n";
    std::cerr << "I4COL_SWAP - Fatal error!\n";
    std::cerr << "  ICOL1 is out of range.\n";
    exit ( 1 );
  }

  if ( icol2 - OFFSET < 0 || n-1 < icol2 - OFFSET )
  {
    std::cerr << "\n";
    std::cerr << "I4COL_SWAP - Fatal error!\n";
    std::cerr << "  ICOL2 is out of range.\n";
    exit ( 1 );
  }

  if ( icol1 == icol2 )
  {
    return;
  }
  for ( i = 0; i < m; i++ )
  {
    t                     = a[i+(icol1-OFFSET)*m];
    a[i+(icol1-OFFSET)*m] = a[i+(icol2-OFFSET)*m];
    a[i+(icol2-OFFSET)*m] = t;
  }

  return;
# undef OFFSET
}
//****************************************************************************80

void SystemRCMReordererProcess::i4vec_reverse ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_REVERSE reverses the elements of an I4VEC.
//
//  Example:
//
//    Input:
//
//      N = 5,
//      A = ( 11, 12, 13, 14, 15 ).
//
//    Output:
//
//      A = ( 15, 14, 13, 12, 11 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input/output, int A(N), the array to be reversed.
//
{
  int i;
  int j;

  for ( i = 0; i < n / 2; i++ )
  {
    j        = a[i];
    a[i]     = a[n-1-i];
    a[n-1-i] = j;
  }

  return;
}
//****************************************************************************80

void SystemRCMReordererProcess::level_set ( int root, int adj_num, int adj_row[], int adj[], int mask[], 
  int *level_num, int level_row[], int level[], int node_num )

//****************************************************************************80
//
//  Purpose:
//
//    LEVEL_SET generates the connected level structure rooted at a given node.
//
//  Discussion:
//
//    Only nodes for which MASK is nonzero will be considered.
//
//    The root node chosen by the user is assigned level 1, and masked.
//    All (unmasked) nodes reachable from a node in level 1 are
//    assigned level 2 and masked.  The process continues until there
//    are no unmasked nodes adjacent to any node in the current level.
//    The number of levels may vary between 2 and NODE_NUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    Original FORTRAN77 version by Alan George, Joseph Liu.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//  Parameters:
//
//    Input, int ROOT, the node at which the level structure
//    is to be rooted.
//
//    Input, int ADJ_NUM, the number of adjacency entries.
//
//    Input, int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    Input, int ADJ[ADJ_NUM], the adjacency structure.
//    For each row, it contains the column indices of the nonzero entries.
//
//    Input/output, int MASK[NODE_NUM].  On input, only nodes with nonzero
//    MASK are to be processed.  On output, those nodes which were included
//    in the level set have MASK set to 1.
//
//    Output, int *LEVEL_NUM, the number of levels in the level
//    structure.  ROOT is in level 1.  The neighbors of ROOT
//    are in level 2, and so on.
//
//    Output, int LEVEL_ROW[NODE_NUM+1], LEVEL[NODE_NUM], the rooted 
//    level structure.
//
//    Input, int NODE_NUM, the number of nodes.
//
{
  int i;
  int iccsze;
  int j;
  int jstop;
  int jstrt;
  int lbegin;
  int lvlend;
  int lvsize;
  int nbr;
  int node;

  mask[root-1] = 0;
  level[0] = root;
  *level_num = 0;
  lvlend = 0;
  iccsze = 1;
//
//  LBEGIN is the pointer to the beginning of the current level, and
//  LVLEND points to the end of this level.
//
  for ( ; ; )
  {
    lbegin = lvlend + 1;
    lvlend = iccsze;
    *level_num = *level_num + 1;
    level_row[*level_num-1] = lbegin;
//
//  Generate the next level by finding all the masked neighbors of nodes
//  in the current level.
//
    for ( i = lbegin; i <= lvlend; i++ )
    {
      node = level[i-1];
      jstrt = adj_row[node-1];
      jstop = adj_row[node] - 1;

      for ( j = jstrt; j <= jstop; j++ )
      {
        nbr = adj[j-1];

        if ( mask[nbr-1] != 0 )
        {
          iccsze = iccsze + 1;
          level[iccsze-1] = nbr;
          mask[nbr-1] = 0;
        }
      }
    }
//
//  Compute the current level width (the number of nodes encountered.)
//  If it is positive, generate the next level.
//
    lvsize = iccsze - lvlend;

    if ( lvsize <= 0 )
    {
      break;
    }
  }

  level_row[*level_num] = lvlend + 1;
//
//  Reset MASK to 1 for the nodes in the level structure.
//
  for ( i = 0; i < iccsze; i++ )
  {
    mask[level[i]-1] = 1;
  }

  return;
}
//****************************************************************************80

bool SystemRCMReordererProcess::perm_check ( int n, int p[], int base )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_CHECK checks that a vector represents a permutation.
//
//  Discussion:
//
//    The routine verifies that each of the integers from BASE to
//    to BASE+N-1 occurs among the N entries of the permutation.
//
//    Set the input quantity BASE to 0, if P is a 0-based permutation,
//    or to 1 if P is a 1-based permutation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries.
//
//    Input, int P[N], the array to check.
//
//    Input, int BASE, the index base.
//
//    Output, bool PERM_CHECK, is TRUE if the permutation is OK.
//
{
  bool found;
  int i;
  int seek;

  for ( seek = base; seek < base + n; seek++ )
  {
    found = false;

    for ( i = 0; i < n; i++ )
    {
      if ( p[i] == seek )
      {
        found = true;
        break;
      }
    }

    if ( !found )
    {
      std::cerr << "\n";
      std::cerr << "PERM_CHECK - Fatal error!\n";
      std::cerr << "  Did not find " << found << "\n";
      exit ( 1 );
    }

  }

  return true;
}
//****************************************************************************80

int* SystemRCMReordererProcess::perm_inverse3 ( int n, int perm[] )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_INVERSE3 produces the inverse of a given permutation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of items permuted.
//
//    Input, int PERM[N], a permutation.
//
//    Output, int PERM_INVERSE3[N], the inverse permutation.
//
{
  int i;
  int *perm_inv;

  perm_inv = new int[n];

  for ( i = 0; i < n; i++ )
  {
    perm_inv[perm[i]] = i;
  }

  return perm_inv;
}
//****************************************************************************80

void SystemRCMReordererProcess::r8col_permute ( int m, int n, int p[], int base, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8COL_PERMUTE permutes an R8COL in place.
//
//  Discussion:
//
//    An R8COL is an M by N array of R8's, regarded as an array of N columns,
//    each of length M.
//
//    This routine permutes an array of real "objects", but the same
//    logic can be used to permute an array of objects of any arithmetic
//    type, or an array of objects of any complexity.  The only temporary
//    storage required is enough to store a single object.  The number
//    of data movements made is N + the number of cycles of order 2 or more,
//    which is never more than N + N/2.
//
//  Example:
//
//    Input:
//
//      M = 2
//      N = 5
//      P = (   2,    4,    5,    1,    3 )
//      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
//          (11.0, 22.0, 33.0, 44.0, 55.0 )
//      BASE = 1
//
//    Output:
//
//      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
//             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the length of objects.
//
//    Input, int N, the number of objects.
//
//    Input, int P[N], the permutation.  P(I) = J means
//    that the I-th element of the output array should be the J-th
//    element of the input array.
//
//    Input, int BASE, is 0 for a 0-based permutation and 1 for a
//    1-based permutation.
//
//    Input/output, double A[M*N], the array to be permuted.
//
{
  double *a_temp;
  int i;
  int iget;
  int iput;
  int istart;
  int j;

  if ( !perm_check ( n, p, base ) )
  {
    std::cerr << "\n";
    std::cerr << "R8COL_PERMUTE - Fatal error!\n";
    std::cerr << "  PERM_CHECK rejects this permutation.\n";
    exit ( 1 );
  }
//
//  In order for the sign negation trick to work, we need to assume that the
//  entries of P are strictly positive.  Presumably, the lowest number is BASE.
//  So temporarily add 1-BASE to each entry to force positivity.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] + 1 - base;
  }

  a_temp = new double[m];
//
//  Search for the next element of the permutation that has not been used.
//
  for ( istart = 1; istart <= n; istart++ )
  {
    if ( p[istart-1] < 0 )
    {
      continue;
    }
    else if ( p[istart-1] == istart )
    {
      p[istart-1] = - p[istart-1];
      continue;
    }
    else
    {
      for ( i = 0; i < m; i++ )
      {
        a_temp[i] = a[i+(istart-1)*m];
      }
      iget = istart;
//
//  Copy the new value into the vacated entry.
//
      for ( ; ; )
      {
        iput = iget;
        iget = p[iget-1];

        p[iput-1] = - p[iput-1];

        if ( iget < 1 || n < iget )
        {
          std::cerr << "\n";
          std::cerr << "R8COL_PERMUTE - Fatal error!\n";
          std::cerr << "  Entry IPUT = " << iput << " of the permutation has\n";
          std::cerr << "  an illegal value IGET = " << iget << ".\n";
          exit ( 1 );
        }

        if ( iget == istart )
        {
          for ( i = 0; i < m; i++ )
          {
            a[i+(iput-1)*m] = a_temp[i];
          }
          break;
        }
        for ( i = 0; i < m; i++ )
        {
          a[i+(iput-1)*m] = a[i+(iget-1)*m];
        }
      }
    }
  }
//
//  Restore the signs of the entries.
//
  for ( j = 0; j < n; j++ )
  {
    p[j] = - p[j];
  }
//
//  Restore the base of the entries.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] - 1 +  base;
  }

  delete [] a_temp;

  return;
}
//****************************************************************************80

void SystemRCMReordererProcess::rcm ( int root, int adj_num, int adj_row[], int adj[], int mask[], 
  int perm[], int *iccsze, int node_num )

//****************************************************************************80
//
//  Purpose:
//
//    RCM renumbers a connected component by the reverse Cuthill McKee algorithm.
//
//  Discussion:
//
//    The connected component is specified by a node ROOT and a mask.
//    The numbering starts at the root node.
//
//    An outline of the algorithm is as follows:
//
//    X(1) = ROOT.
//
//    for ( I = 1 to N-1)
//      Find all unlabeled neighbors of X(I),
//      assign them the next available labels, in order of increasing degree.
//
//    When done, reverse the ordering.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    Original FORTRAN77 version by Alan George, Joseph Liu.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//  Parameters:
//
//    Input, int ROOT, the node that defines the connected component.
//    It is used as the starting point for the RCM ordering.
//
//    Input, int ADJ_NUM, the number of adjacency entries.
//
//    Input, int ADJ_ROW(NODE_NUM+1).  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    Input, int ADJ(ADJ_NUM), the adjacency structure.
//    For each row, it contains the column indices of the nonzero entries.
//
//    Input/output, int MASK(NODE_NUM), a mask for the nodes.  Only 
//    those nodes with nonzero input mask values are considered by the 
//    routine.  The nodes numbered by RCM will have their mask values 
//    set to zero.
//
//    Output, int PERM(NODE_NUM), the RCM ordering.
//
//    Output, int ICCSZE, the size of the connected component
//    that has been numbered.
//
//    Input, int NODE_NUM, the number of nodes.
//
//  Local Parameters:
//
//    Workspace, int DEG[NODE_NUM], a temporary vector used to hold 
//    the degree of the nodes in the section graph specified by mask and root.
//
{
  int *deg;
  int fnbr;
  int i;
  int j;
  int jstop;
  int jstrt;
  int k;
  int l;
  int lbegin;
  int lnbr;
  int lperm;
  int lvlend;
  int nbr;
  int node;
//
//  Find the degrees of the nodes in the component specified by MASK and ROOT.
//
  deg = new int[node_num];

  degree ( root, adj_num, adj_row, adj, mask, deg, iccsze, perm, node_num );

  mask[root-1] = 0;

  if ( *iccsze <= 1 )
  {
    delete [] deg;
    return;
  }

  lvlend = 0;
  lnbr = 1;
//
//  LBEGIN and LVLEND point to the beginning and
//  the end of the current level respectively.
//
  while ( lvlend < lnbr )
  {
    lbegin = lvlend + 1;
    lvlend = lnbr;

    for ( i = lbegin; i <= lvlend; i++ )
    {
//
//  For each node in the current level...
//
      node = perm[i-1];
      jstrt = adj_row[node-1];
      jstop = adj_row[node] - 1;
//
//  Find the unnumbered neighbors of NODE.
//
//  FNBR and LNBR point to the first and last neighbors
//  of the current node in PERM.
//
      fnbr = lnbr + 1;

      for ( j = jstrt; j <= jstop; j++ )
      {
        nbr = adj[j-1];

        if ( mask[nbr-1] != 0 )
        {
          lnbr = lnbr + 1;
          mask[nbr-1] = 0;
          perm[lnbr-1] = nbr;
        }
      }
//
//  If no neighbors, skip to next node in this level.
//
      if ( lnbr <= fnbr )
      {
        continue;
      }
//
//  Sort the neighbors of NODE in increasing order by degree.
//  Linear insertion is used.
//
      k = fnbr;

      while ( k < lnbr )
      {
        l = k;
        k = k + 1;
        nbr = perm[k-1];

        while ( fnbr < l )
        {
          lperm = perm[l-1];

          if ( deg[lperm-1] <= deg[nbr-1] )
          {
            break;
          }

          perm[l] = lperm;
          l = l - 1;
        }
        perm[l] = nbr;
      }
    }
  }
//
//  We now have the Cuthill-McKee ordering.  Reverse it.
//
  i4vec_reverse ( *iccsze, perm );

  delete [] deg;

  return;
}
//****************************************************************************80

void SystemRCMReordererProcess::root_find ( int *root, int adj_num, int adj_row[], int adj[], int mask[], 
  int *level_num, int level_row[], int level[], int node_num )

//****************************************************************************80
//
//  Purpose:
//
//    ROOT_FIND finds a pseudo-peripheral node.
//
//  Discussion:
//
//    The diameter of a graph is the maximum distance (number of edges)
//    between any two nodes of the graph.
//
//    The eccentricity of a node is the maximum distance between that
//    node and any other node of the graph.
//
//    A peripheral node is a node whose eccentricity equals the
//    diameter of the graph.
//
//    A pseudo-peripheral node is an approximation to a peripheral node;
//    it may be a peripheral node, but all we know is that we tried our
//    best.
//
//    The routine is given a graph, and seeks pseudo-peripheral nodes,
//    using a modified version of the scheme of Gibbs, Poole and
//    Stockmeyer.  It determines such a node for the section subgraph
//    specified by MASK and ROOT.
//
//    The routine also determines the level structure associated with
//    the given pseudo-peripheral node; that is, how far each node
//    is from the pseudo-peripheral node.  The level structure is
//    returned as a list of nodes LS, and pointers to the beginning
//    of the list of nodes that are at a distance of 0, 1, 2, ...,
//    NODE_NUM-1 from the pseudo-peripheral node.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    Original FORTRAN77 version by Alan George, Joseph Liu.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//    Norman Gibbs, William Poole, Paul Stockmeyer,
//    An Algorithm for Reducing the Bandwidth and Profile of a Sparse Matrix,
//    SIAM Journal on Numerical Analysis,
//    Volume 13, pages 236-250, 1976.
//
//    Norman Gibbs,
//    Algorithm 509: A Hybrid Profile Reduction Algorithm,
//    ACM Transactions on Mathematical Software,
//    Volume 2, pages 378-387, 1976.
//
//  Parameters:
//
//    Input/output, int *ROOT.  On input, ROOT is a node in the
//    the component of the graph for which a pseudo-peripheral node is
//    sought.  On output, ROOT is the pseudo-peripheral node obtained.
//
//    Input, int ADJ_NUM, the number of adjacency entries.
//
//    Input, int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    Input, int ADJ[ADJ_NUM], the adjacency structure.
//    For each row, it contains the column indices of the nonzero entries.
//
//    Input, int MASK[NODE_NUM], specifies a section subgraph.  Nodes 
//    for which MASK is zero are ignored by FNROOT.
//
//    Output, int *LEVEL_NUM, is the number of levels in the level structure
//    rooted at the node ROOT.
//
//    Output, int LEVEL_ROW(NODE_NUM+1), LEVEL(NODE_NUM), the 
//    level structure array pair containing the level structure found.
//
//    Input, int NODE_NUM, the number of nodes.
//
{
  int iccsze;
  int j;
  int jstrt;
  int k;
  int kstop;
  int kstrt;
  int level_num2;
  int mindeg;
  int nabor;
  int ndeg;
  int node;
//
//  Determine the level structure rooted at ROOT.
//
  level_set ( *root, adj_num, adj_row, adj, mask, level_num, 
    level_row, level, node_num );
//
//  Count the number of nodes in this level structure.
//
  iccsze = level_row[*level_num] - 1;
//
//  Extreme case:
//    A complete graph has a level set of only a single level.
//    Every node is equally good (or bad).
//
  if ( *level_num == 1 )
  {
    return;
  }
//
//  Extreme case:
//    A "line graph" 0--0--0--0--0 has every node in its only level.
//    By chance, we've stumbled on the ideal root.
//
  if ( *level_num == iccsze )
  {
    return;
  }
//
//  Pick any node from the last level that has minimum degree
//  as the starting point to generate a new level set.
//
  for ( ; ; )
  {
    mindeg = iccsze;

    jstrt = level_row[*level_num-1];
    *root = level[jstrt-1];

    if ( jstrt < iccsze )
    {
      for ( j = jstrt; j <= iccsze; j++ )
      {
        node = level[j-1];
        ndeg = 0;
        kstrt = adj_row[node-1];
        kstop = adj_row[node] - 1;

        for ( k = kstrt; k <= kstop; k++ )
        {
          nabor = adj[k-1];
          if ( 0 < mask[nabor-1] )
          {
            ndeg = ndeg + 1;
          }
        }

        if ( ndeg < mindeg )
        {
          *root = node;
          mindeg = ndeg;
        }
      }
    }
//
//  Generate the rooted level structure associated with this node.
//
    level_set ( *root, adj_num, adj_row, adj, mask, &level_num2,
      level_row, level, node_num );
//
//  If the number of levels did not increase, accept the new ROOT.
//
    if ( level_num2 <= *level_num )
    {
      break;
    }

    *level_num = level_num2;
//
//  In the unlikely case that ROOT is one endpoint of a line graph,
//  we can exit now.
//
    if ( iccsze <= *level_num )
    {
      break;
    }
  }

  return;
}
//****************************************************************************80

void SystemRCMReordererProcess::sort_heap_external ( int n, int *indx, int *i, int *j, int isgn )

//****************************************************************************80
//
//  Purpose:
//
//    SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
//
//  Discussion:
//
//    The actual list is not passed to the routine.  Hence it may
//    consist of integers, reals, numbers, names, etc.  The user,
//    after each return from the routine, will be asked to compare or
//    interchange two items.
//
//    The current version of this code mimics the FORTRAN version,
//    so the values of I and J, in particular, are FORTRAN indices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 February 2004
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the length of the input list.
//
//    Input/output, int *INDX.
//    The user must set INDX to 0 before the first call.
//    On return,
//      if INDX is greater than 0, the user must interchange
//      items I and J and recall the routine.
//      If INDX is less than 0, the user is to compare items I
//      and J and return in ISGN a negative value if I is to
//      precede J, and a positive value otherwise.
//      If INDX is 0, the sorting is done.
//
//    Output, int *I, *J.  On return with INDX positive,
//    elements I and J of the user's list should be
//    interchanged.  On return with INDX negative, elements I
//    and J are to be compared by the user.
//
//    Input, int ISGN. On return with INDX negative, the
//    user should compare elements I and J of the list.  If
//    item I is to precede item J, set ISGN negative,
//    otherwise set ISGN positive.
//
{
  static int i_save = 0;
  static int j_save = 0;
  static int k = 0;
  static int k1 = 0;
  static int n1 = 0;
//
//  INDX = 0: This is the first call.
//
  if ( *indx == 0 )
  {

    i_save = 0;
    j_save = 0;
    k = n / 2;
    k1 = k;
    n1 = n;
  }
//
//  INDX < 0: The user is returning the results of a comparison.
//
  else if ( *indx < 0 )
  {
    if ( *indx == -2 ) 
    {
      if ( isgn < 0 ) 
      {
        i_save = i_save + 1;
      }
      j_save = k1;
      k1 = i_save;
      *indx = -1;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( 0 < isgn ) 
    {
      *indx = 2;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( k <= 1 ) 
    {
      if ( n1 == 1 ) 
      {
        i_save = 0;
        j_save = 0;
        *indx = 0;
      }
      else 
      {
        i_save = n1;
        j_save = 1;
        n1 = n1 - 1;
        *indx = 1;
      }
      *i = i_save;
      *j = j_save;
      return;
    }
    k = k - 1;
    k1 = k;
  }
//
//  0 < INDX: the user was asked to make an interchange.
//
  else if ( *indx == 1 ) 
  {
    k1 = k;
  }

  for ( ; ; )
  {

    i_save = 2 * k1;

    if ( i_save == n1 ) 
    {
      j_save = k1;
      k1 = i_save;
      *indx = -1;
      *i = i_save;
      *j = j_save;
      return;
    }
    else if ( i_save <= n1 ) 
    {
      j_save = i_save + 1;
      *indx = -2;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( k <= 1 ) 
    {
      break;
    }

    k = k - 1;
    k1 = k;
  }

  if ( n1 == 1 ) 
  {
    i_save = 0;
    j_save = 0;
    *indx = 0;
    *i = i_save;
    *j = j_save;
  }
  else 
  {
    i_save = n1;
    j_save = 1;
    n1 = n1 - 1;
    *indx = 1;
    *i = i_save;
    *j = j_save;
  }

  return;
}
//****************************************************************************80

void SystemRCMReordererProcess::timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}

}

#endif

