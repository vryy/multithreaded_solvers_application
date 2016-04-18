#if !defined(MESH_RCM_H)
#define MESH_RCM_H

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>

#include "includes/model_part.h"

namespace Kratos
{

class MeshRCM {

public:
    MeshRCM() {}
    ~MeshRCM() {}
    
    KRATOS_CLASS_POINTER_DEFINITION(MeshRCM);
    
    /**
     * Renumber a tetrahedra mesh. The model_part is assumed to contains only tetrahedra
     */
    void Renumber(ModelPart& r_model_part);
    
private:
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
    void add_pair(int type, std::vector<int> nodes, std::vector<int>& pair, int& pair_num);
    void mesh_adj_count ( int node_num, int element_num, std::vector<std::vector<int> > element_nodes, std::vector<int> element_type, int *adj_num, int adj_row[] );
    int *mesh_adj_set ( int node_num, int element_num, std::vector<std::vector<int> > element_nodes, std::vector<int> element_type, int adj_num, int adj_row[] );
    void timestamp ( );

};

}

#endif // MESH_RCM_H

