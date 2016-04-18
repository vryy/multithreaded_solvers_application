/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*             ********   ***                                 SparseLib++    */
/*          *******  **  ***       ***      ***                              */
/*           *****      ***     ******** ********                            */
/*            *****    ***     ******** ********              R. Pozo        */
/*       **  *******  ***   **   ***      ***                 K. Remington   */
/*        ********   ********                                 A. Lumsdaine   */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                                                                           */
/*                                                                           */
/*                     SparseLib++ : Sparse Matrix Library                   */
/*                                                                           */
/*               National Institute of Standards and Technology              */
/*                        University of Notre Dame                           */
/*              Authors: R. Pozo, K. Remington, A. Lumsdaine                 */
/*                                                                           */
/*                                 NOTICE                                    */
/*                                                                           */
/* Permission to use, copy, modify, and distribute this software and         */
/* its documentation for any purpose and without fee is hereby granted       */
/* provided that the above notice appear in all copies and supporting        */
/* documentation.                                                            */
/*                                                                           */
/* Neither the Institutions (National Institute of Standards and Technology, */
/* University of Notre Dame) nor the Authors make any representations about  */
/* the suitability of this software for any purpose.  This software is       */
/* provided ``as is'' without expressed or implied warranty.                 */
/*                                                                           */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#ifndef ILUTPRE_H
#define ILUTPRE_H

#include "vecdefs.h"
#include "pre.h"
#include VECTOR_H
#include "comprow_double.h"
#include "compcol_double.h"

class CompRow_ILUtPreconditioner_double : public Preconditioner_double {

 private:
  iluptr lu;
  int dim_[2];

 public:
  CompRow_ILUtPreconditioner_double(const CompRow_Mat_double &A, int lfil, double tol);
  ~CompRow_ILUtPreconditioner_double() {cleanILU(lu);}
  
  virtual VECTOR_double     solve(const VECTOR_double &x) const;
  virtual VECTOR_double     trans_solve(const VECTOR_double &x) const;
};

#endif
