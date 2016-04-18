/*          
 * =======================================================================*
 * kkkk   kkkk  kkkkkkkkkk   kkkkk    kkkkkkkkkk kkkkkkkkkk kkkkkkkkkK    *
 * kkkk  kkkk   kkkk   kkkk  kkkkkk   kkkkkkkkkk kkkkkkkkkk kkkkkkkkkK    *
 * kkkkkkkkk    kkkk   kkkk  kkkkkkk     kkkk    kkk    kkk  kkkk         *
 * kkkkkkkkk    kkkkkkkkkkk  kkkk kkk    kkkk    kkk    kkk    kkkk       *
 * kkkk  kkkk   kkkk  kkkk   kkkk kkkk   kkkk    kkk    kkk      kkkk     *
 * kkkk   kkkk  kkkk   kkkk  kkkk  kkkk  kkkk    kkkkkkkkkk  kkkkkkkkkk   *
 * kkkk    kkkk kkkk    kkkk kkkk   kkkk kkkk    kkkkkkkkkk  kkkkkkkkkk   *
 *                                                                        *
 * krATos: a fREe opEN sOURce CoDE for mULti-pHysIC aDaptIVe SoLVErS,     *
 * aN extEnsIBLe OBjeCt oRiEnTEd SOlutION fOR fInITe ELemEnt fORmULatIONs *
 * Copyleft by 2003 ciMNe                                                 *
 * Copyleft by 2003 originary authors Copyleft by 2003 your name          *
 * This library is free software; you can redistribute it and/or modify   *
 * it under the terms of the GNU Lesser General Public License as         *
 * published by the Free Software Foundation; either version 2.1 of       *
 * the License, or any later version.                                     *
 *                                                                        *
 * This library is distributed in the hope that it will be useful, but    *
 * WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                   *
 * See the GNU Lesser General Public License for more details.            *
 *                                                                        *
 * You should have received a copy of the GNU Lesser General Public       *
 * License along with this library; if not, write to International Centre *
 * for Numerical Methods in Engineering (CIMNE),                          *
 * Edifici C1 - Campus Nord UPC, Gran Capità s/n, 08034 Barcelona.        *
 *                                                                        *
 * You can also contact us to the following email address:                *
 * kratos@cimne.upc.es                                                    *
 * or fax number: +34 93 401 65 17                                        *
 *                                                                        *
 * Created at Institute for Structural Mechanics                          *
 * Ruhr-University Bochum, Germany                                        *
 * Last modified by:    $Author: hbui $                                   *
 * Date:                $Date: Sep 3, 2014 $                              *
 * Revision:            $Revision: 1.0 $                                  *
 *========================================================================*
 * International Center of Numerical Methods in Engineering - CIMNE       *
 * Barcelona - Spain                                                      *
 *========================================================================*
 */

#if !defined(KRATOS_MULTITHREADED_SOLVERS_APPLICATION_FEAST_SOLVER_H_INCLUDED )
#define  KRATOS_MULTITHREADED_SOLVERS_APPLICATION_FEAST_SOLVER_H_INCLUDED

// System includes
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// External includes
#include <boost/timer.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/numeric/bindings/traits/sparse_traits.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_sparse.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include "feast.h"
#include "feast_sparse.h"

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "utilities/openmp_utils.h"

namespace ublas = boost::numeric::ublas;

namespace Kratos {

template<class TSparseSpaceType, class TDenseSpaceType>
class FeastSolver
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(FeastSolver);

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    /**
     * @param niter number of iterative refinements allowed
     */
    FeastSolver()
    {
    }

    /**
     * Destructor
     */
    virtual ~FeastSolver() {}

    /**
     * Find the nlambda largest eigenvalues of A
     * @param rA        System matrix
     * @param nlambda   initial guess of number of eigenvalues in [a, b]
     * @param rLambda   vector contains eigenvalues
     * @param emin      lower bound to find eigenvalues
     * @param emax      upper bound to find eigenvalues
     */
    bool Solve(SparseMatrixType& rA, int& nlambda, VectorType& rLambda, double emin, double emax)
    {
//        double start_solver = OpenMPUtils::GetCurrentTime();
//        
//        int n = TSparseSpaceType::Size1(rA);
//        
//        /* nonzeros in rA */
//        double* a = rA.value_data().begin();

//        /* manual index vector generation */
//        int* ia = new int[rA.index1_data().size()];
//        int* ja = new int[rA.index2_data().size()];
//        std::cout << "Feast: size of the problem: " << n << std::endl;
//        std::cout << "Feast: size of ia: " << rA.index1_data().size() << std::endl;
//        std::cout << "Feast: size of ja: " << rA.index2_data().size() << std::endl;
//        // ia is rowptr
//        for (unsigned int i = 0; i < rA.index1_data().size(); ++i)
//            ia[i] = (int) (rA.index1_data()[i]) + 1;
//        // ja is colind
//        for (unsigned int i = 0; i < rA.index2_data().size(); ++i)
//            ja[i] = (int) (rA.index2_data()[i]) + 1;
//        
//        /*!!!!!!!!!!!!!!!!! Feast variable */
//        int  feastparam[64];
//        double epsout;
//        int loop;
//        char UPLO = 'F';
//        int  i, k, err;
//        int  M0, M, info;
//        double Emin, Emax, trace;
//        double *E;   // eigenvalues  
//        double *res; // residual
//        double *X;   // eigenvectors
//        
//        Emin = emin;
//        Emax = emax;
//        M0   = nlambda;

//        /*!!!!!!!!!!!!! ALLOCATE VARIABLE */
//        E   = new double[M0];
//        res = new double[M0];
//        X   = new double[n*M0];
//        
//        /*!!!!!!!!!!!!  FEAST */
//        feastinit(feastparam);
//        feastparam[0] = 1;  /* Print runtime comments on screen (0: No, 1: Yes) */
//        dfeast_scsrev(&UPLO, &n, a, ia, ja, feastparam, &epsout, &loop, &Emin, &Emax, &M0, E, X, &M, res, &info);
//        
//        /*!!!!!!!!!! REPORT !!!!!!!!!*/
//        printf("FEAST OUTPUT INFO %d\n", info);
//        if (info == 0)
//        {
//            printf("*************************************************\n");
//            printf("************** REPORT ***************************\n");
//            printf("*************************************************\n");
//            printf("# Search interval [%f, %f]\n", Emin, Emax);
//            printf("# mode found/subspace %d %d \n", M, M0);
//            printf("# iterations %d \n", loop);
//            trace = 0.0;
//            for (i = 0; i <= M - 1; ++i)
//                trace = trace + *(E+i);
//            printf("TRACE %f\n", trace);
//            printf("Relative error on the Trace %f\n",epsout );
//            printf("Eigenvalues/Residuals\n");
//            for (i = 0; i <= M-1; ++i)
//                printf("   %d %f %f\n", i, *(E+i), *(res+i));
//        }
//        else
//            KRATOS_THROW_ERROR(std::logic_error, "Error at FEAST, error code =", info);
//        
//        delete [] E;
//        delete [] X;
//        delete [] res;
//        delete [] ia;
//        delete [] ja;

//        // post-process the results
//        nlambda = M; //number of eigenvalues found
//        TSparseSpaceType::Resize(rLambda, nlambda);
//        for(int i = 0; i < nlambda; ++i)
//            rLambda(i) = E[i];

//        std::cout << "#### EIGENSOLVER TIME: " << OpenMPUtils::GetCurrentTime() - start_solver << " ####" << std::endl;
//        return true;



        /*!!!!!!!!!!!!!!!!! FeastSolver declaration variable */
          int  feastparam[64]; 
          double epsout;
          int loop;
          char UPLO='F'; 

          /*!!!!!!!!!!!!!!!!! Matrix declaration variable */
          FILE *fp;
          char name[]="system2";
          int  N,nnz;
          double *sa;
          int *isa,*jsa;
          /*!!!!!!!!!!!!!!!!! Others */
          struct timeval t1, t2;
          int  i,k,err;
          int  M0,M,info;
          double Emin,Emax,trace,dummy;
          double *X; //! eigenvectors
          double *E,*res; //! eigenvalue+residual


          /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!! read input file in csr format!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

          // !!!!!!!!!! form CSR arrays isa,jsa,sa 
          fp = fopen (name, "r");
          err=fscanf (fp, "%d%d\n",&N,&nnz);
          sa=(double*)calloc(nnz,sizeof(double));
          isa=(int*)calloc(N+1,sizeof(int));
          jsa=(int*)calloc(nnz,sizeof(int));

          for (i=0;i<=N;i++){
            *(isa+i)=0;
          };
          *(isa)=1;
          for (k=0;k<=nnz-1;k++){
            err=fscanf(fp,"%d%d%lf%lf\n",&i,jsa+k,sa+k,&dummy);
            *(isa+i)=*(isa+i)+1;
          };
          fclose(fp);
          for (i=1;i<=N;i++){
            *(isa+i)=*(isa+i)+*(isa+i-1);
          };
          
          /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!! INFORMATION ABOUT MATRIX !!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
          printf("sparse matrix -system1- size %.d\n",N);
          printf("nnz %d \n",nnz);

          /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!! FEAST in sparse format !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

          /*!!! search interval [Emin,Emax] including M eigenpairs*/
          Emin=(double) -0.1;
          Emax=(double) 0.1;
          M0=100; // !! M0>=M

          /*!!!!!!!!!!!!! ALLOCATE VARIABLE */
          E=(double*)calloc(M0,sizeof(double));  // eigenvalues
          res=(double*)calloc(M0,sizeof(double));// eigenvectors 
          X=(double*)calloc(N*M0,sizeof(double));// residual


          /*!!!!!!!!!!!!  FEAST */
          feastinit(feastparam);
          feastparam[0]=1;  /*change from default value */
          dfeast_scsrev(&UPLO,&N,sa,isa,jsa,feastparam,&epsout,&loop,&Emin,&Emax,&M0,E,X,&M,res,&info);

          /*!!!!!!!!!! REPORT !!!!!!!!!*/
          printf("FEAST OUTPUT INFO %d\n",info);
          if (info==0) {
            printf("*************************************************\n");
            printf("************** REPORT ***************************\n");
            printf("*************************************************\n");
            printf("SIMULATION TIME %f\n",(t2.tv_sec-t1.tv_sec)*1.0+(t2.tv_usec-t1.tv_usec)*0.000001);
            printf("# Search interval [Emin,Emax] %.15e %.15e\n",Emin,Emax);
            printf("# mode found/subspace %d %d \n",M,M0);
            printf("# iterations %d \n",loop);
            trace=(double) 0.0;
            for (i=0;i<=M-1;i=i+1){
              trace=trace+*(E+i);
            }	  
            printf("TRACE %.15e\n", trace);
            printf("Relative error on the Trace %.15e\n",epsout );
            printf("Eigenvalues/Residuals\n");
            for (i=0;i<=M-1;i=i+1){
              printf("   %d %.15e %.15e\n",i,*(E+i),*(res+i));
            }
          }






    }

    /// Return information about this object.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "FEAST solver";
        return buffer.str();
    }
    
    /**
     * Print information about this object.
     */
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "FEAST solver finished.";
    }

    /**
     * Print object's data.
     */
    void PrintData(std::ostream& rOStream) const
    {
    }

private:

    /**
     * Assignment operator.
     */
    FeastSolver& operator=(const FeastSolver& Other);

    /**
     * Copy constructor.
     */

};
// Class FeastSolver

} // namespace Kratos.

#endif // KRATOS_MULTITHREADED_SOLVERS_APPLICATION_FEAST_SOLVER_H_INCLUDED  defined
