#include "spaces/ublas_space.h"
#include "feast_solver.h"
#include "feast.h"
#include "feast_sparse.h"

namespace Kratos
{

    template<class TSparseSpaceType, class TDenseSpaceType>
    bool FeastSolver<TSparseSpaceType, TDenseSpaceType>::Solve(FeastSolver::SparseMatrixType& rA,
            int nlambda, FeastSolver::VectorType& rLambda, double emin, double emax)
    {
        double start_solver = OpenMPUtils::GetCurrentTime();

        int n = TSparseSpaceType::Size1(rA);

        /* nonzeros in rA */
        double* a = rA.value_data().begin();

        /* manual index vector generation */
        int* ia = new int[rA.index1_data().size()];
        int* ja = new int[rA.index2_data().size()];
        std::cout << "Feast: size of the problem: " << n << std::endl;
        std::cout << "Feast: size of ia: " << rA.index1_data().size() << std::endl;
        std::cout << "Feast: size of ja: " << rA.index2_data().size() << std::endl;
        // ia is rowptr
        for (unsigned int i = 0; i < rA.index1_data().size(); ++i)
            ia[i] = (int) (rA.index1_data()[i]) + 1;
        // ja is colind
        for (unsigned int i = 0; i < rA.index2_data().size(); ++i)
            ja[i] = (int) (rA.index2_data()[i]) + 1;

        /*!!!!!!!!!!!!!!!!! Feast variable */
        int  feastparam[64];
        double epsout;
        int loop;
        char UPLO = 'F';
        int  i, k, err;
        int  M0, M, info;
        double Emin, Emax, trace;
        double *E;   // eigenvalues
        double *res; // residual
        double *X;   // eigenvectors

        Emin = emin;
        Emax = emax;
        M0   = nlambda; // size of subspae

        /*!!!!!!!!!!!!! ALLOCATE VARIABLE */
        E   = new double[M0];
        res = new double[M0];
        X   = new double[n*M0];

        /*!!!!!!!!!!!!  FEAST */
        feastinit(feastparam);
        feastparam[0] = 1;  /* Print runtime comments on screen (0: No, 1: Yes) */
        dfeast_scsrev(&UPLO, &n, a, ia, ja, feastparam, &epsout, &loop, &Emin, &Emax, &M0, E, X, &M, res, &info);

        /*!!!!!!!!!! REPORT !!!!!!!!!*/
        printf("FEAST OUTPUT INFO %d\n", info);
        if(info == 0)
        {
            printf("*************************************************\n");
            printf("************** REPORT ***************************\n");
            printf("*************************************************\n");
            printf("# Search interval [%f, %f]\n", Emin, Emax);
            printf("# mode found/subspace %d %d \n", M, M0);
            printf("# iterations %d \n", loop);
            trace = 0.0;
            for (i = 0; i <= M - 1; ++i)
                trace = trace + *(E+i);
            printf("TRACE %f\n", trace);
            printf("Relative error on the Trace %f\n",epsout );
            printf("Eigenvalues/Residuals\n");
            for (i = 0; i <= M-1; ++i)
                printf("   %d %f %f\n", i, *(E+i), *(res+i));
        }
        else if(info == 1)
        {
            std::cout << "Warning: No Eigenvalue found in the search interval" << std::endl;
        }
        else if(info == 2)
        {
            std::cout << "Warning: No Convergence (#iteration loops>fpm(4))" << std::endl;
        }
        else if(info == 3)
        {
            std::cout << "Warning: Size of the subspace M0 is too small (M0<=M)" << std::endl;
        }
        else if(info == 4)
        {
            std::cout << "Warning: Only the subspace has been returned using fpm(14)=1" << std::endl;
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "Error at FEAST, error code =", info);

        // post-process the results
        nlambda = M; //number of eigenvalues found
        TSparseSpaceType::Resize(rLambda, nlambda);
        for(int i = 0; i < nlambda; ++i)
        {
            rLambda(i) = E[i];
        }

        delete [] E;
        delete [] X;
        delete [] res;
        delete [] ia;
        delete [] ja;

        std::cout << "#### EIGENSOLVER TIME: " << OpenMPUtils::GetCurrentTime() - start_solver << " ####" << std::endl;
        return true;
    }


    template<class TSparseSpaceType, class TDenseSpaceType>
    bool FeastSolver<TSparseSpaceType, TDenseSpaceType>::SolveGeneralized(FeastSolver::SparseMatrixType& rA,
            FeastSolver::SparseMatrixType& rB, int nlambda, FeastSolver::VectorType& rLambda, boost::python::list EigenVectors,
            double emin, double emax)
    {
        double start_solver = OpenMPUtils::GetCurrentTime();

        int n = TSparseSpaceType::Size1(rA);

        /* nonzeros in rA */
        double* a = rA.value_data().begin();

        /* manual index vector generation */
        int* ia = new int[rA.index1_data().size()];
        int* ja = new int[rA.index2_data().size()];
        std::cout << "Feast: size of the problem: " << n << std::endl;
        std::cout << "Feast: size of ia: " << rA.index1_data().size() << std::endl;
        std::cout << "Feast: size of ja: " << rA.index2_data().size() << std::endl;
        // ia is rowptr
        for (unsigned int i = 0; i < rA.index1_data().size(); ++i)
            ia[i] = (int) (rA.index1_data()[i]) + 1;
        // ja is colind
        for (unsigned int i = 0; i < rA.index2_data().size(); ++i)
            ja[i] = (int) (rA.index2_data()[i]) + 1;

        /* nonzeros in rB */
        double* b = rB.value_data().begin();

        /* manual index vector generation */
        int* ib = new int[rB.index1_data().size()];
        int* jb = new int[rB.index2_data().size()];
        // ib is rowptr
        for (unsigned int i = 0; i < rB.index1_data().size(); ++i)
            ib[i] = (int) (rB.index1_data()[i]) + 1;
        // jb is colind
        for (unsigned int i = 0; i < rB.index2_data().size(); ++i)
            jb[i] = (int) (rB.index2_data()[i]) + 1;

        /*!!!!!!!!!!!!!!!!! Feast variable */
        int  feastparam[64];
        double epsout;
        int loop;
        char UPLO = 'F';
        int  i, k, err;
        int  M0, M, info;
        double Emin, Emax, trace;
        double *E;   // eigenvalues
        double *res; // residual
        double *X;   // eigenvectors

        Emin = emin;
        Emax = emax;
        M0   = nlambda; // size of subspae

        /*!!!!!!!!!!!!! ALLOCATE VARIABLE */
        E   = new double[M0];
        res = new double[M0];
        X   = new double[n*M0];

        /*!!!!!!!!!!!!  FEAST */
        feastinit(feastparam);
        feastparam[0] = 1;  /* Print runtime comments on screen (0: No, 1: Yes) */
        // dfeast_scsrev(&UPLO, &n, a, ia, ja, feastparam, &epsout, &loop, &Emin, &Emax, &M0, E, X, &M, res, &info);
        dfeast_scsrgv(&UPLO, &n, a, ia, ja, b, ib, jb, feastparam, &epsout, &loop, &Emin, &Emax, &M0, E, X, &M, res, &info);

        /*!!!!!!!!!! REPORT !!!!!!!!!*/
        printf("FEAST OUTPUT INFO %d\n", info);
        if(info == 0)
        {
            printf("*************************************************\n");
            printf("************** REPORT ***************************\n");
            printf("*************************************************\n");
            printf("# Search interval [%f, %f]\n", Emin, Emax);
            printf("# mode found/subspace %d %d \n", M, M0);
            printf("# iterations %d \n", loop);
            trace = 0.0;
            for (i = 0; i <= M - 1; ++i)
                trace = trace + *(E+i);
            printf("TRACE %f\n", trace);
            printf("Relative error on the Trace %f\n",epsout );
            printf("Eigenvalues/Residuals\n");
            for (i = 0; i <= M-1; ++i)
                printf("   %d %f %f\n", i, *(E+i), *(res+i));
        }
        else if(info == 1)
        {
            std::cout << "Warning: No Eigenvalue found in the search interval" << std::endl;
        }
        else if(info == 2)
        {
            std::cout << "Warning: No Convergence (#iteration loops>fpm(4))" << std::endl;
        }
        else if(info == 3)
        {
            std::cout << "Warning: Size of the subspace M0 is too small (M0<=M)" << std::endl;
        }
        else if(info == 4)
        {
            std::cout << "Warning: Only the subspace has been returned using fpm(14)=1" << std::endl;
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "Error at FEAST, error code =", info);

        // post-process the results
        nlambda = M; //number of eigenvalues found
        TSparseSpaceType::Resize(rLambda, nlambda);
        for(int i = 0; i < nlambda; ++i)
        {
            rLambda(i) = E[i];

            Vector V(n);
            for(int j = 0; j < n; ++j)
                V(j) = X[j + i*n];
            EigenVectors.append(V);
        }

        delete [] E;
        delete [] X;
        delete [] res;
        delete [] ia;
        delete [] ja;
        delete [] ib;
        delete [] jb;

        std::cout << "#### EIGENSOLVER TIME: " << OpenMPUtils::GetCurrentTime() - start_solver << " ####" << std::endl;
        return true;
    }

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    template class FeastSolver<SparseSpaceType, LocalSpaceType>;

}
