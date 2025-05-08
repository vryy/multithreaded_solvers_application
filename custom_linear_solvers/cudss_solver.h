/*
 * see multithreaded_solvers_application/LICENSE.txt                      *
 *========================================================================*
 * Last modified by:    $Author: hbui $                                   *
 * Date:                $Date: Feb 9, 2025 $                              *
 * Revision:            $Revision: 1.0 $                                  *
 *========================================================================*
 */

#if !defined(KRATOS_MULTITHREADED_SOLVERS_APPLICATION_CUDSS_SOLVER_H_INCLUDED )
#define  KRATOS_MULTITHREADED_SOLVERS_APPLICATION_CUDSS_SOLVER_H_INCLUDED

// System includes

// External includes
#include <cuda_runtime.h>
#include "cudss.h"

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "utilities/openmp_utils.h"

#define ENABLE_PROFILING

#define CUDSS_SOLVER_FREE \
    do { \
        cudaFree(csr_offsets_d); \
        cudaFree(csr_columns_d); \
        cudaFree(csr_values_d); \
        cudaFree(x_values_d); \
        cudaFree(b_values_d); \
    } while(0);

#define CUDA_CALL_AND_CHECK(call, msg) \
    do { \
        cuda_error = call; \
        if (cuda_error != cudaSuccess) { \
            printf("Example FAILED: CUDA API returned error = %d, details: " #msg "\n", cuda_error); \
            CUDSS_SOLVER_FREE; \
            return -1; \
        } \
    } while(0);


#define CUDSS_CALL_AND_CHECK(call, status, msg) \
    do { \
        status = call; \
        if (status != CUDSS_STATUS_SUCCESS) { \
            printf("Example FAILED: CUDSS call ended unsuccessfully with status = %d, details: " #msg "\n", status); \
            CUDSS_SOLVER_FREE; \
            return -2; \
        } \
    } while(0);

namespace Kratos
{

template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class cuDSSSolver: public DirectSolver<TSparseSpaceType, TDenseSpaceType, ModelPart, TReordererType>
{
public:
    /**
     * Counted pointer of SuperLUSolver
     */
    KRATOS_CLASS_POINTER_DEFINITION(cuDSSSolver);

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename BaseType::SparseMatrixType SparseMatrixType;

    typedef typename BaseType::VectorType VectorType;

    typedef typename BaseType::DenseMatrixType DenseMatrixType;

    /**
     * @param niter number of iterative refinements allowed
     */
    cuDSSSolver()
    {
        std::cout << "cuDSSSolver created" << std::endl;
        this->PrintDeviceInfo();
    }

    /**
     * Destructor
     */
    ~cuDSSSolver() override {}

    /**
     * [derived]
     */
    bool AdditionalPhysicalDataIsNeeded() override
    {
        return false;
    }

    /**
     * Normal solve method.
     * Solves the linear system Ax=b and puts the result on SystemVector& rX.
     * rX is also th initial guess for iterative methods.
     * @param rA. System matrix
     * @param rX. Solution vector.
     * @param rB. Right hand side vector.
     */
    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        double start_solver = OpenMPUtils::GetCurrentTime();

        const int n = rA.size1();
        assert(n == rA.size2());
        assert(n == rB.size());
        assert(n == rX.size());

        /* nonzeros in rA */
        double* a = rA.value_data().begin();

        /* manual index vector generation */
        std::vector<int> index1_vector(rA.index1_data().size());
        std::vector<int> index2_vector(rA.index2_data().size());
        std::cout << "Size of the problem: " << n << std::endl;
        std::cout << "Size of index1_vector: " << rA.index1_data().size() << std::endl;
        std::cout << "Size of index2_vector: " << rA.index2_data().size() << std::endl;
        // index1_vector is rowptr
        for (unsigned int i = 0; i < rA.index1_data().size(); ++i)
            index1_vector[i] = (int) (rA.index1_data()[i]);
        // index2_vector is colind
        for (unsigned int i = 0; i < rA.index2_data().size(); ++i)
            index2_vector[i] = (int) (rA.index2_data()[i]);

        const int nnz = (int) index2_vector.size();
        const int nrhs = 1;

        /* transfer the data and solve on the device */

        cudaError_t cuda_error = cudaSuccess;
        cudssStatus_t status = CUDSS_STATUS_SUCCESS;

        int *csr_offsets_h = index1_vector.data();
        int *csr_columns_h = index2_vector.data();
        double *csr_values_h = a;
        double *x_values_h = &rX[0], *b_values_h = &rB[0];

        // create and allocate the data on the device
        int *csr_offsets_d = nullptr;
        int *csr_columns_d = nullptr;
        double *csr_values_d = nullptr;
        double *x_values_d = nullptr, *b_values_d = nullptr;

        #ifdef ENABLE_PROFILING
        double start_device_alloc = OpenMPUtils::GetCurrentTime();
        #endif

        CUDA_CALL_AND_CHECK(cudaMalloc(&csr_offsets_d, (n + 1) * sizeof(int)), "cudaMalloc for csr_offsets");
        CUDA_CALL_AND_CHECK(cudaMalloc(&csr_columns_d, nnz * sizeof(int)), "cudaMalloc for csr_columns");
        CUDA_CALL_AND_CHECK(cudaMalloc(&csr_values_d, nnz * sizeof(double)), "cudaMalloc for csr_values");
        CUDA_CALL_AND_CHECK(cudaMalloc(&b_values_d, nrhs * n * sizeof(double)), "cudaMalloc for b_values");
        CUDA_CALL_AND_CHECK(cudaMalloc(&x_values_d, nrhs * n * sizeof(double)), "cudaMalloc for x_values");

        #ifdef ENABLE_PROFILING
        double end_device_alloc = OpenMPUtils::GetCurrentTime();
        std::cout << "-- Device memory allocation: " << end_device_alloc - start_device_alloc << "s" << std::endl;
        #endif

        #ifdef ENABLE_PROFILING
        double start_device_copy = OpenMPUtils::GetCurrentTime();
        #endif

        /* Copy host memory to device for A and b */
        CUDA_CALL_AND_CHECK(cudaMemcpy(csr_offsets_d, csr_offsets_h, (n + 1) * sizeof(int), cudaMemcpyHostToDevice), "cudaMemcpy for csr_offsets");
        CUDA_CALL_AND_CHECK(cudaMemcpy(csr_columns_d, csr_columns_h, nnz * sizeof(int), cudaMemcpyHostToDevice), "cudaMemcpy for csr_columns");
        CUDA_CALL_AND_CHECK(cudaMemcpy(csr_values_d, csr_values_h, nnz * sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy for csr_values");
        CUDA_CALL_AND_CHECK(cudaMemcpy(b_values_d, b_values_h, nrhs * n * sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy for b_values");

        #ifdef ENABLE_PROFILING
        double end_device_copy = OpenMPUtils::GetCurrentTime();
        std::cout << "-- Copying data to device: " << end_device_copy - start_device_copy << "s" << std::endl;
        #endif

        /* Create a CUDA stream */
        cudaStream_t stream = NULL;
        CUDA_CALL_AND_CHECK(cudaStreamCreate(&stream), "cudaStreamCreate");

        /* Creating the cuDSS library handle */
        cudssHandle_t handle;

        CUDSS_CALL_AND_CHECK(cudssCreate(&handle), status, "cudssCreate");

        /* (optional) Setting the custom stream for the library handle */
        CUDSS_CALL_AND_CHECK(cudssSetStream(handle, stream), status, "cudssSetStream");

        /* Creating cuDSS solver configuration and data objects */
        cudssConfig_t solverConfig;
        cudssData_t solverData;

        CUDSS_CALL_AND_CHECK(cudssConfigCreate(&solverConfig), status, "cudssConfigCreate");
        CUDSS_CALL_AND_CHECK(cudssDataCreate(handle, &solverData), status, "cudssDataCreate");

        /* Create matrix objects for the right-hand side b and solution x (as dense matrices). */
        cudssMatrix_t x, b;

        int64_t nrows = n, ncols = n;
        int ldb = ncols, ldx = nrows;
        CUDSS_CALL_AND_CHECK(cudssMatrixCreateDn(&b, ncols, nrhs, ldb, b_values_d, CUDA_R_64F,
                             CUDSS_LAYOUT_COL_MAJOR), status, "cudssMatrixCreateDn for b");
        CUDSS_CALL_AND_CHECK(cudssMatrixCreateDn(&x, nrows, nrhs, ldx, x_values_d, CUDA_R_64F,
                             CUDSS_LAYOUT_COL_MAJOR), status, "cudssMatrixCreateDn for x");

        /* Create a matrix object for the sparse input matrix. */
        cudssMatrix_t A;
        cudssMatrixType_t mtype     = CUDSS_MTYPE_SPD;
        cudssMatrixViewType_t mview = CUDSS_MVIEW_UPPER;
        cudssIndexBase_t base       = CUDSS_BASE_ZERO;
        CUDSS_CALL_AND_CHECK(cudssMatrixCreateCsr(&A, nrows, ncols, nnz, csr_offsets_d, NULL,
                             csr_columns_d, csr_values_d, CUDA_R_32I, CUDA_R_64F, mtype, mview,
                             base), status, "cudssMatrixCreateCsr");

        #ifdef ENABLE_PROFILING
        double start_device_symb = OpenMPUtils::GetCurrentTime();
        #endif

        /* Symbolic factorization */
        CUDSS_CALL_AND_CHECK(cudssExecute(handle, CUDSS_PHASE_ANALYSIS, solverConfig, solverData,
                             A, x, b), status, "cudssExecute for analysis");

        #ifdef ENABLE_PROFILING
        double end_device_symb = OpenMPUtils::GetCurrentTime();
        std::cout << "-- Symbolic factorization: " << end_device_symb - start_device_symb << "s" << std::endl;
        double start_device_fact = OpenMPUtils::GetCurrentTime();
        #endif

        /* Factorization */
        CUDSS_CALL_AND_CHECK(cudssExecute(handle, CUDSS_PHASE_FACTORIZATION, solverConfig,
                             solverData, A, x, b), status, "cudssExecute for factor");

        #ifdef ENABLE_PROFILING
        double end_device_fact = OpenMPUtils::GetCurrentTime();
        std::cout << "-- Numerical factorization: " << end_device_fact - start_device_fact << "s" << std::endl;
        double start_device_solve = OpenMPUtils::GetCurrentTime();
        #endif

        /* Solving */
        CUDSS_CALL_AND_CHECK(cudssExecute(handle, CUDSS_PHASE_SOLVE, solverConfig, solverData,
                             A, x, b), status, "cudssExecute for solve");

        #ifdef ENABLE_PROFILING
        double end_device_solve = OpenMPUtils::GetCurrentTime();
        std::cout << "-- Back solve: " << end_device_solve - start_device_solve << "s" << std::endl;
        #endif

        /* Destroying opaque objects, matrix wrappers and the cuDSS library handle */
        CUDSS_CALL_AND_CHECK(cudssMatrixDestroy(A), status, "cudssMatrixDestroy for A");
        CUDSS_CALL_AND_CHECK(cudssMatrixDestroy(b), status, "cudssMatrixDestroy for b");
        CUDSS_CALL_AND_CHECK(cudssMatrixDestroy(x), status, "cudssMatrixDestroy for x");
        CUDSS_CALL_AND_CHECK(cudssDataDestroy(handle, solverData), status, "cudssDataDestroy");
        CUDSS_CALL_AND_CHECK(cudssConfigDestroy(solverConfig), status, "cudssConfigDestroy");
        CUDSS_CALL_AND_CHECK(cudssDestroy(handle), status, "cudssHandleDestroy");

        CUDA_CALL_AND_CHECK(cudaStreamSynchronize(stream), "cudaStreamSynchronize");

        #ifdef ENABLE_PROFILING
        double start_host_copy = OpenMPUtils::GetCurrentTime();
        #endif

        /* Copying back the solution to host */
        CUDA_CALL_AND_CHECK(cudaMemcpy(x_values_h, x_values_d, nrhs * n * sizeof(double),
                            cudaMemcpyDeviceToHost), "cudaMemcpy for x_values");

        #ifdef ENABLE_PROFILING
        double end_host_copy = OpenMPUtils::GetCurrentTime();
        std::cout << "-- Copying solution to host: " << end_host_copy - start_host_copy << "s" << std::endl;
        #endif

        /* Release the data allocated on the user side */
        CUDSS_SOLVER_FREE;

        std::cout << "#### SOLVER TIME: " << OpenMPUtils::GetCurrentTime() - start_solver << " ####" << std::endl;
        return (status == CUDSS_STATUS_SUCCESS);
    }

    /**
     * Multi solve method for solving a set of linear systems with same coefficient matrix.
     * Solves the linear system Ax=b and puts the result on SystemVector& rX.
     * rX is also th initial guess for iterative methods.
     * @param rA. System matrix
     * @param rX. Solution vector.
     * @param rB. Right hand side vector.
     */
    bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB) override
    {
        KRATOS_ERROR << "This solver can be used for single RHS only";
        return false;
    }

    /// Return information about this object.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "cuDSS solver";
        return buffer.str();
    }

    /**
     * Print information about this object.
     */
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "cuDSS solver finished.";
    }

    /**
     * Print object's data.
     */
    void PrintData(std::ostream& rOStream) const override
    {
    }

private:

    void PrintDeviceInfo() const
    {
        int deviceCount;
        cudaGetDeviceCount(&deviceCount);
        if (deviceCount == 0)
        {
            std::cout << "No CUDA devices found.\n";
            return;
        }

        for (int i = 0; i < deviceCount; ++i)
        {
            cudaDeviceProp prop;
            cudaGetDeviceProperties(&prop, i);
            std::cout << "Device " << i << ": " << prop.name << "\n";
            std::cout << "  Compute Capability: " << prop.major << "." << prop.minor << "\n";
            std::cout << "  Total Global Memory: " << prop.totalGlobalMem / (1024 * 1024) << " MB\n";
            std::cout << "  Shared Memory per Block: " << prop.sharedMemPerBlock / 1024 << " KB\n";
            std::cout << "  Max Threads per Block: " << prop.maxThreadsPerBlock << "\n";
            std::cout << "  Multi-Processor Count: " << prop.multiProcessorCount << "\n";
            std::cout << "  Max Grid Size: " << prop.maxGridSize[0] << " x " << prop.maxGridSize[1] << " x " << prop.maxGridSize[2] << "\n";
            std::cout << "  Max Threads per Multiprocessor: " << prop.maxThreadsPerMultiProcessor << "\n\n";
        }
    }

    /**
     * Assignment operator.
     */
    cuDSSSolver& operator=(const cuDSSSolver& Other);
}; // Class cuDSSSolver

} // namespace Kratos.

#undef ENABLE_PROFILING
#undef CUDA_CALL_AND_CHECK
#undef CUDSS_CALL_AND_CHECK

#endif // KRATOS_MULTITHREADED_SOLVERS_APPLICATION_CUDSS_SOLVER_H_INCLUDED  defined
