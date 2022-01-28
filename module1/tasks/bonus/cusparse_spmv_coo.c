#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <cuda_runtime_api.h>
#include <cusparse.h>         
#include "loadmm.h"
#include "SparseMatrix.h"

#define CHECK_CUDA(func)                                                       \
{                                                                              \
    cudaError_t status = (func);                                               \
    if (status != cudaSuccess) {                                               \
        printf("CUDA API failed at line %d with error: %s (%d)\n",             \
               __LINE__, cudaGetErrorString(status), status);                  \
        return EXIT_FAILURE;                                                   \
    }                                                                          \
}

#define CHECK_CUSPARSE(func)                                                   \
{                                                                              \
    cusparseStatus_t status = (func);                                          \
    if (status != CUSPARSE_STATUS_SUCCESS) {                                   \
        printf("CUSPARSE API failed at line %d with error: %s (%d)\n",         \
               __LINE__, cusparseGetErrorString(status), status);              \
        return EXIT_FAILURE;                                                   \
    }                                                                          \
}

int main(int argc, char** argv) {

    char filename[60];
    clock_t t;

    if(argv[1] == NULL){
        fprintf(stderr, "Usage: ./1/SpMV.exe YourMatrixMarketFile.mtx\n" );
        return 1;
    }else{
        strcpy(filename, argv[1]);
    }


    int niter;
    if(argv[2] == NULL){
   	niter = 1 ;
    }else{
	niter = atoi(argv[2]);    
    }


    int failed;
    //mm load
        mm_file_t *mm_file = loadmm(filename);

    //show basic information of loaded matrix
    printf("\nLoaded Matrix: \n");
    printf("    name: %s\n", filename);
    printf("    dim: %d x %d\n", mm_file->nrow, mm_file->ncol);
    printf("    nnz: %d\n\n", mm_file->data_size);

    double alpha = 1.0, beta = 0.0;

    double *u = malloc(mm_file->nrow * sizeof(double));

    for(int i = 0; i < mm_file->nrow; i++){
        u[i] = 1.0;
    }

    double *v = calloc(mm_file->nrow, sizeof(double));

    //coo
    coo_t *coo;
    failed = coo_from_mm_buffer(&coo, mm_file);

    if(failed){
	printf("coo    : failed for allocation\n");
    }else{
        // Device memory management
	int *d_rowind, *d_colind;
	double *d_values, *d_u, *d_v;
        
	t = clock();

	CHECK_CUDA( cudaMalloc((void**) &d_rowind,  (coo->nnz) * sizeof(int)) );
        CHECK_CUDA( cudaMalloc((void**) &d_colind,  coo->nnz * sizeof(int)) );
        CHECK_CUDA( cudaMalloc((void**) &d_values,  coo->nnz * sizeof(double)) );
        CHECK_CUDA( cudaMalloc((void**) &d_u,  coo->m * sizeof(double)) );
        CHECK_CUDA( cudaMalloc((void**) &d_v,  coo->m * sizeof(double)) );

        CHECK_CUDA( cudaMemcpy(d_rowind, coo->rowind, (coo->nnz) * sizeof(int), cudaMemcpyHostToDevice) );
        CHECK_CUDA( cudaMemcpy(d_colind, coo->colind, (coo->nnz) * sizeof(int), cudaMemcpyHostToDevice) );
        CHECK_CUDA( cudaMemcpy(d_values, coo->val, (coo->nnz) * sizeof(double), cudaMemcpyHostToDevice) );
        CHECK_CUDA( cudaMemcpy(d_u, u, (coo->m) * sizeof(double), cudaMemcpyHostToDevice) );
        CHECK_CUDA( cudaMemcpy(d_v, v, (coo->m) * sizeof(double), cudaMemcpyHostToDevice) );

    	cusparseHandle_t     handle = NULL;
    	cusparseSpMatDescr_t coo_descr;
    	cusparseDnVecDescr_t u_descr, v_descr;

    	void*                dBuffer    = NULL;
    	size_t               bufferSize = 0;
    	
	CHECK_CUSPARSE( cusparseCreate(&handle) );

    	CHECK_CUSPARSE( cusparseCreateCoo(&coo_descr, coo->m, coo->m, coo->nnz,
                                      d_rowind, d_colind, d_values,
                                      CUSPARSE_INDEX_32I,
                                      CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F) );

    	CHECK_CUSPARSE( cusparseCreateDnVec(&u_descr, coo->m, d_u, CUDA_R_64F) );
        CHECK_CUSPARSE( cusparseCreateDnVec(&v_descr, coo->m, d_v, CUDA_R_64F) );	


        CHECK_CUSPARSE( cusparseSpMV_bufferSize(
              			 handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                 &alpha, coo_descr, u_descr, &beta, v_descr, CUDA_R_64F,
                                 CUSPARSE_MV_ALG_DEFAULT, &bufferSize) );

        CHECK_CUDA( cudaMalloc(&dBuffer, bufferSize) )
    	
	CHECK_CUSPARSE( cusparseSpMV(handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                 &alpha, coo_descr, u_descr, &beta, v_descr, CUDA_R_64F,
                                 CUSPARSE_MV_ALG_DEFAULT, dBuffer) );

        CHECK_CUSPARSE( cusparseDestroySpMat(coo_descr) );
        CHECK_CUSPARSE( cusparseDestroyDnVec(u_descr) );
        CHECK_CUSPARSE( cusparseDestroyDnVec(v_descr) );
        CHECK_CUSPARSE( cusparseDestroy(handle) );

        CHECK_CUDA( cudaMemcpy(v, d_v, coo->m * sizeof(double), cudaMemcpyDeviceToHost) );

        t = clock() - t;
	printf("cuSparse SpMV (COO) time: %fs\n", ((double) t )/CLOCKS_PER_SEC);
	
        CHECK_CUDA( cudaFree(dBuffer) );
    	CHECK_CUDA( cudaFree(d_rowind) );
    	CHECK_CUDA( cudaFree(d_colind) );
   	CHECK_CUDA( cudaFree(d_values) );
    	CHECK_CUDA( cudaFree(d_u) );
    	CHECK_CUDA( cudaFree(d_v) );

	coo_free(coo);	
    }

    free(u); free(v);
    freemm(mm_file);

    return 0;
}
