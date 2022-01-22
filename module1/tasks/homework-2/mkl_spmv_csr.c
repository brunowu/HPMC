#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include "loadmm.h"
#include "SparseMatrix.h"

#include "mkl_spblas.h"
#include "omp.h"

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

    //csr
    csr_t *csr;
    failed = csr_from_mm_buffer(&csr, mm_file);

    if(failed){
	printf("coo    : failed for allocation\n");
    }else{
    
        //MKL part
	struct matrix_descr descr = {SPARSE_MATRIX_TYPE_GENERAL};
	sparse_matrix_t mkl_csr_handle;
	mkl_sparse_d_create_csr(&mkl_csr_handle, SPARSE_INDEX_BASE_ZERO, csr->m, csr->m, csr->rowoffs, csr->rowoffs + 1, csr->colind, csr->val);

	if(niter > 1){
            mkl_sparse_set_mv_hint(mkl_csr_handle, SPARSE_OPERATION_NON_TRANSPOSE, descr, niter);
            mkl_sparse_set_memory_hint ( mkl_csr_handle, SPARSE_MEMORY_AGGRESSIVE );
            mkl_sparse_optimize(mkl_csr_handle);
        }

        for(int i = 0; i < niter; i++){
            t = clock();
            mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, mkl_csr_handle, descr, u, beta, v);
            t = clock() - t;
            printf("(%d): MKL SpMV (CSR) time: %fs\n", i + 1, ((double) t )/CLOCKS_PER_SEC);
        }

        mkl_sparse_destroy(mkl_csr_handle);
        csr_free(csr);	
    }

    free(u); free(v);
    freemm(mm_file);

    return 0;
}
