#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include "loadmm.h"
#include "SparseMatrix.h"

int main(int argc, char** argv) {

    char *filename = "../data/4x4.mtx";

    int failed;
    //mm load
	mm_file_t *mm_file = loadmm(filename);
    
    //show basic information of loaded matrix
    printf("\nLoaded Matrix: \n");
    printf("    name: %s\n", filename);
    printf("    dim: %d x %d\n", mm_file->nrow, mm_file->ncol);
    printf("    nnz: %d\n\n", mm_file->data_size);


    double *u = malloc(mm_file->nrow * sizeof(double));

    for(int i = 0; i < mm_file->nrow; i++){
    	u[i] = 1.0;
    }

    double *v_dense = calloc(mm_file->nrow, sizeof(double)); 
    double *v_coo = calloc(mm_file->nrow, sizeof(double)); 
    double *v_csr = calloc(mm_file->nrow, sizeof(double)); 
    double *v_dia = calloc(mm_file->nrow, sizeof(double)); 
    double *v_ell = calloc(mm_file->nrow, sizeof(double)); 

    //dense
    dense_t *dense;
    failed = dense_from_mm_buffer(&dense, mm_file); 
    assert(failed == 0);
    
    dense_spmv(dense, u, v_dense);

    //coo
    coo_t *coo;
    failed = coo_from_mm_buffer(&coo, mm_file); 
    assert(failed == 0);
    
    coo_spmv(coo, u, v_coo);
    
    for(int i = 0; i < mm_file->nrow; i++){
        assert(v_dense[i] == v_coo[i]);
    }

    printf("COO SpMV: passed\n");

    //csr
    csr_t *csr;
    failed = csr_from_mm_buffer(&csr, mm_file); 
    assert(failed == 0);
    
    csr_spmv(csr, u, v_csr);

    for(int i = 0; i < mm_file->nrow; i++){
        assert(v_dense[i] == v_csr[i]);
    }

    printf("CSR SpMV: passed\n");
    
    //dia
    dia_t *dia;
    failed = dia_from_mm_buffer(&dia, mm_file); 
    assert(failed == 0);
    
    dia_spmv(dia, u, v_dia);
    
    for(int i = 0; i < mm_file->nrow; i++){
        assert(v_dense[i] == v_coo[i]);
    }

    printf("DIA SpMV: passed\n");
    

    //ell
    ell_t *ell;
    failed = ell_from_mm_buffer(&ell, mm_file); 
    assert(failed == 0);
    
    ell_spmv(ell, u, v_ell);

    for(int i = 0; i < mm_file->nrow; i++){
        assert(v_dense[i] == v_coo[i]);
    }

    printf("ELL SpMV: passed\n");
    
    dense_free(dense);
    coo_free(coo);
    csr_free(csr);
    dia_free(dia);
    ell_free(ell);

    free(u); free(v_dense); free(v_coo); free(v_csr); free(v_dia); free(v_ell);
    freemm(mm_file);

    return 0;
}
