#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include "loadmm.h"
#include "SparseMatrix.h"

#define GRN   "\x1B[32m"
#define RED   "\x1B[31m"
#define RESET "\x1B[0m"
#define BLU   "\x1B[34m"

int main(int argc, char** argv) {

    //char *filename = "../data/4x4.mtx";
    char filename[60];
    clock_t t;

    if(argv[1] == NULL){
        fprintf(stderr, RED "Usage: ./1/SpMV.exe YourMatrixMarketFile.mtx\n" RESET);
        return 1;
    }else{
        strcpy(filename, argv[1]);
    }

    int failed;
    //mm load
	mm_file_t *mm_file = loadmm(filename);
    
    //show basic information of loaded matrix
    printf(BLU "\nLoaded Matrix: \n");
    printf(BLU "    name: %s\n", filename);
    printf(BLU "    dim: %d x %d\n", mm_file->nrow, mm_file->ncol);
    printf(BLU "    nnz: %d\n\n", mm_file->data_size);


    double *u = malloc(mm_file->nrow * sizeof(double));

    for(int i = 0; i < mm_file->nrow; i++){
    	u[i] = 1.0;
    }

    double *v = malloc(mm_file->nrow * sizeof(double)); 

    //dense
    dense_t *dense;
    failed = dense_from_mm_buffer(&dense, mm_file); 
    
    if(failed)
    {
        printf(RED"dense    : failed for allocation\n"RESET);
    }else{
        
        for(int i = 0; i < mm_file->nrow; i++){
            v[i] = 0.0;
        }
        t = clock();
        dense_spmv(dense, u, v);
        t = clock() - t;
        printf(GRN "MV (DENSE) time: %fs\n", ((double) t )/CLOCKS_PER_SEC);
        dense_free(dense);
    }

    //coo
    coo_t *coo;
    failed = coo_from_mm_buffer(&coo, mm_file); 
    
    if(failed)
    {
        printf(RED"coo    : failed for allocation\n"RESET);
    }else{
        
        for(int i = 0; i < mm_file->nrow; i++){
            v[i] = 0.0;
        }
        t = clock();
        coo_spmv(coo, u, v);
        t = clock() - t;
        printf(GRN "SpMV (COO) time: %fs\n", ((double) t )/CLOCKS_PER_SEC);
        coo_free(coo);
    }

    //csr
    csr_t *csr;
    failed = csr_from_mm_buffer(&csr, mm_file); 
    
    if(failed)
    {
        printf(RED"csr    : failed for allocation\n"RESET);
    }else{
        
        for(int i = 0; i < mm_file->nrow; i++){
            v[i] = 0.0;
        }
        t = clock();
        csr_spmv(csr, u, v);
        t = clock() - t;
        printf(GRN "SpMV (CSR) time: %fs\n", ((double) t )/CLOCKS_PER_SEC);
        csr_free(csr);
    }


    //ell
    ell_t *ell;
    failed = ell_from_mm_buffer(&ell, mm_file); 
    
    if(failed)
    {
        printf(RED"ell    : failed for allocation\n"RESET);
    }else{
        
        for(int i = 0; i < mm_file->nrow; i++){
            v[i] = 0.0;
        }
        t = clock();
        ell_spmv(ell, u, v);
        t = clock() - t;
        printf(GRN "SpMV (ELL) time: %fs\n", ((double) t )/CLOCKS_PER_SEC);
        ell_free(ell);
    }


    //dia
    dia_t *dia;
    failed = dia_from_mm_buffer(&dia, mm_file); 
    
    if(failed)
    {
        printf(RED"dia    : failed for allocation\n"RESET);
    }else{
        
        for(int i = 0; i < mm_file->nrow; i++){
            v[i] = 0.0;
        }
        t = clock();
        dia_spmv(dia, u, v);
        t = clock() - t;
        printf(GRN "SpMV (DIA) time: %fs\n", ((double) t )/CLOCKS_PER_SEC);
        dia_free(dia);
    }


    free(u); free(v); freemm(mm_file);

    return 0;
}