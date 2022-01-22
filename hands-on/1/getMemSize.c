#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "loadmm.h"
#include "SparseMatrix.h"

#define GRN   "\x1B[32m"
#define RED   "\x1B[31m"
#define RESET "\x1B[0m"
#define BLU   "\x1B[34m"

int main(int argc, char** argv) {

    // char *filename = "../data/4x4.mtx";
	char filename[60];

    if(argv[1] == NULL){
        fprintf(stderr, RED "Usage: ./1/getMemSize.exe YourMatrixMarketFile.mtx\n" RESET);
        return 1;
    }else{
        strcpy(filename, argv[1]);
    }

    int failed;
    // load matrices from local in MatrixMarket format
	mm_file_t *mm_file = loadmm(filename);

    //show basic information of loaded matrix
    printf(BLU "\nLoaded Matrix: \n");
    printf(BLU "    name: %s\n", filename);
    printf(BLU "    dim: %d x %d\n", mm_file->nrow, mm_file->ncol);
    printf(BLU "    nnz: %d\n\n"RESET, mm_file->data_size);

    printf("\nBytes / non-zeros Entry\n");
    printf("-----------------------------------------\n");
    printf("Matrix : %s\n", filename);
    printf("-----------------------------------------\n");

    //dense
    dense_t *dense;
    // convert loaded matrix into dense format
    failed = dense_from_mm_buffer(&dense, mm_file); 
    if(!failed){
        printf("dense  : %.2f bytes\n", (double)dense_getMemSize(dense) / (double)mm_file->nnz);
        printf("-----------------------------------------\n");
        dense_free(dense);  
    }else{
        printf(RED"dense    : failed for allocation\n"RESET);
    }

    //coo
    coo_t *coo;
    // convert loaded matrix into COO format
    failed = coo_from_mm_buffer(&coo, mm_file); 
    if(!failed){
        printf("coo  : %.2f bytes\n", (double)coo_getMemSize(coo) / (double)mm_file->nnz);
        printf("-----------------------------------------\n");
        coo_free(coo);  
    }else{
        printf(RED"coo    : failed for allocation\n"RESET);
    }

    //csr
    csr_t *csr;
    // convert loaded matrix into CSR format
    failed = csr_from_mm_buffer(&csr, mm_file); 
    if(!failed){
        printf("csr  : %.2f bytes\n", (double)csr_getMemSize(csr) / (double)mm_file->nnz);
        printf("-----------------------------------------\n");
        csr_free(csr);  
    }else{
        printf(RED"csr    : failed for allocation\n"RESET);
    }

    //ell
    ell_t *ell;
    // convert loaded matrix into ELL format
    failed = ell_from_mm_buffer(&ell, mm_file);
    if(!failed){
        printf("ell  : %.2f bytes\n", (double)ell_getMemSize(ell) / (double)mm_file->nnz);
        printf("-----------------------------------------\n");
        ell_free(ell);  
    }else{
        printf(RED"ell    : failed for allocation\n"RESET);
    }

    //dia
    dia_t *dia;
    // convert loaded matrix into DIA format
    failed = dia_from_mm_buffer(&dia, mm_file);
    if(!failed){
        printf("dia    : %.2f bytes\n", (double)dia_getMemSize(dia) / (double)mm_file->nnz);
        printf("-----------------------------------------\n");
        dia_free(dia);  
    }else{
        printf(RED"dia    : failed for allocation\n"RESET);
    }
    
    freemm(mm_file);

    return 0;
}