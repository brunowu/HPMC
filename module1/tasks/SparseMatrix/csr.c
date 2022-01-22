#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "loadmm.h"
#include "csr.h"
#include <string.h>
#include "cmp.h"
#if defined(USE_OPENMP)
#include <omp.h>
#endif

void csr_alloc(csr_t **csr, int m, int nnz){
	
	*csr = calloc(1, sizeof(csr_t));

	(*csr)->m = m;
	(*csr)->nnz = nnz;

	(*csr)->rowoffs = calloc( (m + 1 + nnz), sizeof(int));
	(*csr)->colind = &((*csr)->rowoffs[m+1]);	
	(*csr)->val = calloc(nnz, sizeof(double));

}


void csr_init(csr_t **csr, int m, int nnz, int* rowoffs, int *colind, double *val){
	*csr = calloc(1, sizeof(csr_t));

	(*csr)->m = m;
	(*csr)->nnz = nnz;

	(*csr)->rowoffs = calloc( (m + 1 + nnz), sizeof(int));
	(*csr)->colind = &((*csr)->rowoffs[m+1]);	
	(*csr)->val = calloc(nnz, sizeof(double));

	memcpy((*csr)->rowoffs, rowoffs, (m + 1) * sizeof(int));
	memcpy((*csr)->colind, colind, nnz * sizeof(int));
	memcpy((*csr)->val, val, nnz * sizeof(double));

}

void csr_free(csr_t *csr){
	free(csr->val);
	free(csr->rowoffs);
	free(csr);
}

void csr_spmv(csr_t *a, double *x, double *y){
#pragma omp parallel for
	for(int i = 0; i < a->m; i++){
		for(int j = a->rowoffs[i]; j < a->rowoffs[i+1]; j++){
			y[i] += a->val[j] * x[a->colind[j]];
		}
	}
}

int csr_from_mm(csr_t **csr, const char *filename){

	mm_file_t *mm_file = loadmm(filename);
    
    if(((unsigned long long)mm_file->nnz + (unsigned long long)mm_file->nrow) > (int)((unsigned int)~0 >> 1)){
    	fprintf(stderr, "[CSR]: this matrix is out of range for CSR format\n");
    	return 1;
    }

	qsort(mm_file->data, mm_file->nnz, sizeof(mm_item_t), cmp_row_major);

	csr_alloc(csr, mm_file->nrow, mm_file->nnz);

	for(int i = 0; i < mm_file->nnz; i++){
		(*csr)->val[i] =  mm_file->data[i].value;
		(*csr)->colind[i] = mm_file->data[i].col;
		(*csr)->rowoffs[mm_file->data[i].row + 1]++;
	}


	for(int i = 0; i < mm_file->nrow; i++){
		(*csr)->rowoffs[i + 1] += (*csr)->rowoffs[i];
	}

	freemm(mm_file);
	return 0;
}


int csr_from_mm_buffer(csr_t **csr, const mm_file_t *mm_file){

	qsort(mm_file->data, mm_file->nnz, sizeof(mm_item_t), cmp_row_major);

    if(((unsigned long long)mm_file->nnz + (unsigned long long)mm_file->nrow) > (int)((unsigned int)~0 >> 1)){
    	fprintf(stderr, "[CSR]: this matrix is out of range for CSR format\n");
    	return 1;
    }

	csr_alloc(csr, mm_file->nrow, mm_file->nnz);

	for(int i = 0; i < mm_file->nnz; i++){
		(*csr)->val[i] =  mm_file->data[i].value;
		(*csr)->colind[i] = mm_file->data[i].col;
		(*csr)->rowoffs[mm_file->data[i].row + 1]++;
	}


	for(int i = 0; i < mm_file->nrow; i++){
		(*csr)->rowoffs[i + 1] += (*csr)->rowoffs[i];
	}

	return 0;

}

unsigned int csr_getMemSize(csr_t *csr){
	return (csr->m + 1 + csr->nnz) * sizeof(int) + csr->nnz * sizeof(double);
}
