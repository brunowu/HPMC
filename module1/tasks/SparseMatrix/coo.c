#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "loadmm.h"
#include <string.h>
#include "coo.h"
#include "cmp.h"

#if defined(USE_OPENMP)
#include <omp.h>
#endif

void coo_alloc(coo_t **coo, int m, int nnz){
	
	*coo = calloc(1, sizeof(coo_t));

	(*coo)->m = m;
	(*coo)->nnz = nnz;

	(*coo)->rowind = calloc(2 * nnz, sizeof(int));
	(*coo)->colind = &((*coo)->rowind[nnz]);
	(*coo)->val = calloc(nnz, sizeof(double));

}

void coo_init(coo_t **coo, int m, int nnz, int* rowind, int* colind, double *val){
	*coo = calloc(1, sizeof(coo_t));

	(*coo)->m = m;
	(*coo)->nnz = nnz;

	(*coo)->rowind = calloc(2 * nnz, sizeof(int));
	(*coo)->colind = &((*coo)->rowind[nnz]);
	(*coo)->val = calloc(nnz, sizeof(double));

	memcpy((*coo)->rowind, rowind, nnz * sizeof(int));
	memcpy((*coo)->colind, colind, nnz * sizeof(int));
	memcpy((*coo)->val, val, nnz * sizeof(double));

}


void coo_free(coo_t *coo){
	free(coo->val);
	free(coo->rowind);
	free(coo);
}

void coo_spmv(coo_t *a, double *x, double *y){
#pragma omp paralllel for
	for(int i =0; i < a->nnz; i++){
		y[a->rowind[i]] += a->val[i] * x[a->colind[i]];  
	}

}

int coo_from_mm(coo_t **coo, const char *filename){

	mm_file_t *mm_file = loadmm(filename);

    if((unsigned long long)(2 * mm_file->nnz) > (int)((unsigned int)~0 >> 1)){
    	fprintf(stderr, "[COO]: this matrix is out of range for COO format\n");
    	return 1;
    }

	assert(mm_file->nrow = mm_file->ncol);
	qsort(mm_file->data, mm_file->nnz, sizeof(mm_item_t), cmp_row_major);
	coo_alloc(coo, mm_file->nrow, mm_file->nnz);

	for (int i = 0; i < mm_file->nnz; i++) {
		(*coo)->rowind[i] = mm_file->data[i].row;
		(*coo)->colind[i] = mm_file->data[i].col;
		(*coo)->val[i] = mm_file->data[i].value;
	}

	freemm(mm_file);
	return 0;

}

int coo_from_mm_buffer(coo_t **coo, const mm_file_t *mm_file){

    if((unsigned long long)(2 * mm_file->nnz) > (int)((unsigned int)~0 >> 1)){
    	fprintf(stderr, "[COO]: this matrix is out of range for COO format\n");
    	return 1;
    }

	qsort(mm_file->data, mm_file->nnz, sizeof(mm_item_t), cmp_row_major);
	coo_alloc(coo, mm_file->nrow, mm_file->nnz);

	for (int i = 0; i < mm_file->nnz; i++) {
		(*coo)->rowind[i] = mm_file->data[i].row;
		(*coo)->colind[i] = mm_file->data[i].col;
		(*coo)->val[i] = mm_file->data[i].value;
	}

	return 0;
}

unsigned int coo_getMemSize(coo_t *coo){
	return 2 * coo->nnz * sizeof(int) + coo->nnz * sizeof(double); 
}
