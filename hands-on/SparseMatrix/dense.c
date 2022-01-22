#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "loadmm.h"
#include <string.h>
#include "dense.h"
#include "cmp.h"

void dense_alloc(dense_t **dense, int m){
	
	*dense = calloc(1, sizeof(dense_t));

	(*dense)->m = m;
	(*dense)->val = calloc(m * m, sizeof(double));

}

void dense_init(dense_t **dense, int m, double *val){
	*dense = calloc(1, sizeof(dense_t));

	(*dense)->m = m;
	(*dense)->val = calloc(m * m, sizeof(double));

	memcpy((*dense)->val, val, m * m * sizeof(double));
}


void dense_free(dense_t *dense){
	free(dense->val);
	free(dense);
}

void dense_spmv(dense_t *a, double *x, double *y){

	for(int i = 0; i < a->m; i++){
		for(int j = 0; j < a->m; j++){
			y[i] += a->val[i * a->m + j] * x[j];
		}
	}
}

int dense_from_mm(dense_t **dense, const char *filename){

	mm_file_t *mm_file = loadmm(filename);
    
    if((unsigned long long)mm_file->nrow * (unsigned long long)mm_file->nrow > (int)((unsigned int)~0 >> 1)){
    	fprintf(stderr, "[DENSE]: this matrix is out of range for DENSE format\n");
    	return 1;
    }
	dense_alloc(dense, mm_file->nrow);


	for (int i = 0; i < mm_file->nnz; i++) {
		int row = mm_file->data[i].row;
		int col = mm_file->data[i].col;
		(*dense)->val[row * (*dense)->m + col] = mm_file->data[i].value;
	}

	freemm(mm_file);
	return 0;

}

int dense_from_mm_buffer(dense_t **dense, const mm_file_t *mm_file){

    if((unsigned long long)mm_file->nrow * (unsigned long long)mm_file->nrow > (int)((unsigned int)~0 >> 1)){
    	fprintf(stderr, "[DENSE]: this matrix is out of range for DENSE format\n");
    	return 1;
    }

	dense_alloc(dense, mm_file->nrow);

	for (int i = 0; i < mm_file->nnz; i++) {
		int row = mm_file->data[i].row;
		int col = mm_file->data[i].col;
		(*dense)->val[row * (*dense)->m + col] = mm_file->data[i].value;
	}

	return 0;
}

unsigned int dense_getMemSize(dense_t *dense){

	return dense->m * dense->m * sizeof(double); 

}