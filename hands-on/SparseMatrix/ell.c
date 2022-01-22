#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "loadmm.h"
#include "ell.h"
#include <string.h>
#include "cmp.h"


void ell_alloc(ell_t **ell, int m, int max_row){
	
	*ell = calloc(1, sizeof(ell_t));

	(*ell)->m = m;
	(*ell)->max_row = max_row;

	(*ell)->colind = calloc(m * max_row, sizeof(int));
	(*ell)->val = calloc(m * max_row, sizeof(double));

}

void ell_init(ell_t **ell, int m, int max_row, int *colind, double *val){
	*ell = calloc(1, sizeof(ell_t));

	(*ell)->m = m;
	(*ell)->max_row = max_row;

	(*ell)->colind = calloc(m * max_row, sizeof(int));
	(*ell)->val = calloc(m * max_row, sizeof(double));

	memcpy((*ell)->colind, colind, m * max_row * sizeof(int));
	memcpy((*ell)->val, val, m * max_row * sizeof(double));	
}


void ell_free(ell_t *ell){
	free(ell->colind);
	free(ell->val);
	free(ell);
}

void ell_spmv(ell_t *a, double *x, double *y){
	//column major format
	for(int i = 0; i < a->m; i++){
		for(int j = 0; j < a->max_row; j++){
			const int idx = i + j * a->m;
			y[i] += a->val[idx] * x[a->colind[idx]];
		}
	}
}

int ell_from_mm(ell_t **ell, const char *filename){

	mm_file_t *mm_file = loadmm(filename);

	qsort(mm_file->data, mm_file->nnz, sizeof(mm_item_t), cmp_row_major);

	int *rowoffs = calloc((mm_file->nrow + 1), sizeof(int));
	
	int max_row = 0;
	
	for(int i = 0; i < mm_file->nnz; i++)
	{
		rowoffs[mm_file->data[i].row + 1]++;
		max_row = rowoffs[mm_file->data[i].row + 1] > max_row ? rowoffs[mm_file->data[i].row + 1] : max_row;
	}

	for (int i = 0; i < mm_file->nrow; i++)
	{
		rowoffs[i + 1] += rowoffs[i];	
	}

    if((unsigned long long)mm_file->nrow * (unsigned long long)max_row > (int)((unsigned int)~0 >> 1)){
    	fprintf(stderr, "[ELL]: this matrix is out of range for ELL format\n");
	    freemm(mm_file);
	    free(rowoffs);
    	return 1;
    }	

	ell_alloc(ell, mm_file->nrow, max_row);

	for(int i = 0; i < mm_file->nrow; i++)
	{
		for (int j = rowoffs[i]; j < rowoffs[i + 1]; j++)
		{
			(*ell)->colind[i + (j - rowoffs[i]) * mm_file->nrow] = mm_file->data[j].col;			
			(*ell)->val[i + (j - rowoffs[i]) * mm_file->nrow] = mm_file->data[j].value;
		}
	}

	freemm(mm_file);
	free(rowoffs);

	return 0;

}

int ell_from_mm_buffer(ell_t **ell, const mm_file_t *mm_file){
	qsort(mm_file->data, mm_file->nnz, sizeof(mm_item_t), cmp_row_major);

	int *rowoffs = calloc((mm_file->nrow + 1), sizeof(int));
	
	int max_row = 0;
	
	for(int i = 0; i < mm_file->nnz; i++)
	{
		rowoffs[mm_file->data[i].row + 1]++;
		max_row = rowoffs[mm_file->data[i].row + 1] > max_row ? rowoffs[mm_file->data[i].row + 1] : max_row;
	}

	for (int i = 0; i < mm_file->nrow; i++)
	{
		rowoffs[i + 1] += rowoffs[i];	
	}

    if((unsigned long long)mm_file->nrow * (unsigned long long)max_row > (int)((unsigned int)~0 >> 1)){
    	fprintf(stderr, "[ELL]: this matrix is out of range for ELL format\n");
	    freemm(mm_file);
	    free(rowoffs);
    	return 1;
    }	

	ell_alloc(ell, mm_file->nrow, max_row);

	for(int i = 0; i < mm_file->nrow; i++)
	{
		for (int j = rowoffs[i]; j < rowoffs[i + 1]; j++)
		{
			(*ell)->colind[i + (j - rowoffs[i]) * mm_file->nrow] = mm_file->data[j].col;			
			(*ell)->val[i + (j - rowoffs[i]) * mm_file->nrow] = mm_file->data[j].value;
		}
	}

	free(rowoffs);

	return 0;
}

unsigned int ell_getMemSize(ell_t *ell){
	return ell->max_row * ell->m * (sizeof(int) + sizeof(double));
}