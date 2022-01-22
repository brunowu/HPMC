#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "loadmm.h"
#include "dia.h"
#include <string.h>
#include "cmp.h"

#define max(x,y) (((x) >= (y)) ? (x) : (y))
#define min(x,y) (((x) <= (y)) ? (x) : (y))

void dia_alloc(dia_t **dia, int m, int ndiag){
	
	*dia = calloc(1, sizeof(dia_t));

	(*dia)->m = m;
	(*dia)->ndiag = ndiag;

	(*dia)->diag = calloc(ndiag, sizeof(int));
	(*dia)->val = calloc(m * ndiag, sizeof(double));

}

void dia_init(dia_t **dia, int m, int ndiag, int *diag, double *val){
	
	*dia = calloc(1, sizeof(dia_t));

	(*dia)->m = m;
	(*dia)->ndiag = ndiag;

	(*dia)->diag = calloc(ndiag, sizeof(int));
	(*dia)->val = calloc(m * ndiag, sizeof(double));

	memcpy((*dia)->diag, diag, ndiag * sizeof(int));
	memcpy((*dia)->val, val, m * ndiag * sizeof(double));		

}

void dia_free(dia_t *dia){
	free(dia->diag);
	free(dia->val);
	free(dia);
}

void dia_spmv(dia_t *a, double *x, double *y){

	for(int i = 0; i < a->m; i++){
		for(int j = 0; j < a->ndiag; j++){
			int start = max(-a->diag[j], 0);
			int end = a->diag[j] > 0 ?  a->m - a->diag[j] : a->m;
			if((i >= start) && (i < end)){
				y[i] += a->val[j * a->m + i] * x[i + a->diag[j]];
			}
		}
	}
}

int dia_from_mm(dia_t **dia, const char *filename){
	
	mm_file_t *mm_file = loadmm(filename);

	qsort(mm_file->data, mm_file->nnz, sizeof(mm_item_t), cmp_diag_major);

	//counting ndiag;
	int tmp = mm_file->data[0].col - mm_file->data[0].row;

	int *diag = calloc(mm_file->nrow, sizeof(int));	
	int *diag_nnz = calloc(mm_file->nrow, sizeof(int));

	int ndiag = 1;
	diag[0] = tmp;
	diag_nnz[0] = 1;

	for(int i = 1; i < mm_file->nnz; i++){
		if((mm_file->data[i].col - mm_file->data[i].row) != tmp)
		{
			ndiag++;
			tmp = mm_file->data[i].col - mm_file->data[i].row;
			diag[ndiag-1] = tmp;
		}
		diag_nnz[ndiag-1]++;

	}

    if((unsigned long long)mm_file->nrow * (unsigned long long)ndiag > (int)((unsigned int)~0 >> 1)){
    	
    	fprintf(stderr, "[DIA]: this matrix is out of range for DIA format\n");
  		free(diag);
		free(diag_nnz);
		freemm(mm_file);

    	return 1;

    }
    
	dia_alloc(dia, mm_file->nrow, ndiag);


	int cnt = 0;
	for(int i = 0; i < ndiag; i++){
		(*dia)->diag[i] = diag[i];
		for(int j = 0; j < diag_nnz[i]; j++){
			(*dia)->val[i * mm_file->nrow + mm_file->data[cnt].row] = mm_file->data[cnt].value;
			cnt++;
		}
	}

    free(diag);
	free(diag_nnz);

	return 0;

}


int dia_from_mm_buffer(dia_t **dia, const mm_file_t *mm_file){
	qsort(mm_file->data, mm_file->nnz, sizeof(mm_item_t), cmp_diag_major);

	//counting ndiag;
	int tmp = mm_file->data[0].col - mm_file->data[0].row;

    int total_ndiag = 2 * mm_file->nrow - 1;
	int *diag = calloc(total_ndiag, sizeof(int));	
	int *diag_nnz = calloc(total_ndiag, sizeof(int));

	int ndiag = 1;
	diag[0] = tmp;
	diag_nnz[0] = 1;

	for(int i = 1; i < mm_file->nnz; i++){
		if((mm_file->data[i].col - mm_file->data[i].row) != tmp)
		{
			ndiag++;
			tmp = mm_file->data[i].col - mm_file->data[i].row;
			diag[ndiag-1] = tmp;
		}
		diag_nnz[ndiag-1]++;

	}

    if((unsigned long long)mm_file->nrow * (unsigned long long)ndiag > (int)((unsigned int)~0 >> 1)){
    	fprintf(stderr, "[DIA]: this matrix is out of range for DIA format\n");
  		free(diag);
		free(diag_nnz);

    	return 1;

    }

	dia_alloc(dia, mm_file->nrow, ndiag);


	int cnt = 0;
	for(int i = 0; i < ndiag; i++){
		(*dia)->diag[i] = diag[i];
		for(int j = 0; j < diag_nnz[i]; j++){
			(*dia)->val[i * mm_file->nrow + mm_file->data[cnt].row] = mm_file->data[cnt].value;
			cnt++;
		}
	}

    free(diag);
	free(diag_nnz);

	return 0;
	

}



unsigned int dia_getMemSize(dia_t *dia){
	return dia->ndiag * sizeof(int) + dia->ndiag * dia->m * sizeof(double);
}
