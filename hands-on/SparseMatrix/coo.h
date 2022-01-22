#pragma once

typedef struct coo
{
	int m;
	int nnz;

	int *rowind;
	int *colind;
	double *val;

} coo_t;

void coo_alloc(coo_t **, int, int);
void coo_init(coo_t **, int, int, int*, int*, double *);
void coo_free(coo_t *);
void coo_spmv(coo_t *, double*, double*);
int coo_from_mm(coo_t **, const char *);
int coo_from_mm_buffer(coo_t **, const mm_file_t *);
unsigned int coo_getMemSize(coo_t *);