#pragma once

typedef struct csr
{
	int m;
	int nnz;

	int *rowoffs;
	int *colind;
	double *val;

} csr_t;

void csr_alloc(csr_t **, int, int);
void csr_init(csr_t **, int, int, int*, int *, double *);
void csr_free(csr_t *);
void csr_spmv(csr_t *, double*, double*);
int csr_from_mm(csr_t **, const char *);
int csr_from_mm_buffer(csr_t **, const mm_file_t *);
unsigned int csr_getMemSize(csr_t *);