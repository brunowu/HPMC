#pragma once

/**
 * A structure to represent sparse matrix in CSR format
 */
typedef struct csr
{
	int m; //size of matrix
	int nnz; //number of non-zero elements in the sparse matrix

	int *rowoffs; //a pointer storing the offsets of pointer of each row
	int *colind; //a pointer storing the column indices of non-zeros elements
	double *val; //a pointer storing the non-zeros elements

} csr_t;

//allocation of sparse matrix in CSR format, with empty entries
void csr_alloc(csr_t **, int, int);

//initialization of sparse matrix in CSR format with given entries
void csr_init(csr_t **, int, int, int*, int *, double *);

//free
void csr_free(csr_t *);

//SpMV with CSR format
void csr_spmv(csr_t *, double*, double*);

//load a matrix from a given matrix market file into CSR format
int csr_from_mm(csr_t **, const char *);

//load a matrix from a loaded buffer of a given matrix market file into CSR format
int csr_from_mm_buffer(csr_t **, const mm_file_t *);

//get the bytes of memory required per non-zeros elements
unsigned int csr_getMemSize(csr_t *);