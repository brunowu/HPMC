#pragma once

/**
 * A structure to represent sparse matrix in COO format
 */
typedef struct coo
{
	int m; //size of matrix
	int nnz; //number of non-zero elements in the sparse matrix

	int *rowind; //a pointer storing the row indices of non-zeros elements
	int *colind; //a pointer storing the column indices of non-zeros elements
	double *val; //a pointer storing the non-zeros elements

} coo_t;

//allocation of sparse matrix in COO format, with empty entries
void coo_alloc(coo_t **, int, int);

//initialization of sparse matrix in COO format with given entries
void coo_init(coo_t **, int, int, int*, int*, double *);

//free
void coo_free(coo_t *);

//SpMV with COO format
void coo_spmv(coo_t *, double*, double*);

//load a matrix from a loaded buffer of a given matrix market file into COO format
int coo_from_mm(coo_t **, const char *);

//load a matrix from a loaded buffer of a given matrix market file into COO format
int coo_from_mm_buffer(coo_t **, const mm_file_t *);

//get the bytes of memory required per non-zeros elements
unsigned int coo_getMemSize(coo_t *);