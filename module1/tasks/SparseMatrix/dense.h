#pragma once

/**
 * A structure to represent dense matrix
 */
typedef struct dense
{
	int m; //size of matrix 
	double *val; //a pointer which stores the elements of matrix

} dense_t;

//allocation of dense matrix, with empty entries
void dense_alloc(dense_t **, int);

//initialization of dense matrix with given entries
void dense_init(dense_t **, int, double *);

//free
void dense_free(dense_t *);

//function which computes matrix-vector product
void dense_spmv(dense_t *, double*, double*);

//load a matrix from local matrix market file
int dense_from_mm(dense_t **, const char *);

//load a matrix from a loaded buffer of a given matrix market file
int dense_from_mm_buffer(dense_t **, const mm_file_t *);

//get the bytes of memory required per non-zeros elements
unsigned int dense_getMemSize(dense_t *);