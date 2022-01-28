#pragma once

/**
 * A structure to represent sparse matrix in DIA format
 */
typedef struct dia
{
	int m; //size of matrix

	int *diag; // a pointer storing the offsets of each non-null diagonal
	int ndiag; // number of diagonals with at least one non-zero element
	double *val; //a pointer storing the non-zeros elements

} dia_t;

//allocation of sparse matrix in DIA format, with empty entries
void dia_alloc(dia_t **, int, int);

//initialization of sparse matrix in DIA format with given entries
void dia_init(dia_t **, int, int, int*, double*);

//free
void dia_free(dia_t *);

//SpMV with DIA format
void dia_spmv(dia_t *, double*, double*);

//load a matrix from a given matrix market file into DIA format
int dia_from_mm(dia_t **, const char *);

//load a matrix from a loaded buffer of a given matrix market file into DIA format
int dia_from_mm_buffer(dia_t **, const mm_file_t *);

//get the bytes of memory required per non-zeros elements
unsigned int dia_getMemSize(dia_t *);
