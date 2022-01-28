#pragma once

/**
 * A structure to represent sparse matrix in ELL format
 */
typedef struct ell
{
	int m; //size of matrix

	int *colind; //a pointer storing the column indices of non-zeros elements
	double *val; //a pointer storing the non-zeros elements
	int max_row; // the maximum number of elements of all rows

} ell_t;

//allocation of sparse matrix in ELL format, with empty entries
void ell_alloc(ell_t **, int, int);

//initialization of sparse matrix in ELL format with given entries
void ell_init(ell_t **, int, int, int *, double*);

//free
void ell_free(ell_t *);

//SpMV with ELL format
void ell_spmv(ell_t *, double*, double*);

//load a matrix from a given matrix market file into ELL format
int ell_from_mm(ell_t **, const char *);

//load a matrix from a loaded buffer of a given matrix market file into ELL format
int ell_from_mm_buffer(ell_t **, const mm_file_t *);

//get the bytes of memory required per non-zeros elements
unsigned int ell_getMemSize(ell_t *);
