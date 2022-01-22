#pragma once

typedef struct ell
{
	int m;

	int *colind;
	double *val;
	int max_row;

} ell_t;

void ell_alloc(ell_t **, int, int);
void ell_init(ell_t **, int, int, int *, double*);
void ell_free(ell_t *);
void ell_spmv(ell_t *, double*, double*);
int ell_from_mm(ell_t **, const char *);
int ell_from_mm_buffer(ell_t **, const mm_file_t *);
unsigned int ell_getMemSize(ell_t *);