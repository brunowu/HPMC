#pragma once

typedef struct dia
{
	int m;

	int *diag;
	int ndiag;
	double *val;

} dia_t;

void dia_alloc(dia_t **, int, int);
void dia_init(dia_t **, int, int, int*, double*);
void dia_free(dia_t *);
void dia_spmv(dia_t *, double*, double*);
int dia_from_mm(dia_t **, const char *);
int dia_from_mm_buffer(dia_t **, const mm_file_t *);
unsigned int dia_getMemSize(dia_t *);