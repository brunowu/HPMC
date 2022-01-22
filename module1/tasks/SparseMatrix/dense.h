#pragma once

typedef struct dense
{
	int m;
	double *val;

} dense_t;

void dense_alloc(dense_t **, int);
void dense_init(dense_t **, int, double *);
void dense_free(dense_t *);
void dense_spmv(dense_t *, double*, double*);
int dense_from_mm(dense_t **, const char *);
int dense_from_mm_buffer(dense_t **, const mm_file_t *);
unsigned int dense_getMemSize(dense_t *);