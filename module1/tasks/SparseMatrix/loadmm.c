#include <stdio.h>
#include <stdlib.h>
#include "loadmm.h"
#include "mmio.h"
#include "cmp.h"

mm_file_t *loadmm(const char *filename)
{
	mm_file_t *mm_file;

	MM_typecode matcode;
	FILE *f;

	int is_symmetric;
	int rounded_size;

	if ((f = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "[ERROR]: File %s doesn't exists\n", filename);
		fflush(stderr);
		abort();
	}

	if (mm_read_banner(f, &matcode) != 0) {
		fprintf(stderr, "[ERROR]: Could not process Matrix Market banner.\n");
		fflush(stderr);
		abort();
	}

	if (mm_is_complex(matcode) && mm_is_matrix(matcode) && mm_is_sparse(matcode)) {
		fprintf(stderr, "[ERROR]: Current IO doesn't support Matrix Market type: [%s]\n", mm_typecode_to_str(matcode));
		fflush(stderr);
		abort();;
	}

	//printf("matcode = %d->%s\n", matcode, mm_typecode_to_str(matcode));

	if (mm_is_symmetric(matcode))
		is_symmetric = 1;
	else
		is_symmetric = 0;

	mm_file = malloc(sizeof(mm_file_t));

	if (mm_read_mtx_crd_size(f, &mm_file->nrow, &mm_file->ncol, &mm_file->nnz) != 0) {
		fprintf(stderr, "[ERROR]: wrong mm_read_mtx_crd_size\n");
		fflush(stderr);
		abort();
	}

	if (is_symmetric)
		mm_file->data_size = mm_file->nnz * 2;
	else
		mm_file->data_size = mm_file->nnz;

	mm_file->data = malloc(mm_file->data_size * sizeof(mm_item_t));

	for (int i = 0; i < mm_file->nnz; i++) {
		if (fscanf(f, "%d %d %lf\n", &mm_file->data[i].row, &mm_file->data[i].col, &mm_file->data[i].value) != 3) {
			fprintf(stderr, "[ERROR]: fscanf failed for file %s and item no. %i\n", filename, i);
			fflush(stderr);
			abort();
		}

		//zero-indexing based
		mm_file->data[i].row--;
		mm_file->data[i].col--;

		if (is_symmetric && mm_file->data[i].row < mm_file->data[i].col) {
			fprintf(stderr, "[ERROR]: Bad symmetric value at %d\n", i);
			fflush(stderr);
			abort();
		}


		if (is_symmetric && mm_file->data[i].row != mm_file->data[i].col) {
			i++;
			mm_file->data[i].row = mm_file->data[i - 1].col;
			mm_file->data[i].col = mm_file->data[i - 1].row;
			mm_file->data[i].value = mm_file->data[i - 1].value;
			mm_file->nnz++;
		}

	}

	if (is_symmetric) {
		qsort(mm_file->data, mm_file->nnz, sizeof(mm_item_t), cmp_col_major);
	}

	return mm_file;

}

void freemm(mm_file_t *mm_file) {
	if (mm_file != NULL)
		free(mm_file->data);

	free(mm_file);
}

