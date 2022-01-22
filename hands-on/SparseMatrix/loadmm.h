#pragma once

typedef struct mm_item
{
	int row;
	int col;
	double value;	
} mm_item_t;


typedef struct mm_file
{
	int nrow;
	int ncol;
	int nnz;

	mm_item_t *data;
	int data_size;
} mm_file_t;

mm_file_t *loadmm(const char *);
void freemm(mm_file_t *);