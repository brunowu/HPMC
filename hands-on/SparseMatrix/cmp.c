#include <stdio.h>
#include <stdlib.h>
#include "loadmm.h"

int cmp_row_major(const void * a, const void * b) {
	int comp = ((mm_item_t*) a)->row - ((mm_item_t*) b)->row;

	if(comp < 0)
		return -1;

	if(comp > 0)
		return 1;

	return (((mm_item_t*) a)->col - ((mm_item_t*) b)->col);
}


int cmp_col_major(const void * a, const void * b) {
    int comp = ((mm_item_t*) a)->col - ((mm_item_t*) b)->col;

    if(comp < 0)
        return -1;

    if(comp > 0)
        return 1;

    return (((mm_item_t*) a)->row - ((mm_item_t*) b)->row);
}


int cmp_diag_major(const void * a, const void * b) {

	int off_diag_a = ((mm_item_t*) a)->col - ((mm_item_t*) a)->row;
	int off_diag_b = ((mm_item_t*) b)->col - ((mm_item_t*) b)->row;

    int comp = off_diag_a - off_diag_b;

    if(comp < 0)
        return -1;

    if(comp > 0)
        return 1;

    return (((mm_item_t*) a)->row - ((mm_item_t*) b)->row);
}

