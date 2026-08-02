#ifndef LAB_1_MATRIX_SLICE_H
#define LAB_1_MATRIX_SLICE_H

#include "matrix.h"

typedef struct
{
	int from_row;
	int from_col;
	int rows_count;
	int cols_count;
	Matrix *matrix;
} MatrixSlice;

double MSGet(const MatrixSlice *sm, int i, int j);

void MSSet(MatrixSlice *sm, int i, int j, double val);

void FillMatrixSlice(MatrixSlice *slice, double val);

double GetMatrixSliceNorm(const MatrixSlice *slice);

void NormalizeMatrixSlice(MatrixSlice *sm);

void AddMultipliedMatrixSlice(MatrixSlice *a, const MatrixSlice *b, double multiple);

void MultiplyMatrixSlices(const MatrixSlice *a, const MatrixSlice *b, MatrixSlice *out, _Bool transpose_first, _Bool transpose_second);

#endif	  // LAB_1_MATRIX_SLICE_H
