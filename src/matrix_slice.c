#include "../include/matrix_slice.h"

#include <math.h>

double MSGet(const MatrixSlice *slice, int i, int j)
{
	return slice->matrix->elements[(slice->from_row + i) * slice->matrix->m + slice->from_col + j];
}

void MSSet(MatrixSlice *slice, int i, int j, double val)
{
	slice->matrix->elements[(slice->from_row + i) * slice->matrix->m + slice->from_col + j] = val;
}

void FillMatrixSlice(MatrixSlice *slice, double val)
{
	for (int i = 0; i < slice->rows_count; i++)
	{
		for (int j = 0; j < slice->cols_count; j++)
		{
			MSSet(slice, i, j, val);
		}
	}
}

/*
 * "slice" argument is treated as matrix slice with only one column
 */
double GetMatrixSliceNorm(const MatrixSlice *slice)
{
	double result = 0;

	for (int i = 0; i < slice->rows_count; i++)
	{
		double el = MSGet(slice, i, 0);
		result += el * el;
	}

	return sqrt(result);
}

/*
 * "slice" argument is treated as matrix slice with only one column
 */
void NormalizeMatrixSlice(MatrixSlice *slice)
{
	double norm = GetMatrixSliceNorm(slice);

	if (norm == 0)
	{
		return;
	}

	for (int i = 0; i < slice->rows_count; i++)
	{
		MSSet(slice, i, 0, MSGet(slice, i, 0) / norm);
	}
}

void AddMultipliedMatrixSlice(MatrixSlice *a, const MatrixSlice *b, double multiple)
{
	for (int i = 0; i < a->rows_count; i++)
	{
		for (int j = 0; j < a->cols_count; j++)
		{
			MSSet(a, i, j, MSGet(a, i, j) + MSGet(b, i, j) * multiple);
		}
	}
}

void MultiplyMatrixSlices(const MatrixSlice *a, const MatrixSlice *b, MatrixSlice *out, _Bool transpose_first, _Bool transpose_second)
{
	int left_rows = transpose_first ? a->cols_count : a->rows_count;
	int left_cols = transpose_first ? a->rows_count : a->cols_count;
	int right_cols = transpose_second ? b->rows_count : b->cols_count;

	out->rows_count = left_rows;
	out->cols_count = right_cols;

	for (int i = 0; i < out->rows_count; i++)
	{
		for (int j = 0; j < out->cols_count; j++)
		{
			double result = 0;

			for (int k = 0; k < left_cols; k++)
			{
				double left_el = transpose_first ? MSGet(a, k, i) : MSGet(a, i, k);
				double right_el = transpose_second ? MSGet(b, j, k) : MSGet(b, k, j);
				result += left_el * right_el;
			}

			MSSet(out, i, j, result);
		}
	}
}
