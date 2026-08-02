#include "../include/matrix.h"
#include "../include/matrix_slice.h"

#include <math.h>
#include <stdio.h>

/*
 * Calculates Frobenius norm of given matrix
 */
double FrobeniusNorm(Matrix *matrix)
{
	double result = 0;

	for (int i = 0; i < matrix->n; i++)
	{
		for (int j = 0; j < matrix->m; j++)
		{
			result += matrix->elements[i * matrix->m + j] * matrix->elements[i * matrix->m + j];
		}
	}

	return sqrt(result);
}

/*
 * Function computes the Householder vector which will be stored into "out" argument
 * Col argument is treated as matrix slice with only one column
 */
void GetU(MatrixSlice *col, MatrixSlice *out)
{
	double norm = GetMatrixSliceNorm(col);
	double sign = MSGet(col, 0, 0) >= 0 ? 1 : -1;

	FillMatrixSlice(out, 0);

	if (norm == 0)
	{
		return;
	}

	MSSet(out, 0, 0, sign * norm);
	AddMultipliedMatrixSlice(out, col, 1);
	NormalizeMatrixSlice(out);
}

/*
 * Applies Householder transformation of the given matrix[n x m] using given Householder vector[n x 1]
 * tmp1[m x m], tmp2[m x m] are used for intermediate results
 */
void HouseholderTransformation(MatrixSlice *a_slice, MatrixSlice *u, MatrixSlice *tmp1, MatrixSlice *tmp2)
{
	MultiplyMatrixSlices(u, a_slice, tmp1, 1, 0);
	MultiplyMatrixSlices(u, tmp1, tmp2, 0, 0);
	AddMultipliedMatrixSlice(a_slice, tmp2, -2);

	int t = a_slice->from_row;
	a_slice->from_row = a_slice->from_col;
	a_slice->from_col = t;

	t = a_slice->rows_count;
	a_slice->rows_count = a_slice->cols_count;
	a_slice->cols_count = t;

	MultiplyMatrixSlices(a_slice, u, tmp1, 0, 0);
	MultiplyMatrixSlices(tmp1, u, tmp2, 0, 1);
	AddMultipliedMatrixSlice(a_slice, tmp2, -2);
}

/*
 * Function transforms the given matrix into Hessenberg form
 */
void ToHessenbergForm(Matrix *matrix, Matrix *u_matrix, Matrix *tmp1_matrix, Matrix *tmp2_matrix)
{
	MatrixSlice tmp1 = { 0, 0, 0, 0, tmp1_matrix };
	MatrixSlice tmp2 = { 0, 0, 0, 0, tmp2_matrix };

	for (int i = 0; i < matrix->n - 2; i++)
	{
		MatrixSlice col = { i + 1, i, matrix->n - i - 1, 1, matrix };
		MatrixSlice u = { 0, 0, col.rows_count, 1, u_matrix };
		MatrixSlice a_slice = { i + 1, 0, matrix->n - i - 1, matrix->n, matrix };

		GetU(&col, &u);

		HouseholderTransformation(&a_slice, &u, &tmp1, &tmp2);
	}
}

/*
 * Function calculates coefficients for polynomial associated with double shift
 * Matrix slice is supposed to be correct
 */
void GetDoubleShiftQRPolynomial(MatrixSlice *h, double *b, double *c)
{
	int n = h->rows_count;

	double trace = MSGet(h, n - 2, n - 2) + MSGet(h, n - 1, n - 1);
	double det = MSGet(h, n - 2, n - 2) * MSGet(h, n - 1, n - 1) - MSGet(h, n - 2, n - 1) * MSGet(h, n - 1, n - 2);

	if (trace * trace > det * 4)
	{
		double v1 = (trace + sqrt(trace * trace - 4 * det)) / 2;
		double v2 = (trace - sqrt(trace * trace - 4 * det)) / 2;

		if (fabs(v1 - MSGet(h, n - 1, n - 1)) < fabs(v2 - MSGet(h, n - 1, n - 1)))
		{
			v2 = v1;
		}
		else
		{
			v1 = v2;
		}

		(*b) = -v1 - v2;
		(*c) = v1 * v2;
	}
	else
	{
		(*b) = -trace;
		(*c) = det;
	}
}

void DoubleShiftQRStep(MatrixSlice *h, Matrix *u_matrix, Matrix *tmp1_matrix, Matrix *tmp2_matrix)
{
	int n = h->rows_count;

	double b, c;
	GetDoubleShiftQRPolynomial(h, &b, &c);

	double col_elements[3] = {
		MSGet(h, 0, 0) * MSGet(h, 0, 0) + MSGet(h, 0, 1) * MSGet(h, 1, 0),
		MSGet(h, 1, 0) * MSGet(h, 0, 0) + MSGet(h, 1, 1) * MSGet(h, 1, 0),
		MSGet(h, 2, 0) * MSGet(h, 0, 0) + MSGet(h, 2, 1) * MSGet(h, 1, 0)
	};
	col_elements[0] += MSGet(h, 0, 0) * b + c;
	col_elements[1] += MSGet(h, 1, 0) * b;

	Matrix col_matrix = { 3, 1, col_elements };

	MatrixSlice col_slice = { 0, 0, 3, 1, &col_matrix };
	MatrixSlice u_slice = { 0, 0, 3, 1, u_matrix };
	MatrixSlice tmp1_slice = { 0, 0, 0, 0, tmp1_matrix };
	MatrixSlice tmp2_slice = { 0, 0, 0, 0, tmp2_matrix };
	MatrixSlice h_slice = { h->from_row, h->from_col, 3, n, h->matrix };

	GetU(&col_slice, &u_slice);
	HouseholderTransformation(&h_slice, &u_slice, &tmp1_slice, &tmp2_slice);

	for (int j = 1; j < n - 1; j++)
	{
		int k = j + 3 > n ? n : j + 3;

		MatrixSlice h_col = { h->from_row + j, h->from_col + j - 1, k - j, 1, h->matrix };

		u_slice.rows_count = k - j;

		h_slice.from_row = h->from_row + j;
		h_slice.from_col = h->from_col;
		h_slice.rows_count = k - j;
		h_slice.cols_count = n;

		GetU(&h_col, &u_slice);
		HouseholderTransformation(&h_slice, &u_slice, &tmp1_slice, &tmp2_slice);

		MSSet(h, k - 1, j - 1, 0);
	}
}

void DoubleShiftQR(Matrix *hess, Matrix *u_matrix, Matrix *tmp1_matrix, Matrix *tmp2_matrix)
{
	int n = hess->n;

	double tolerance = FrobeniusNorm(hess) * 1e-8;

	int start = 0;
	int end = n;

	while (end - start > 2)
	{
		if (fabs(hess->elements[(end - 1) * n + (end - 2)]) < tolerance)
		{
			hess->elements[(end - 1) * n + (end - 2)] = 0;
			end -= 1;
		}
		else if (fabs(hess->elements[(start + 1) * n + start]) < tolerance)
		{
			hess->elements[(start + 1) * n + start] = 0;
			start += 1;
		}
		else if (fabs(hess->elements[(end - 2) * n + (end - 3)]) < tolerance)
		{
			hess->elements[(end - 2) * n + (end - 3)] = 0;
			end -= 2;
		}
		else if (fabs(hess->elements[(start + 2) * n + start + 1]) < tolerance)
		{
			hess->elements[(start + 2) * n + start + 1] = 0;
			start += 2;
		}
		else
		{
			MatrixSlice h_sm = { start, start, end - start, end - start, hess };
			DoubleShiftQRStep(&h_sm, u_matrix, tmp1_matrix, tmp2_matrix);
		}
	}
}

/*
 * Calculates eigenvalues and prints them into out
 */
void PrintEigenvalues(FILE *out, Matrix *h, Matrix *u_matrix, Matrix *tmp1_matrix, Matrix *tmp2_matrix)
{
	int n = h->n;

	DoubleShiftQR(h, u_matrix, tmp1_matrix, tmp2_matrix);

	int i = 0;
	while (i < n)
	{
		if (i == n - 1 || h->elements[(i + 1) * n + i] == 0)
		{
			fprintf(out, "%f\n", h->elements[i * n + i]);
			i += 1;
		}
		else
		{
			double a = (h->elements[i * n + i] + h->elements[(i + 1) * n + i + 1]) / 2;
			double b = h->elements[i * n + i] * h->elements[(i + 1) * n + i + 1] -
					   h->elements[i * n + i + 1] * h->elements[(i + 1) * n + i];

			if (a * a >= b)
			{
				double c = sqrt(a * a - b);
				fprintf(out, "%g\n%g\n", a + c, a - c);
			}
			else
			{
				double c = sqrt(b - a * a);
				fprintf(out, "%g +%gi\n%g -%gi\n", a, c, a, c);
			}

			i += 2;
		}
	}
}
