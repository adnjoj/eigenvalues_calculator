#include "../include/matrix.h"

#include <stdio.h>
#include <stdlib.h>

void InitializeMatrix(Matrix *matrix, int n, int m)
{
	matrix->n = n;
	matrix->m = m;

	matrix->elements = malloc(sizeof(double) * n * m);
}

int ReadMatrix(FILE *in, Matrix *matrix)
{
	for (int i = 0; i < matrix->n * matrix->m; i++)
	{
		if (fscanf(in, "%lf", matrix->elements + i) != 1)
		{
			return 1;
		}
	}

	return 0;
}

void ClearMatrix(Matrix *matrix)
{
	if (matrix->elements != NULL)
	{
		free(matrix->elements);
	}
}
