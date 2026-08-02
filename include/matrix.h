#ifndef LAB_1_MATRIX_H
#define LAB_1_MATRIX_H

#include <stdio.h>

typedef struct
{
	int n, m;
	double *elements;
} Matrix;

void InitializeMatrix(Matrix *matrix, int n, int m);

int ReadMatrix(FILE *in, Matrix *matrix);

void ClearMatrix(Matrix *matrix);

#endif	  // LAB_1_MATRIX_H
