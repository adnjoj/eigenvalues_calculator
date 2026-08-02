#ifndef LAB_1_ALGORITHMS_H
#define LAB_1_ALGORITHMS_H

#include "matrix.h"

#include <stdio.h>

void ToHessenbergForm(Matrix *matrix, Matrix *u_matrix, Matrix *tmp1_matrix, Matrix *tmp2_matrix);

void PrintEigenvalues(FILE *out, Matrix *h, Matrix *u_matrix, Matrix *tmp1_matrix, Matrix *tmp2_matrix);

#endif	  // LAB_1_ALGORITHMS_H
