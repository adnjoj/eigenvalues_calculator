#include "include/algorithms.h"
#include "include/matrix.h"
#include "return_codes.h"

#include <stdio.h>

int main(int argc, char **argv)
{
	if (argc < 3)
	{
		fprintf(stderr, "Error: expected 2 arguments, got %i\n", argc);
		return ERROR_PARAMETER_INVALID;
	}

	int status_code = SUCCESS;

	FILE *in = fopen(argv[1], "r");
	FILE *out = NULL;

	if (in == NULL)
	{
		fprintf(stderr, "Error: failed to open file \"%s\"\n", argv[1]);
		return ERROR_CANNOT_OPEN_FILE;
	}

	int n;
	fscanf(in, "%i\n", &n);

	Matrix matrix, u_matrix, tmp1_matrix, tmp2_matrix;
	InitializeMatrix(&matrix, n, n);
	InitializeMatrix(&u_matrix, n, n);
	InitializeMatrix(&tmp1_matrix, n, n);
	InitializeMatrix(&tmp2_matrix, n, n);

	if (matrix.elements == NULL || u_matrix.elements == NULL || tmp1_matrix.elements == NULL || tmp2_matrix.elements == NULL)
	{
		fprintf(stderr, "Error: failed to allocate enough memory\n");
		status_code = ERROR_OUT_OF_MEMORY;
		goto CleanUp;
	}

	if (ReadMatrix(in, &matrix) != 0)
	{
		fprintf(stderr, "Error: failed to read matrix, input file has incorrect format\n");
		status_code = ERROR_DATA_INVALID;
		goto CleanUp;
	}

	fclose(in);
	in = NULL;

	out = fopen(argv[2], "w");

	if (out == NULL)
	{
		fprintf(stderr, "Error: failed to open file \"%s\"\n", argv[2]);
		status_code = ERROR_CANNOT_OPEN_FILE;
		goto CleanUp;
	}

	ToHessenbergForm(&matrix, &u_matrix, &tmp1_matrix, &tmp2_matrix);
	PrintEigenvalues(out, &matrix, &u_matrix, &tmp1_matrix, &tmp2_matrix);

CleanUp:
	if (in != NULL)
		fclose(in);
	if (out != NULL)
		fclose(out);

	ClearMatrix(&matrix);
	ClearMatrix(&u_matrix);
	ClearMatrix(&tmp1_matrix);
	ClearMatrix(&tmp2_matrix);

	return status_code;
}
