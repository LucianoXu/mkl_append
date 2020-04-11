#include "mkl_append.h"

void zpma(const lapack_complex_double* H1, const lapack_complex_double* H2, int n, lapack_complex_double* Hr) {
	int N = n * (n + 1) / 2;
	for (int i = 0; i < N; i++)
	{
		Hr[i] = H1[i] + H2[i];
	}
}
