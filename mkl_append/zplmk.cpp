#include "mkl_append.h"

//Compute Hr(mn*mn):=H1(m*m)(DirectProduct)H2(n*n)
//packed storage, L, row major, a[j+i*(i+1)/2]->A[i][j]
void zplmk(const lapack_complex_double* H1, int m, const lapack_complex_double* H2, int n, lapack_complex_double* Hr)
{
	for (int i1 = 0; i1 < m; i1++)
	{
		int n_ahead_H1 = i1 * (i1 + 1) / 2;
		for (int i2 = 0; i2 < n; i2++)
		{
			int n_ahead_H2 = i2 * (i2 + 1) / 2;
			int ki = i1 * n + i2;
			int n_ahead_Hr = ki * (ki + 1) / 2;
			for (int j1 = 0; j1 <= i1; j1++)
			{
				for (int j2 = 0; j2 <= i2; j2++)
				{
					int kj = j1 * n + j2;
					Hr[n_ahead_Hr + kj] = H1[n_ahead_H1 + j1] * H2[n_ahead_H2 + j2];
				}
			}
		}
	}
	for (int i1 = 1; i1 < m; i1++)
	{
		int n_ahead_H1 = i1 * (i1 + 1) / 2;
		for (int i2 = 0; i2 < n - 1; i2++)
		{
			int ki = i1 * n + i2;
			int n_ahead_Hr = ki * (ki + 1) / 2;
			for (int j1 = 0; j1 < i1; j1++)
			{
				for (int j2 = i2 + 1; j2 < n; j2++)
				{
					int kj = j1 * n + j2;
					Hr[n_ahead_Hr + kj] = H1[n_ahead_H1 + j1] * conj(H2[j2 * (j2 + 1) / 2 + i2]);
				}
			}
		}
	}
}
