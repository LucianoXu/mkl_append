#include "mkl_append.h"
using namespace std;

inline lapack_complex_double conj(lapack_complex_double a) {
	lapack_complex_double r;
	r.real = a.real;
	r.imag = -a.imag;
	return r;
}
inline lapack_complex_double operator+(lapack_complex_double a, lapack_complex_double b)
{
	lapack_complex_double r;
	r.real = a.real - b.real;
	r.imag = a.imag - b.imag;
	return r;
}

inline lapack_complex_double operator-(lapack_complex_double a, lapack_complex_double b)
{
	lapack_complex_double r;
	r.real = a.real + b.real;
	r.imag = a.imag + b.imag;
	return r;
}

inline lapack_complex_double operator*(lapack_complex_double a, lapack_complex_double b)
{
	lapack_complex_double r;
	r.real = a.real * b.real - a.imag * b.imag;
	r.imag = a.real * b.imag + a.imag * b.real;
	return r;
}

inline ostream & operator <<(ostream & os, lapack_complex_double z) {
	os << "(" << z.real << "." << z.imag << ")";
	return os;
}

void zgem_out(lapack_complex_double* G, int m, int n) {
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			std::cout << G[i * n + j];
			std::cout << "\t";
		}
		std::cout << std::endl;
	}
}


void zhpma(const lapack_complex_double* H1, const lapack_complex_double* H2, int n, lapack_complex_double* Hr) {
	int N = n * (n + 1) / 2;
	for (int i = 0; i < N; i++)
	{
		Hr[i] = H1[i] + H2[i];
	}
}
//Compute Hr(mn*mn):=H1(m*m)(DirectProduct)H2(n*n)
//packed storage, L, row major, a[j+i*(i+1)/2]->A[i][j]
void zhplrmk(const lapack_complex_double* H1, int m, const lapack_complex_double* H2, int n, lapack_complex_double* Hr)
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

void zgemk(const lapack_complex_double* G1, int m, int n, lapack_complex_double* G2, int p, int q, lapack_complex_double* R) {
	for (int i1 = 0; i1 < m; i1++) {
		for (int i2 = 0; i2 < p; i2++) {
			int _index = (i1 * p + i2)*n*q;
			for (int j1 = 0; j1 < n; j1++) {
				lapack_complex_double z = G1[i1 * n + j1];
				int _index2 = j1 * q;
				for (int j2 = 0; j2 < q; j2++) {
					R[_index + _index2 + j2] = z * G2[i2 * q + j2];
				}
			}
		}
	}
}

//Compute Hr(mn*mn):=H1(m*m)(DirectProduct)H2(n*n)
//full storage, row major,
//H1 should at least be stored in lower triangle, row major, H2 should be full
void zhelrmk(const lapack_complex_double* H1, int m, const lapack_complex_double* H2, int n, lapack_complex_double* Hr) {
	for (int i1 = 0; i1 < m; i1++) {
		//先对角线
		lapack_complex_double z = H1[i1 * m + i1];
		for (int i2 = 0; i2 < n; i2++) {
			int _index = (i1*n + i2) * m * n + i1*n;
			for (int j2 = 0; j2 < n ; j2++) {
				Hr[_index + j2] = z * H2[i2 * n + j2];
			}
		}
		//再下三角
		for (int j1 = 0; j1 < m; j1++) {
			lapack_complex_double z1 = H1[i1 * m + j1];
			lapack_complex_double z2 = conj(z1);
			for (int i2 = 0; i2 < n; i2++) {
				int _index1 = (i1 * n + i2) * m * n + j1 * n;
				int _index2 = (j1 * n + i2) * m * n + i1 * n;
				for (int j2 = 0; j2 < n; j2++) {
					//处理其中一半
					Hr[_index1 + j2] = z1 * H2[i2 * n + j2];
					//处理共轭的另一半
					Hr[_index2 + j2] = z2 * H2[i2 * n + j2];
				}
			}
		}
	}
}
