#ifndef __MKL_APPEND_H__
#define __MKL_APPEND_H__
#include <mkl.h>
#include <iostream>

inline lapack_complex_double conj(lapack_complex_double a);
inline lapack_complex_double operator+(lapack_complex_double a, lapack_complex_double b);
inline lapack_complex_double operator-(lapack_complex_double a, lapack_complex_double b);
inline lapack_complex_double operator*(lapack_complex_double a, lapack_complex_double b);

inline std::ostream & operator <<(std::ostream& os, lapack_complex_double z);

//output the matrix(vector) G(n*m) with ostream
void zgem_out(lapack_complex_double* G, int m, int n);

//Hermitan matrix, packed storage, add
void zhpma(const lapack_complex_double* H1, const lapack_complex_double* H2, int n, lapack_complex_double* Hr);

//Compute Hr(mn*mn):=H1(m*m)(DirectProduct)H2(n*n)
//packed storage, L, row major, a[j+i*(i+1)/2]->A[i][j]
void zhplrmk(const lapack_complex_double* H1, int m, const lapack_complex_double* H2, int n, lapack_complex_double* Hr);

//DirectProduct for general matrixs(vectors)
void zgemk(const lapack_complex_double* G1, int m, int n, lapack_complex_double* G2, int p, int q, lapack_complex_double* R);

//Compute Hr(mn*mn):=H1(m*m)(DirectProduct)H2(n*n)
//full storage, L, row major,
//H1,H2 should be stored in lower triangle, row major
void zhelrmk(const lapack_complex_double* H1, int m, const lapack_complex_double* H2, int n, lapack_complex_double* Hr);

#endif