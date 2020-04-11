#ifndef __MKL_APPEND_H__
#define __MKL_APPEND_H__
#include <mkl.h>

lapack_complex_double conj(lapack_complex_double a);
lapack_complex_double operator+(lapack_complex_double a, lapack_complex_double b);
lapack_complex_double operator-(lapack_complex_double a, lapack_complex_double b);
lapack_complex_double operator*(lapack_complex_double a, lapack_complex_double b);

void zpma(const lapack_complex_double* H1, const lapack_complex_double* H2, int n, lapack_complex_double* Hr);

void zplmk(const lapack_complex_double* H1, int m, const lapack_complex_double* H2, int n, lapack_complex_double* Hr);

#endif