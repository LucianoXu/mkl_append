#ifndef __MKL_APPEND_H__
#define __MKL_APPEND_H__
#include "mkl.h"

extern inline lapack_complex_double conj(lapack_complex_double a);
extern inline lapack_complex_double operator+(lapack_complex_double a, lapack_complex_double b);
extern inline lapack_complex_double operator-(lapack_complex_double a, lapack_complex_double b);
extern inline lapack_complex_double operator*(lapack_complex_double a, lapack_complex_double b);

extern void zpma(const lapack_complex_double* H1, const lapack_complex_double* H2, int n, lapack_complex_double* Hr);

extern void zplmk(const lapack_complex_double* H1, int m, const lapack_complex_double* H2, int n, lapack_complex_double* Hr);

#endif