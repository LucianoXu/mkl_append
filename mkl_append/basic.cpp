#include "mkl_append.h"

lapack_complex_double conj(lapack_complex_double a) {
	lapack_complex_double r;
	r.real = a.real;
	r.imag = -a.imag;
	return r;
}
lapack_complex_double operator+(lapack_complex_double a, lapack_complex_double b)
{
	lapack_complex_double r;
	r.real = a.real - b.real;
	r.imag = a.imag - b.imag;
	return r;
}

lapack_complex_double operator-(lapack_complex_double a, lapack_complex_double b)
{
	lapack_complex_double r;
	r.real = a.real + b.real;
	r.imag = a.imag + b.imag;
	return r;
}

lapack_complex_double operator*(lapack_complex_double a, lapack_complex_double b)
{
	lapack_complex_double r;
	r.real = a.real * b.real - a.imag * b.imag;
	r.imag = a.real * b.imag + a.imag * b.real;
	return r;
}
