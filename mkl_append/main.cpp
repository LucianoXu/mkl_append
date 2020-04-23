#include "mkl_append.h"
#include <iostream>
using namespace std;
int main() {
	lapack_complex_double* a = (lapack_complex_double*)mkl_calloc(6, sizeof(lapack_complex_double), 64);
	a[0] = { 1.,0. };
	a[1] = { 3.,0. };
	a[3] = { 3.,0. };
	a[4] = { 2.,0. };
	zgem_out(a, 2, 2, 3);

	lapack_complex_double* b = (lapack_complex_double*)mkl_calloc(6, sizeof(lapack_complex_double), 64);
	b[0] = { 1.,0. };
	b[1] = { 3.,0. };
	b[3] = { 3.,0. };
	b[4] = { 2.,0. };
	zgem_out(b, 2, 2, 3);

	lapack_complex_double* c = (lapack_complex_double*)mkl_calloc(20, sizeof(lapack_complex_double), 64);
	zhelrmk(a, 2, 3, b, 2, 3, c, 5);
	zgem_out(c, 4, 4, 5);
	return 0;
}