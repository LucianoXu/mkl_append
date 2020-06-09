#include "mkl_append.h"
#include <iostream>
using namespace std;


int main() {
	lapack_complex_double* a = (lapack_complex_double*)mkl_malloc(9 * sizeof(lapack_complex_double), 64);
	a[0] = { 1,0 };
	a[3] = { 2,0 };
	a[4] = { 3,0 };
	a[6] = { 4,0 };
	a[7] = { 5,0 };
	a[8] = { 6,0 };

	lapack_complex_double* b = (lapack_complex_double*)mkl_calloc(9, sizeof(lapack_complex_double), 64);
	b[0] = { 1,0 };
	b[3] = { 2,0 };
	b[4] = { 3,0 };
	b[6] = { 4,0 };
	b[7] = { 5,0 };
	b[8] = { 6,0 };

	lapack_complex_double* c = (lapack_complex_double*)mkl_calloc(81,sizeof(lapack_complex_double), 64);

	zhelrmk(a, 3, 3, b, 3, 3, c, 9);
	zgem_out(c, 9, 9, 9);
	return 0;
}