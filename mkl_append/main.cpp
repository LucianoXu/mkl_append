#include "mkl_append.h"
#include <iostream>
using namespace std;
int main() {
	lapack_complex_double a = { 1.,-1. };
	cout << conj(a).imag;
	return 0;
}