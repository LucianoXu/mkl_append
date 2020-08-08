#include "mkl_append.h"
#include <iostream>
using namespace std;


int main() {
	int N;
	cin >> N;
	int* p = new int[N];
	for (int i = 0; i < N; i++) {
		cin >> p[i];
	}
	cout << Epsilon(N, p);
	return 0;
}