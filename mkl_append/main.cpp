#include "mkl_append.h"
#include <iostream>
using namespace std;
int main() {
	Zgem_I i;
	i.Prepare_I(3);
	zgem_out(i.Get_I(), 2, 2, i.Get_N());
	return 0;
}