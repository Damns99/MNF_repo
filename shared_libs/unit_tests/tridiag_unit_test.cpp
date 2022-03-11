#include <iostream>
#include "tridiag.h"

#define N 5
#define EPS 5e-16

int main() {
	std::cout << "Testing triDiag: ";
	
	double a[N] = {1, 1, 1, 1, 1};
	double b[N] = {4, 4, 4, 4, 4};
	double c[N] = {1, 1, 1, 1, 1};
	double d[N] = {4, 3, 11, 14, 3};
	double x[N];
	double pred_x[N] = {1, 0, 2, 3, 0};
	
	triDiag::triDiagSolve(N, a, b, c, d, x);
	
	bool result = true;
	for(unsigned ii = 0; ii < N; ii++) if((x[ii] - pred_x[ii]) * (x[ii] - pred_x[ii]) > EPS * EPS) result = false;
	
	std::string res = result ? "passed" : "failed";
	std::cout << res << std::endl << std::endl;
}