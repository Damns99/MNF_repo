#include <iostream>
#include <vector>
#include <math.h>
#include "myfft.h"

#define N 1024
#define M N / 2 + 1
#define FREQ 1. / N * (N / 16)
#define OMEGA 2 * M_PI * FREQ

int main() {
	std::cout << "Testing myfft: ";
	
	std::vector<double> a, b, c, d;
	a.reserve(N);
	b.reserve(M);
	c.reserve(M);
	d.reserve(N);
	for(unsigned ii = 0; ii < N; ii++) a[ii] = std::sin(ii * OMEGA);
	myFFT::my_dft_r2c_1d(N, a.data(), b.data(), c.data());
	myFFT::my_dft_c2r_1d(N, b.data(), c.data(), d.data());
	
	bool result = true;
	for(unsigned ii = 0; ii < N; ii++) if(std::abs(a[ii] - d[ii]) > MYFFT_EPS) {
		std::cout << ii << ":\t" << a[ii] << "\t" << d[ii] << std::endl;
		result = false;
	}
	
	std::string res = result ? "passed" : "failed";
	std::cout << res << std::endl << std::endl;
}