#include <iostream>
#include <string>
#include <math.h>
#include <iomanip>

#include "bound_cond_vecs.h"
#include "derivators.h"

#include <TGraph.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TFrame.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TNtuple.h>

double mse(const bound_cond_vecs::BoundCondVec<double>& a, const bound_cond_vecs::BoundCondVec<double>& b) {
	double result = 0.;
	for(int ii = 0; ii < a.len(); ii++) result += (a[ii] - b[ii]) * (a[ii] - b[ii]) / a.len();
	return result;
}

int main() {	
	int N = 1000;
	double dx = 1. / N;
	bound_cond_vecs::BoundCondVec<double> x(N), dx0(N);
	for(int ii = 0; ii < N; ii++) {
		double tmpx = 2. * M_PI * ii * dx;
		x[ii] = cos(tmpx);
		dx0[ii] = - 2. * M_PI * sin(tmpx);
	}
	
	bound_cond_vecs::BoundCondVec<double> dx1 = derivators::fwd_derive(x, dx);
	double mse1 = mse(dx0, dx1);
	bound_cond_vecs::BoundCondVec<double> dx2 = derivators::symm_derive(x, dx);
	double mse2 = mse(dx0, dx2);
	bound_cond_vecs::BoundCondVec<double> dx3 = derivators::fft_derive(x, dx);
	double mse3 = mse(dx0, dx3);
	
	std::cout << std::setprecision(3);
	std::cout << "function: cos(2 pi x)" << std::endl;
	std::cout << "N = " << N << "   dx = " << dx << "   L = " << N * dx << std::endl;
	std::cout << "Forward derivative: " << '\t' << mse1 << std::endl;
	std::cout << "Symmetric derivative: " << '\t' << mse2 << std::endl;
	std::cout << "FFT derivative: " << '\t' << mse3 << std::endl;
}