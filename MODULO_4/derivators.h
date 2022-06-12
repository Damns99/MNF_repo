#ifndef DERIVATORS_H
#define DERIVATORS_H

#include <math.h>

#include "bound_cond_vecs.h"
#include "myfft.h"

namespace derivators {
	
	typedef bound_cond_vecs::BoundCondVec<double> (*Derivator)(const bound_cond_vecs::BoundCondVec<double>&, const double);
	
	namespace utilities {
		
		bound_cond_vecs::BoundCondVec<double> fwd_der_function(const bound_cond_vecs::BoundCondVec<double>& input, const double dx) {
			bound_cond_vecs::BoundCondVec<double> output(input.len(), input.getMode());
			for(int ii = 0; ii < input.len(); ii++) output[ii] = (input[ii + 1] - input[ii]) / dx;
			return output;
		}
		
		bound_cond_vecs::BoundCondVec<double> symm_der_function(const bound_cond_vecs::BoundCondVec<double>& input, const double dx) {
			bound_cond_vecs::BoundCondVec<double> output(input.len(), input.getMode());
			for(int ii = 0; ii < input.len(); ii++) output[ii] = (input[ii + 1] - input[ii - 1]) / (2. * dx);
			return output;
		}
		
		bound_cond_vecs::BoundCondVec<double> fft_der_function(const bound_cond_vecs::BoundCondVec<double>& input, const double dx) {
			int N = input.len();
			int M = N / 2 + 1;
			double areal[M], aimag[M], outputp[N];
			myFFT::my_dft_r2c_1d(N, input.data(), areal, aimag);
			for(int ii = 0; ii < M; ii++) {
				double ff = 2. * M_PI * ii / (N * dx);
				areal[ii] *= ff;
				aimag[ii] *= - ff;
			}
			if (N % 2 == 0) {
				areal[M-1] = 0.;
				aimag[M-1] = 0.;
			}
			myFFT::my_dft_c2r_1d(N, aimag, areal, outputp);
			bound_cond_vecs::BoundCondVec<double> output(N, outputp, input.getMode());
			return output;
		}
		
		bound_cond_vecs::BoundCondVec<double> fwd_der_i_function(const bound_cond_vecs::BoundCondVec<double>& input, const double dx) {
			bound_cond_vecs::BoundCondVec<double> output(input.len(), input.getMode());
			for(int ii = 0; ii < input.len(); ii++) output[ii] = (- 1 * input[ii + 2] + 4 * input[ii + 1] - 3 * input[ii]) / (2 * dx);
			return output;
		}
		
		
		bound_cond_vecs::BoundCondVec<double> fwd_der_2_function(const bound_cond_vecs::BoundCondVec<double>& input, const double dx) {
			bound_cond_vecs::BoundCondVec<double> output(input.len(), input.getMode());
			for(int ii = 0; ii < input.len(); ii++) output[ii] = (- 1 * input[ii + 3] + 4 * input[ii + 2] - 5 * input[ii + 1] + 2 * input[ii]) / (dx * dx);
			return output;
		}
		
		bound_cond_vecs::BoundCondVec<double> symm_der_2_function(const bound_cond_vecs::BoundCondVec<double>& input, const double dx) {
			bound_cond_vecs::BoundCondVec<double> output(input.len(), input.getMode());
			for(int ii = 0; ii < input.len(); ii++) output[ii] = (1 * input[ii + 1] - 2 * input[ii] + 1 * input[ii - 1]) / (dx * dx);
			return output;
		}
		
		bound_cond_vecs::BoundCondVec<double> fft_der_2_function(const bound_cond_vecs::BoundCondVec<double>& input, const double dx) {
			int N = input.len();
			int M = N / 2 + 1;
			double areal[M], aimag[M], outputp[N];
			myFFT::my_dft_r2c_1d(N, input.data(), areal, aimag);
			for(int ii = 0; ii < M; ii++) {
				double ff = 2. * M_PI * ii / (N * dx);
				areal[ii] *= - ff * ff;
				aimag[ii] *= - ff * ff;
			}
			if (N % 2 == 0) {
				areal[M-1] = 0.;
				aimag[M-1] = 0.;
			}
			myFFT::my_dft_c2r_1d(N, areal, aimag, outputp);
			bound_cond_vecs::BoundCondVec<double> output(N, outputp, input.getMode());
			return output;
		}
	}
	
	Derivator fwd_derive = &(utilities::fwd_der_function);
	Derivator symm_derive = &(utilities::symm_der_function);
	Derivator fft_derive = &(utilities::fft_der_function);
	Derivator fwd_derive_i = &(utilities::fwd_der_i_function);
	Derivator fwd_derive_2 = &(utilities::fwd_der_2_function);
	Derivator symm_derive_2 = &(utilities::symm_der_2_function);
	Derivator fft_derive_2 = &(utilities::fft_der_2_function);
	
}

#endif