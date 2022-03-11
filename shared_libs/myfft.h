#ifndef MY_FFT_H
#define MY_FFT_H

#define MYFFT_EPS 1e-15

namespace myFFT {
	
	void my_dft_r2c_1d(const unsigned int & N, double* in, double* outreal, double* outimag);
	void my_dft_c2r_1d(const unsigned int & N, double* inreal, double* inimag, double* out);
	void my_dft_c2c_1d(const unsigned int & N, double* inreal, double* inimag, double* outreal, double* outimag, int flag);
	
}

#endif