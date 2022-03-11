#include <fftw3.h>
#include <math.h>
#include "myfft.h"

void myFFT::my_dft_r2c_1d(const unsigned int & N, double* in, double* outreal, double* outimag) {
    unsigned int M = N / 2 + 1;
	double norm = std::sqrt(N);
    fftw_plan p;
    double* in_ = fftw_alloc_real(2*M);
    fftw_complex* out_ = fftw_alloc_complex(2*M);
    p = fftw_plan_dft_r2c_1d(N, in_, out_, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
    for (unsigned int i = 0; i < N; i++) {
        in_[i] = in[i];
    }
    fftw_execute(p);
    for (unsigned int j = 0; j < M; j++) {
        outreal[j] = out_[j][0] / norm;
        outimag[j] = out_[j][1] / norm;
    }
    fftw_destroy_plan(p);
    fftw_free(in_);
    fftw_free(out_);
}

void myFFT::my_dft_c2r_1d(const unsigned int & N, double* inreal, double* inimag, double* out) {
    unsigned int M = N / 2 + 1;
	double norm = std::sqrt(N);
    fftw_plan p;
    fftw_complex* in_ = fftw_alloc_complex(2*M);
    double* out_ = fftw_alloc_real(2*M);
    p = fftw_plan_dft_c2r_1d(N, in_, out_, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
    for (unsigned int i = 0; i < M; i++) {
        in_[i][0] = inreal[i];
        in_[i][1] = inimag[i];
    }
    fftw_execute(p);
    for (unsigned int j = 0; j < N; j++) {
        out[j] = out_[j] / norm;
    }
    fftw_destroy_plan(p);
    fftw_free(in_);
    fftw_free(out_);
}

void myFFT::my_dft_c2c_1d(const unsigned int & N, double* inreal, double* inimag, double* outreal, double* outimag, int flag) {
    fftw_plan p;
    fftw_complex* in_ = fftw_alloc_complex(N);
    fftw_complex* out_ = fftw_alloc_complex(N);
    p = fftw_plan_dft_1d(N, in_, out_, flag, FFTW_ESTIMATE);
    for (unsigned int i = 0; i < N; i++) {
        in_[i][0] = inreal[i];
        in_[i][1] = inimag[i];
    }
    fftw_execute(p);
    for (unsigned int j = 0; j < N; j++) {
        outreal[j] = out_[j][0];
        outimag[j] = out_[j][1];
    }
    fftw_destroy_plan(p);
    fftw_free(in_);
    fftw_free(out_);
}