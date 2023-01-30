#ifndef CUDA_LATTICE_H
#define CUDA_LATTICE_H

#include "lattice.h"

#define cuda_index1d(I, L) (((I) >= 0) ? ((I) % (L)) : ((L) - (((L) - (I)) % (L))))

inline unsigned int iDivUp(const unsigned int &a, const unsigned int &b) {
	return (a%b != 0) ? (a/b+1) : (a/b);
}

__device__ double myAtomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val + __longlong_as_double(assumed)));
    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);
    return __longlong_as_double(old);
}

constexpr int BLOCK_SIDE = 1024; // so that 1024 = 1024
constexpr int MAX_THREADS = 1024;
const dim3 thread_block(BLOCK_SIDE);
dim3 grid;

__device__ int d_length;
__device__ int d_nparticles;
__device__ int d_p_length;
__device__ double d_beta;

__device__ double d_y[MAX_LENGTH];
__device__ int d_links[MAX_LENGTH];
__device__ double d_rr[MAX_LENGTH];

/* __device__ Rule d_rules[MAX_RULES];
__device__ int d_repetitions[MAX_RULES];
__device__ int d_nrules; */

__device__ double d_obs1[MAX_PARTICLES];
__device__ double d_obs2[MAX_PARTICLES];
__device__ double d_delta_obs1[MAX_LENGTH];
__device__ double d_delta_obs2[MAX_LENGTH];

__device__ double d_pars[MAX_PARS];

__host__ void cudaInitFromLattice();

__host__ void cudaMeasureEnergyMagnetization();

__host__ void cudaUpdateMetropolis();

__host__ void cudaRetrieveLattice();

__host__ void cudaDestroyLattice();

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, std::string file, int line, bool abort=false)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file.c_str(), line);
      if (abort) exit(code);
   }
}

__global__ void calculateObsGPU();

__global__ void metropolisStepGPU(int, int, int, Function_ds);

unsigned int nextPow2(unsigned int x) {
	x--;
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;
	return x++;
}

typedef __device__ double Function_ds(double, double, double, double, double*, double*);

__device__ double harm_pot(double y0, double y1, double y2, double yp, double* pars, double outputs[2]) {
	double tmp0 = (yp - y0), tmp1 = (yp + y0), tmp2 = (y1 + y2);
	double ds = tmp0 * (tmp1 - tmp2) / pars[0] + tmp0 * tmp1 * pars[0] / 2;
	outputs[0] = tmp0 * tmp1;
	outputs[1] = 2. * (tmp0 * (tmp1 - tmp2));
	return ds;
}

__device__ double double_hole(double y0, double y1, double y2, double yp, double* pars, double outputs[2]) {
	double tmp0 = (yp - y0), tmp1 = (yp + y0), tmp2 = (y1 + y2);
	double tmp3 = (yp * yp + y0 * y0 - 2.);
	double ds = tmp0 * (tmp1 - tmp2) / pars[0] + tmp0 * tmp1 * tmp3 * pars[0] * pars[1];
	outputs[0] = tmp0 * tmp1;
	outputs[1] = 2. * (tmp0 * (tmp1 - tmp2));
	return ds;
}

Function_ds func_ds = nullptr;



#endif