#ifndef CUDA_LATTICE_H
#define CUDA_LATTICE_H

#include "lattice.h"

#define cuda_index1d(I, L) (((I) >= 0) ? ((I) % (L)) : ((L) - (((L) - (I)) % (L))))
#define cuda_index2d(I, J, L) ((L) * cuda_index1d((I), (L)) + cuda_index1d((J), (L)))

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

constexpr int BLOCK_SIDE = 32; // so that 32 * 32 = 1024
constexpr int MAX_THREADS = 1024;
const dim3 thread_block(BLOCK_SIDE, BLOCK_SIDE);
dim3 grid;
constexpr int shared_size = (2*BLOCK_SIDE + 2) * (BLOCK_SIDE + 2);

__device__ int d_length;
__device__ double d_beta;
__device__ double d_extrafield;

__device__ int d_spin[MAX_LENGTH];
__device__ double d_rr[MAX_LENGTH];
__device__ double d_energy;
__device__ double d_magnetization;

__device__ double d_delta_energy[MAX_LENGTH];
__device__ double d_delta_magnetization[MAX_LENGTH];

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

__global__ void calculateEnergyMagnetizatonGPU();

__global__ void squareMetropolisStepGPU(int, int);

__global__ void sumDeltaReduction(int);

unsigned int nextPow2(unsigned int x) {
	x--;
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;
	return x++;
}

#endif