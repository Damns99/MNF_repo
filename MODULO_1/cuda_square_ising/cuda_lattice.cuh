#ifndef CUDA_LATTICE_H
#define CUDA_LATTICE_H

#include "lattice.h"

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
constexpr int shared_size = (BLOCK_SIDE + 2) * (BLOCK_SIDE + 2);

__device__ int d_length;
__device__ double d_beta;
__device__ double d_extrafield;

__device__ int d_spin[MAX_LENGTH];
__device__ double d_rr[MAX_LENGTH];
__device__ double d_energy;
__device__ double d_magnetization;

__device__ int d_delta_energy[MAX_LENGTH];
__device__ int d_delta_magnetization[MAX_LENGTH];

/* double delta_energy[MAX_LENGTH];
double delta_magnetization[MAX_LENGTH]; */

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

__global__ void squareMetropolisStepGPU(int bw, int dead_spin);

unsigned int nextPow2(unsigned int x) {
	x--;
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;
	return x++;
}



// 1) atomicAdd
__global__ void squareMetropolisStepGPU1(int bw, int dead_spin, int length_);

__host__ void cudaUpdateMetropolis1();

__host__ void cudaMeasureEnergyMagnetization1();

// 2) delta + GPU
__global__ void squareMetropolisStepGPU2(int bw, int dead_spin, int length_);

__global__ void sumDeltaReduction2(int n);

__host__ void cudaUpdateMetropolis2();

__host__ void cudaMeasureEnergyMagnetization2();

// 3) delta + CPU
__host__ void sumDeltaSerial3();

__host__ void cudaUpdateMetropolis3();

__host__ void cudaMeasureEnergyMagnetization3();

// 4) old

__global__ void calculateEnergyMagnetizatonGPU4(int length_);

__global__ void squareMetropolisStepGPU4(int bw, int dead_spin, int length_);

__host__ void cudaUpdateMetropolis4();

__host__ void cudaMeasureEnergyMagnetization4();

#endif