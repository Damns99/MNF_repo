#ifndef CUDA_LATTICE_H
#define CUDA_LATTICE_H

#include "lattice.h"

__device__ int d_length;
__device__ int d_links_per_spin;
__device__ double d_beta;
__device__ double d_extrafield;

int nthreads;
int nblocks;

__device__ int d_spin[MAX_LENGTH];
__device__ int d_links[MAX_LENGTH];
__device__ double d_rr[MAX_LENGTH];
__device__ double d_energy;
__device__ double d_magnetization;

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

#endif