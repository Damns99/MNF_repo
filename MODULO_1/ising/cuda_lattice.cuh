#ifndef CUDA_LATTICE_H
#define CUDA_LATTICE_H

#include "lattice.h"

class CudaLattice2D {
	public:
		int length;
		int links_per_spin;
		double beta;
		double extrafield;
		int nbw;
		int nthreads;
		int nblocks;
		ran2::RandomGenerator gen;
		
		int* d_spin;
		int* d_links;
		double* d_rr;
		double* d_energy;
		double* d_magnetization;
		
		__host__ CudaLattice2D();
		__host__ ~CudaLattice2D();
		
		__host__ void cudaInitFromLattice(Lattice2D* h_lattice);
		
		__host__ void cudaMeasureEnergyMagnetization(double* h_energy, double* h_magnetization);
		
		__host__ void cudaUpdateMetropolis();
		
		__host__ void cudaRetrieveLattice(Lattice2D* h_lattice);
		
		__host__ void cudaDestroyLattice();
};

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, std::string file, int line, bool abort=false)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file.c_str(), line);
      if (abort) exit(code);
   }
}

__global__ void calculateEnergyMagnetizatonGPU(CudaLattice2D* lattice);

__global__ void squareMetropolisStepGPU(CudaLattice2D* lattice, int bw, int dead_spin);

#endif