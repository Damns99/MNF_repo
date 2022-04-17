#include "cuda_lattice.cuh"
#include <cuda.h>
#include <stdio.h>

__global__ void calculateEnergyMagnetizatonGPU() {
	d_energy = 0.;
	d_magnetization = 0.;
	for (int i = 0; i < d_length * d_length; i++) {
		for (int j = 0; j < d_links_per_spin; j++) {
			d_energy += - JACC / 2. * d_spin[i] * d_spin[d_links[d_links_per_spin * i + j]] / (1. * d_length * d_length);
		}
		d_energy += - d_spin[i] * d_extrafield / (1. * d_length * d_length);
		d_magnetization += d_spin[i] / (1. * d_length * d_length);
	}
}

__global__ void squareMetropolisStepGPU(int bw, int dead_spin) {
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	if(index >= d_length * d_length / 2) return;
	int x = 2 * index + ((2 * index / d_length + bw) % 2);
	if(x == dead_spin) return; // avoid loops!
	int newspin = - d_spin[x];
	double denergy = -2. * newspin * d_extrafield;
	for (int j = 0; j < d_links_per_spin; j++) {
		denergy += -2. * JACC * newspin * d_spin[d_links[d_links_per_spin * x + j]];
	}
	double r = exp(- d_beta * denergy);
	if (d_rr[2 * index + bw] < r) d_spin[x] = newspin;
}

__host__ void cudaInitFromLattice() {
	assert(length % 2 == 0 && links_per_spin == 4);
    
	int nbw = length * length / 2;
	nthreads = nbw > 256 ? 256 : nbw;
	nblocks = (nbw - 1) / nthreads + 1;
	
	/* int device_number;
	gpuErrchk(cudaGetDevice(&device_number));
	std::cout << "device:" << device_number << std::endl; */
	
	gpuErrchk(cudaMemcpyToSymbol(d_length, &length, sizeof(int)));
	gpuErrchk(cudaMemcpyToSymbol(d_links_per_spin, &links_per_spin, sizeof(int)));
	gpuErrchk(cudaMemcpyToSymbol(d_beta, &beta, sizeof(double)));
	gpuErrchk(cudaMemcpyToSymbol(d_extrafield, &extrafield, sizeof(double)));
	gpuErrchk(cudaMemcpyToSymbol(d_energy, &energy, sizeof(double)));
	gpuErrchk(cudaMemcpyToSymbol(d_spin, spin, sizeof(spin)));
	gpuErrchk(cudaMemcpyToSymbol(d_links, links, sizeof(links)));
	gpuErrchk(cudaMemcpyToSymbol(d_energy, &energy, sizeof(double)));
	gpuErrchk(cudaMemcpyToSymbol(d_magnetization, &magnetization, sizeof(double)));
}

__host__ void cudaRetrieveLattice() {	
	gpuErrchk(cudaMemcpyFromSymbol(spin, d_spin, sizeof(spin)));
	gpuErrchk(cudaMemcpyFromSymbol(&energy, d_energy, sizeof(double)));
	gpuErrchk(cudaMemcpyFromSymbol(&magnetization, d_magnetization, sizeof(double)));
}

__host__ void cudaDestroyLattice() {
	cudaDeviceReset();
}

__host__ void cudaMeasureEnergyMagnetization() {
	calculateEnergyMagnetizatonGPU<<<1, 1>>>();
	gpuErrchk(cudaDeviceSynchronize());
	gpuErrchk(cudaMemcpyFromSymbol(&energy, d_energy, sizeof(double)));
	gpuErrchk(cudaMemcpyFromSymbol(&magnetization, d_magnetization, sizeof(double)));
}

__host__ void cudaUpdateMetropolis() {
	double rr[length * length];
	for(int ii = 0; ii < length * length; ii++) rr[ii] = gen.randF();
	gpuErrchk(cudaMemcpyToSymbol(d_rr, rr, sizeof(rr)));
	int dead_spin = gen.randL(0, length * length);
	
	squareMetropolisStepGPU<<<nblocks, nthreads>>>(0, dead_spin);
	gpuErrchk(cudaDeviceSynchronize());
	squareMetropolisStepGPU<<<nblocks, nthreads>>>(1, dead_spin);
	gpuErrchk(cudaDeviceSynchronize());
}