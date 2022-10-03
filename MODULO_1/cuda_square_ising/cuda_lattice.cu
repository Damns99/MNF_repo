#include "cuda_lattice.cuh"
#include <cuda.h>
#include <stdio.h>

// sum d_delta_energy
__global__ void calculateEnergyMagnetizatonGPU() {
	d_magnetization = 0.;
	for (int i = 0; i < d_length * d_length; i++) {
        int curr_spin = d_spin[i];
		d_energy += d_delta_energy[i] / (1. * d_length * d_length);
        d_delta_energy[i] = 0.;
		d_energy += - curr_spin * d_extrafield / (1. * d_length * d_length);
		d_magnetization += curr_spin / (1. * d_length * d_length);
	}
}//..........................

// the new one takes half the time
__global__ void calculateEnergyMagnetizatonGPU_old() {
	d_energy = 0.;
	d_magnetization = 0.;
	for (int i = 0; i < d_length; i++) for (int j = 0; j < d_length; j++) {
        int curr_spin = d_spin[index2d(i,j,d_length)];
		int neighbours = d_spin[index2d(i+1,j,d_length)] + d_spin[index2d(i,j+1,d_length)];
		d_energy += - curr_spin * neighbours;
		d_magnetization += curr_spin;
	}
	double normalization = 1. / (d_length * d_length);
	d_magnetization *= normalization;
	d_energy *= JACC * normalization;
	d_energy += - d_extrafield * magnetization;
}

__global__ void squareMetropolisStepGPU(int bw, int dead_spin) {
    extern __shared__ int s_sub_spin[];
    
    int index = threadIdx.x + blockIdx.x * blockDim.x;
	if(index >= d_length * d_length / 2) return;
    int idx = 2 * index + ((2 * index / d_length + bw) % 2);
    int idx2 = idx - ((2 * blockIdx.x * blockDim.x) / d_length) * d_length;
	
    int x = idx2 % d_length, y = idx2 / d_length;
    //printf("th = %d, bl = %d, index = %d, idx = %d, idx2 = %d, x = %d, y = %d\n", threadIdx.x, blockIdx.x, index, idx, idx2, x, y);
    //printf("to fill: %d %d %d %d\n", (x + 1) + (d_length + 2) * (y + 1), (x) + (d_length + 2) * (y + 1), (x + 1) + (d_length + 2) * (y), (x + 2) + (d_length + 2) * (y + 1));
    //printf("to use: %d %d %d %d %d\n", (x + 1) + (d_length + 2) * (y + 1), (x) + (d_length + 2) * (y + 1), (x + 2) + (d_length + 2) * (y + 1), (x + 1) + (d_length + 2) * (y), (x + 1) + (d_length + 2) * (y + 2));
    s_sub_spin[(x + 1) + (d_length + 2) * (y + 1)] = d_spin[idx]; //centre
    s_sub_spin[(x) + (d_length + 2) * (y + 1)] = d_spin[d_links[d_links_per_spin * idx + 2]]; //left
    if(threadIdx.x < d_length) s_sub_spin[(x + 1) + (d_length + 2) * (y)] = d_spin[d_links[d_links_per_spin * idx + 0]]; //up
    if(threadIdx.x >= blockDim.x - d_length - 1) s_sub_spin[(x + 1) + (d_length + 2) * (y + 2)] = d_spin[d_links[d_links_per_spin * idx + 1]]; //down
    if(x % d_length >= d_length - 2) s_sub_spin[(x + 2) + (d_length + 2) * (y + 1)] = d_spin[d_links[d_links_per_spin * idx + 3]]; //right
    
    __syncthreads();
    
	if(idx == dead_spin) return; // avoid loops!
	int newspin = - s_sub_spin[(x + 1) + (d_length + 2) * (y + 1)];
	double denergy = -2. * newspin * d_extrafield;
	denergy += -2. * JACC * newspin * s_sub_spin[(x) + (d_length + 2) * (y + 1)];
	denergy += -2. * JACC * newspin * s_sub_spin[(x + 2) + (d_length + 2) * (y + 1)];
	denergy += -2. * JACC * newspin * s_sub_spin[(x + 1) + (d_length + 2) * (y)];
	denergy += -2. * JACC * newspin * s_sub_spin[(x + 1) + (d_length + 2) * (y + 2)];
	double r = exp(- d_beta * denergy);
	if (d_rr[2 * index + bw] < r) {
        d_spin[idx] = newspin;
        d_delta_energy[idx] += denergy;
    }
}

// old but takes almost same time...
__global__ void squareMetropolisStepGPU_old(int bw, int dead_spin) {    
    int index = threadIdx.x + blockIdx.x * blockDim.x;
	if(index >= d_length * d_length / 2) return;
    int idx = 2 * index + ((2 * index / d_length + bw) % 2);    
	if(idx == dead_spin) return; // avoid loops!
	int newspin = - d_spin[idx];
	double denergy = -2. * newspin * d_extrafield;
    for (int j = 0; j < d_links_per_spin; j++) {
        d_energy += -2. * JACC * newspin * d_spin[d_links[d_links_per_spin * idx + j]] / (1. * d_length * d_length);
    }
	double r = exp(- d_beta * denergy);
	if (d_rr[2 * index + bw] < r) {
        d_spin[idx] = newspin;
        d_delta_energy[idx] += denergy;
    }
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
    
    double* delta_energy;
    gpuErrchk(cudaGetSymbolAddress((void **)&delta_energy, d_delta_energy));
    gpuErrchk(cudaMemset(delta_energy, 0., sizeof(double) * sizeof(spin) / sizeof(int)));
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
	int shared_size = (length + 2) * ((2 * nthreads - 1) / length + 3);
    
	squareMetropolisStepGPU<<<nblocks, nthreads, shared_size * sizeof(int)>>>(0, dead_spin);
	gpuErrchk(cudaDeviceSynchronize());
	squareMetropolisStepGPU<<<nblocks, nthreads, shared_size * sizeof(int)>>>(1, dead_spin);
	gpuErrchk(cudaDeviceSynchronize());
}