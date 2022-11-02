#include "cuda_lattice.cuh"
#include <cuda.h>
#include <stdio.h>

__host__ void cudaInitFromLattice() {
	grid = dim3(iDivUp(length, 2 * BLOCK_SIDE), iDivUp(length, BLOCK_SIDE));
	
	gpuErrchk(cudaMemcpyToSymbol(d_length, &length, sizeof(int)));
	gpuErrchk(cudaMemcpyToSymbol(d_beta, &beta, sizeof(double)));
	gpuErrchk(cudaMemcpyToSymbol(d_extrafield, &extrafield, sizeof(double)));
	gpuErrchk(cudaMemcpyToSymbol(d_spin, spin, sizeof(spin)));
	gpuErrchk(cudaMemcpyToSymbol(d_energy, &energy, sizeof(double)));
	gpuErrchk(cudaMemcpyToSymbol(d_magnetization, &magnetization, sizeof(double)));
    
    double* temp_delta_energy;
    gpuErrchk(cudaGetSymbolAddress((void **)&temp_delta_energy, d_delta_energy));
    gpuErrchk(cudaMemset(temp_delta_energy, 0, sizeof(spin) * sizeof(double) / sizeof(int)));
    double* temp_delta_magnetization;
    gpuErrchk(cudaGetSymbolAddress((void **)&temp_delta_magnetization, d_delta_magnetization));
    gpuErrchk(cudaMemset(temp_delta_magnetization, 0, sizeof(spin) * sizeof(double) / sizeof(int)));
}

__host__ void cudaRetrieveLattice() {
	gpuErrchk(cudaMemcpyFromSymbol(spin, d_spin, sizeof(spin)));
	//gpuErrchk(cudaMemcpyFromSymbol(&energy, d_energy, sizeof(double)));
	//gpuErrchk(cudaMemcpyFromSymbol(&magnetization, d_magnetization, sizeof(double)));
}

__host__ void cudaDestroyLattice() {
	cudaDeviceReset();
}

// 2) delta + GPU
__global__ void squareMetropolisStepGPU(int bw, int dead_spin, int length_) {
    extern __shared__ int s_sub_spin[];
	
	// thread position
	int i = blockIdx.y * blockDim.y + threadIdx.y;
	int j = blockIdx.x * blockDim.x + threadIdx.x;
	// spin position
	int ii = i;
	int jj = 2 * j + ((i + bw) % 2);
	if((ii >= length_) || (jj >= length_)) return;
	// shared position
	int x = 2 * threadIdx.x + ((threadIdx.y + bw) % 2) + 1;
	int y = threadIdx.y + 1;
	
	s_sub_spin[(2*BLOCK_SIDE+2)*y+x] = d_spin[cuda_index2d(ii,jj,length_)]; //centre
	s_sub_spin[(2*BLOCK_SIDE+2)*y+(x-1)] = d_spin[cuda_index2d(ii,jj-1,length_)]; //left
    if(threadIdx.y == 0) s_sub_spin[(2*BLOCK_SIDE+2)*(y-1)+x] = d_spin[cuda_index2d(ii-1,jj,length_)]; //up
    if((threadIdx.y == blockDim.y - 1) || (ii == length_ - 1)) s_sub_spin[(2*BLOCK_SIDE+2)*(y+1)+x] = d_spin[cuda_index2d(ii+1,jj,length_)]; //down
    if((threadIdx.x == blockDim.x - 1) || (jj == length_ - 1)) s_sub_spin[(2*BLOCK_SIDE+2)*y+(x+1)] = d_spin[cuda_index2d(ii,jj+1,length_)]; //right
    __syncthreads();
	
	if(cuda_index2d(ii,jj,length_) == dead_spin) return; // avoid loops!
	int newspin = - s_sub_spin[(2*BLOCK_SIDE+2)*y+x];
	int neighbours = s_sub_spin[(2*BLOCK_SIDE+2)*(y-1)+x] + s_sub_spin[(2*BLOCK_SIDE+2)*(y+1)+x] + s_sub_spin[(2*BLOCK_SIDE+2)*y+(x-1)] + s_sub_spin[(2*BLOCK_SIDE+2)*y+(x+1)];
	double denergy = -2. * newspin * (JACC * neighbours + d_extrafield);
	
	double r = exp(- d_beta * denergy);
	if (d_rr[cuda_index2d(ii,jj,length_)] < r) {
        d_spin[cuda_index2d(ii,jj,length_)] = newspin;
		d_delta_energy[cuda_index2d(ii,jj,length_)] += denergy;
		d_delta_magnetization[cuda_index2d(ii,jj,length_)] += 2. * newspin;
    }
}

__global__ void sumDeltaReduction(int n) {
	extern __shared__ double data[];
	
	int tid = threadIdx.x;
	int index = threadIdx.x + (blockIdx.x * 2) * blockDim.x;
	
	//if(index < n) printf("block %d: adding d_d_m[%d] = %f to data[%d]\n", blockIdx.x, index, d_delta_magnetization[index], tid + blockDim.x);
	double esum = (index < n) ? d_delta_energy[index] : 0;
	double msum = (index < n) ? d_delta_magnetization[index] : 0;
	if(index + blockDim.x < n) {
		//printf("block %d: adding d_d_m[%d] = %f to data[%d]\n", blockIdx.x, index + blockDim.x, d_delta_magnetization[index + blockDim.x], tid + blockDim.x);
		esum += d_delta_energy[index + blockDim.x];
		msum += d_delta_magnetization[index + blockDim.x];
	}
	data[tid] = esum;
	data[tid + blockDim.x] = msum;
	__syncthreads();
	
	for(int s = blockDim.x / 2; s > 0; s /= 2) {
		if(tid < s) {
			//printf("block %d: adding data[%d] = %f to data[%d]\n", blockIdx.x, blockDim.x + tid + s, data[blockDim.x + tid + s], blockDim.x + tid);
			data[tid] += data[tid + s];
			data[blockDim.x + tid] += data[blockDim.x + tid + s];
		}
		__syncthreads();
	}
	
	if(tid == 0) {
		//printf("block %d: setting d_d_m[%d] to data[%d] = %f\n", blockIdx.x, blockIdx.x, blockDim.x, data[blockDim.x]);
		d_delta_energy[blockIdx.x] = data[0];
		d_delta_magnetization[blockIdx.x] = data[blockDim.x];
	}
}

__host__ void cudaUpdateMetropolis() {
	double rr[length * length];
	for(int ii = 0; ii < length * length; ii++) rr[ii] = gen.randF();
	gpuErrchk(cudaMemcpyToSymbol(d_rr, rr, sizeof(rr)));
	int dead_spin = gen.randL(0, length * length);
    
	squareMetropolisStepGPU<<<grid, thread_block, shared_size * sizeof(int)>>>(0, dead_spin, length);
	gpuErrchk(cudaDeviceSynchronize());
	squareMetropolisStepGPU<<<grid, thread_block, shared_size * sizeof(int)>>>(1, dead_spin, length);
	gpuErrchk(cudaDeviceSynchronize());
}

__host__ void cudaMeasureEnergyMagnetization() {
	int n = length * length;
	int threads = (n < MAX_THREADS) ? nextPow2((n + 1) / 2) + 1 : MAX_THREADS;
	int blocks = (n + (threads * 2 - 1)) / (threads * 2);
	//int smem_size = (threads <= 32) ? 4 * threads * sizeof(double) : 2 * threads * sizeof(double);
	int smem_size = 2 * threads * sizeof(double);
	
	sumDeltaReduction<<<blocks,threads,smem_size>>>(n);
	
	double h_delta_energy[blocks], h_delta_magnetization[blocks];
	gpuErrchk(cudaMemcpyFromSymbol(h_delta_energy, d_delta_energy, blocks * sizeof(double)));
	gpuErrchk(cudaMemcpyFromSymbol(h_delta_magnetization, d_delta_magnetization, blocks * sizeof(double)));
	for(int i = 1; i < blocks; i++) {
		//printf("adding h_d_m[%d] = %f to h_d_m[0] = %f\n", i, h_delta_magnetization[i], h_delta_magnetization[0]);
		h_delta_energy[0] += h_delta_energy[i];
		h_delta_magnetization[0] += h_delta_magnetization[i];
	}
	
	double* temp_delta_energy;
    gpuErrchk(cudaGetSymbolAddress((void **)&temp_delta_energy, d_delta_energy));
    gpuErrchk(cudaMemset(temp_delta_energy, 0, sizeof(spin) * sizeof(double) / sizeof(int)));
    double* temp_delta_magnetization;
    gpuErrchk(cudaGetSymbolAddress((void **)&temp_delta_magnetization, d_delta_magnetization));
    gpuErrchk(cudaMemset(temp_delta_magnetization, 0, sizeof(spin) * sizeof(double) / sizeof(int)));
	
	double normalization = 1. / (length * length);
	energy += h_delta_energy[0] * normalization;
	magnetization += h_delta_magnetization[0] * normalization;
}
