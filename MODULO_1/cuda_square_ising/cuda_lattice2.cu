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
    
    int* temp_delta_energy;
    gpuErrchk(cudaGetSymbolAddress((void **)&temp_delta_energy, d_delta_energy));
    gpuErrchk(cudaMemset(temp_delta_energy, 0, sizeof(spin)));
    int* temp_delta_magnetization;
    gpuErrchk(cudaGetSymbolAddress((void **)&temp_delta_magnetization, d_delta_magnetization));
    gpuErrchk(cudaMemset(temp_delta_magnetization, 0, sizeof(spin)));
}

__host__ void cudaRetrieveLattice() {
	gpuErrchk(cudaMemcpyFromSymbol(spin, d_spin, sizeof(spin)));
	gpuErrchk(cudaMemcpyFromSymbol(&energy, d_energy, sizeof(double)));
	gpuErrchk(cudaMemcpyFromSymbol(&magnetization, d_magnetization, sizeof(double)));
}

__host__ void cudaDestroyLattice() {
	cudaDeviceReset();
}

// .......................... temporary .......................
//__host__ void cudaMeasureEnergyMagnetization() {
	// .................................. calculateEnergyMagnetizatonGPU<<<1, 1>>>();
//	gpuErrchk(cudaDeviceSynchronize());
//	gpuErrchk(cudaMemcpyFromSymbol(&energy, d_energy, sizeof(double)));
//	gpuErrchk(cudaMemcpyFromSymbol(&magnetization, d_magnetization, sizeof(double)));
//}

//__host__ void cudaUpdateMetropolis() {
//	double rr[length * length];
//	for(int ii = 0; ii < length * length; ii++) rr[ii] = gen.randF();
//	gpuErrchk(cudaMemcpyToSymbol(d_rr, rr, sizeof(rr)));
//	int dead_spin = gen.randL(0, length * length);
	// ....................................int shared_size = (length + 2) * ((2 * nthreads - 1) / length + 3);
    /*
	squareMetropolisStepGPU<<<nblocks, nthreads, shared_size * sizeof(int)>>>(0, dead_spin);
	gpuErrchk(cudaDeviceSynchronize());
	squareMetropolisStepGPU<<<nblocks, nthreads, shared_size * sizeof(int)>>>(1, dead_spin);
	gpuErrchk(cudaDeviceSynchronize());
	*/
//}
// .......................... end temporary .....................

// 1) atomicAdd
__global__ void squareMetropolisStepGPU1(int bw, int dead_spin, int length_) {
    extern __shared__ int s_sub_spin[];
	
	// thread position
	int i = blockIdx.y * blockDim.y + threadIdx.y;
	int j = blockIdx.x * blockDim.x + threadIdx.x;
	// spin position
	int ii = i;
	int jj = 2 * j + ((i + bw) % 2);
	if((ii >= length_) || (jj >= length_)) return;
	// shared position
	int x = jj + 1;
	int y = ii + 1;
	
	s_sub_spin[index2d(y,x,BLOCK_SIDE+2)] = d_spin[index2d(ii,jj,length_)]; //centre
	s_sub_spin[index2d(y,x-1,BLOCK_SIDE+2)] = d_spin[index2d(ii,jj-1,length_)]; //left
    if(threadIdx.y == 0) s_sub_spin[index2d(y+1,x,BLOCK_SIDE+2)] = d_spin[index2d(ii+1,jj,length_)]; //up
    if((threadIdx.y == blockDim.y - 1) || (ii == length_ - 1)) s_sub_spin[index2d(y-1,x,BLOCK_SIDE+2)] = d_spin[index2d(ii-1,jj,length_)]; //down
    if((threadIdx.x == blockDim.x - 1) || (jj == length_ - 1)) s_sub_spin[index2d(y,x+1,BLOCK_SIDE+2)] = d_spin[index2d(ii,jj+1,length_)]; //right
	
    __syncthreads();
    
	if(index2d(ii,jj,length_) == dead_spin) return; // avoid loops!
	int newspin = - s_sub_spin[index2d(y,x,BLOCK_SIDE+2)];
	int neighbours = s_sub_spin[index2d(y-1,x,BLOCK_SIDE+2)] + s_sub_spin[index2d(y+1,x,BLOCK_SIDE+2)] + s_sub_spin[index2d(y,x-1,BLOCK_SIDE+2)] + s_sub_spin[index2d(y,x+1,BLOCK_SIDE+2)];
	
	double denergy = -2. * newspin * d_extrafield;
	denergy += JACC * -2. * newspin * neighbours;
	
	double r = exp(- d_beta * denergy);
	if (d_rr[index2d(ii,jj,length_)] < r) {
        d_spin[index2d(ii,jj,length_)] = newspin;
		myAtomicAdd(&d_energy, denergy);
		myAtomicAdd(&d_magnetization, -2. * newspin);
    }
}

__host__ void cudaUpdateMetropolis1() {
	double rr[length * length];
	for(int ii = 0; ii < length * length; ii++) rr[ii] = gen.randF();
	gpuErrchk(cudaMemcpyToSymbol(d_rr, rr, sizeof(rr)));
	int dead_spin = gen.randL(0, length * length);
    
	squareMetropolisStepGPU1<<<grid, thread_block, shared_size * sizeof(int)>>>(0, dead_spin, length);
	gpuErrchk(cudaDeviceSynchronize());
	squareMetropolisStepGPU1<<<grid, thread_block, shared_size * sizeof(int)>>>(1, dead_spin, length);
	gpuErrchk(cudaDeviceSynchronize());
}

__host__ void cudaMeasureEnergyMagnetization1() {
	gpuErrchk(cudaMemcpyFromSymbol(&energy, d_energy, sizeof(double)));
	gpuErrchk(cudaMemcpyFromSymbol(&magnetization, d_magnetization, sizeof(double)));
	
	double normalization = 1. / (length * length);
	energy *= normalization;
	magnetization *= normalization;
}

// 2) delta + GPU
__global__ void squareMetropolisStepGPU2(int bw, int dead_spin, int length_) {
    extern __shared__ int s_sub_spin[];
	
	// thread position
	int i = blockIdx.y * blockDim.y + threadIdx.y;
	int j = blockIdx.x * blockDim.x + threadIdx.x;
	// spin position
	int ii = i;
	int jj = 2 * j + ((i + bw) % 2);
	if((ii >= length_) || (jj >= length_)) return;
	// shared position
	int x = jj + 1;
	int y = ii + 1;
	
	s_sub_spin[index2d(y,x,BLOCK_SIDE+2)] = d_spin[index2d(ii,jj,length_)]; //centre
	s_sub_spin[index2d(y,x-1,BLOCK_SIDE+2)] = d_spin[index2d(ii,jj-1,length_)]; //left
    if(threadIdx.y == 0) s_sub_spin[index2d(y-1,x,BLOCK_SIDE+2)] = d_spin[index2d(ii-1,jj,length_)]; //up
    if((threadIdx.y == blockDim.y - 1) || (ii == length_ - 1)) s_sub_spin[index2d(y+1,x,BLOCK_SIDE+2)] = d_spin[index2d(ii+1,jj,length_)]; //down
    if((threadIdx.x == blockDim.x - 1) || (jj == length_ - 1)) s_sub_spin[index2d(y,x+1,BLOCK_SIDE+2)] = d_spin[index2d(ii,jj+1,length_)]; //right
	
    __syncthreads();
    
	if(index2d(ii,jj,length_) == dead_spin) return; // avoid loops!
	int newspin = - s_sub_spin[index2d(y,x,BLOCK_SIDE+2)];
	int neighbours = s_sub_spin[index2d(y-1,x,BLOCK_SIDE+2)] + s_sub_spin[index2d(y+1,x,BLOCK_SIDE+2)] + s_sub_spin[index2d(y,x-1,BLOCK_SIDE+2)] + s_sub_spin[index2d(y,x+1,BLOCK_SIDE+2)];
	
	int denergy = -2. * newspin; // double ... * d_extrafield
	denergy += JACC * -2. * newspin * neighbours;
	
	double r = exp(- d_beta * denergy);
	if (d_rr[index2d(ii,jj,length_)] < r) {
        d_spin[index2d(ii,jj,length_)] = newspin;
		d_delta_energy[index2d(ii,jj,length_)] += denergy;
		d_delta_magnetization[index2d(ii,jj,length_)] += -2. * newspin;
    }
}

__global__ void sumDeltaReduction2(int n) {
	extern __shared__ int data[];
	
	int tid = threadIdx.x;
	int index = threadIdx.x + (blockIdx.x * 2) * blockDim.x;
	
	int esum = (index < n) ? d_delta_energy[index] : 0;
	int msum = (index < n) ? d_delta_magnetization[index] : 0;
	if(index + blockDim.x < n) {
		esum += d_delta_energy[index + blockDim.x];
		msum += d_delta_magnetization[index + blockDim.x];
	}
	data[tid] = esum;
	data[tid + blockDim.x] = msum;
	__syncthreads();
	
	for(int s = blockDim.x / 2; s > 0; s /= 2) {
		if(tid < s) {
			data[tid] += data[tid + s];
			data[blockDim.x + tid] += data[blockDim.x + tid + s];
		}
		__syncthreads();
	}
	
	if(tid == 0) {
		d_delta_energy[blockIdx.x] = data[0];
		d_delta_magnetization[blockIdx.x] = data[blockDim.x];
	}
}

__host__ void cudaUpdateMetropolis2() {
	double rr[length * length];
	for(int ii = 0; ii < length * length; ii++) rr[ii] = gen.randF();
	gpuErrchk(cudaMemcpyToSymbol(d_rr, rr, sizeof(rr)));
	int dead_spin = gen.randL(0, length * length);
    
	squareMetropolisStepGPU2<<<grid, thread_block, shared_size * sizeof(int)>>>(0, dead_spin, length);
	gpuErrchk(cudaDeviceSynchronize());
	squareMetropolisStepGPU2<<<grid, thread_block, shared_size * sizeof(int)>>>(1, dead_spin, length);
	gpuErrchk(cudaDeviceSynchronize());
}

__host__ void cudaMeasureEnergyMagnetization2() {
	int n = length;
	int threads = (n < MAX_THREADS) ? nextPow2((n + 1) / 2) : MAX_THREADS;
	int blocks = (n + (threads * 2 - 1)) / (threads * 2);
	int smem_size = (threads <= 32) ? 4 * threads * sizeof(int) : 2 * threads * sizeof(int);
	
	sumDeltaReduction2<<<blocks,threads,smem_size>>>(n);
	
	int h_delta_energy[blocks], h_delta_magnetization[blocks];
	gpuErrchk(cudaMemcpyFromSymbol(h_delta_energy, d_delta_energy, blocks * sizeof(int)));
	gpuErrchk(cudaMemcpyFromSymbol(h_delta_magnetization, d_delta_magnetization, blocks * sizeof(int)));
	for(int i = 1; i < blocks; i++) {
		h_delta_energy[0] += h_delta_energy[i];
		h_delta_magnetization[0] += h_delta_magnetization[i];
	}
	
	int* temp_delta_energy;
    gpuErrchk(cudaGetSymbolAddress((void **)&temp_delta_energy, d_delta_energy));
    gpuErrchk(cudaMemset(temp_delta_energy, 0, sizeof(spin)));
    int* temp_delta_magnetization;
    gpuErrchk(cudaGetSymbolAddress((void **)&temp_delta_magnetization, d_delta_magnetization));
    gpuErrchk(cudaMemset(temp_delta_magnetization, 0, sizeof(spin)));
	
	double normalization = 1. / (length * length);
	energy += h_delta_energy[0] * normalization;
	magnetization += h_delta_magnetization[0] * normalization;
}

// 3) delta + CPU
__host__ void sumDeltaSerial3() {
	int h_delta_energy[length * length], h_delta_magnetization[length * length];
	gpuErrchk(cudaMemcpyFromSymbol(h_delta_energy, d_delta_energy, length * length * sizeof(int)));
	gpuErrchk(cudaMemcpyFromSymbol(h_delta_magnetization, d_delta_magnetization, length * length * sizeof(int)));
	
	double normalization = 1. / (length * length);
	for (int i = 0; i < length * length; i++) {
		energy += h_delta_energy[i] * normalization;
		magnetization += h_delta_magnetization[i] * normalization;
	}
}

__host__ void cudaUpdateMetropolis3() {
	double rr[length * length];
	for(int ii = 0; ii < length * length; ii++) rr[ii] = gen.randF();
	gpuErrchk(cudaMemcpyToSymbol(d_rr, rr, sizeof(rr)));
	int dead_spin = gen.randL(0, length * length);
    
	squareMetropolisStepGPU2<<<grid, thread_block, shared_size * sizeof(int)>>>(0, dead_spin, length);
	gpuErrchk(cudaDeviceSynchronize());
	squareMetropolisStepGPU2<<<grid, thread_block, shared_size * sizeof(int)>>>(1, dead_spin, length);
	gpuErrchk(cudaDeviceSynchronize());
}

__host__ void cudaMeasureEnergyMagnetization3() {
	sumDeltaSerial3();
	
	int* temp_delta_energy;
    gpuErrchk(cudaGetSymbolAddress((void **)&temp_delta_energy, d_delta_energy));
    gpuErrchk(cudaMemset(temp_delta_energy, 0, sizeof(spin)));
    int* temp_delta_magnetization;
    gpuErrchk(cudaGetSymbolAddress((void **)&temp_delta_magnetization, d_delta_magnetization));
    gpuErrchk(cudaMemset(temp_delta_magnetization, 0, sizeof(spin)));
}

// 4) old

__global__ void calculateEnergyMagnetizatonGPU4(int length_) {
	d_energy = 0.;
	d_magnetization = 0.;
	for (int i = 0; i < length_; i++) for (int j = 0; j < length_; j++) {
        int curr_spin = d_spin[index2d(i,j,length_)];
		int neighbours = d_spin[index2d(i+1,j,length_)] + d_spin[index2d(i,j+1,length_)];
		d_energy += - JACC * curr_spin * neighbours;
		d_magnetization += curr_spin;
	}
	d_energy += - d_extrafield * d_magnetization;
}

__global__ void squareMetropolisStepGPU4(int bw, int dead_spin, int length_) {
    extern __shared__ int s_sub_spin[];
	
	// thread position
	int i = blockIdx.y * blockDim.y + threadIdx.y;
	int j = blockIdx.x * blockDim.x + threadIdx.x;
	// spin position
	int ii = i;
	int jj = 2 * j + ((i + bw) % 2);
	if((ii >= length_) || (jj >= length_)) return;
	// shared position
	int x = jj + 1;
	int y = ii + 1;
	
	s_sub_spin[index2d(y,x,BLOCK_SIDE+2)] = d_spin[index2d(ii,jj,length_)]; //centre
	s_sub_spin[index2d(y,x-1,BLOCK_SIDE+2)] = d_spin[index2d(ii,jj-1,length_)]; //left
    if(threadIdx.y == 0) s_sub_spin[index2d(y+1,x,BLOCK_SIDE+2)] = d_spin[index2d(ii+1,jj,length_)]; //up
    if((threadIdx.y == blockDim.y - 1) || (ii == length_ - 1)) s_sub_spin[index2d(y-1,x,BLOCK_SIDE+2)] = d_spin[index2d(ii-1,jj,length_)]; //down
    if((threadIdx.x == blockDim.x - 1) || (jj == length_ - 1)) s_sub_spin[index2d(y,x+1,BLOCK_SIDE+2)] = d_spin[index2d(ii,jj+1,length_)]; //right
	
    __syncthreads();
    
	if(index2d(ii,jj,length_) == dead_spin) return; // avoid loops!
	int newspin = - s_sub_spin[index2d(y,x,BLOCK_SIDE+2)];
	int neighbours = s_sub_spin[index2d(y-1,x,BLOCK_SIDE+2)] + s_sub_spin[index2d(y+1,x,BLOCK_SIDE+2)] + s_sub_spin[index2d(y,x-1,BLOCK_SIDE+2)] + s_sub_spin[index2d(y,x+1,BLOCK_SIDE+2)];
	
	double denergy = -2. * newspin * d_extrafield;
	denergy += JACC * -2. * newspin * neighbours;
	
	double r = exp(- d_beta * denergy);
	if (d_rr[index2d(ii,jj,length_)] < r) {
        d_spin[index2d(ii,jj,length_)] = newspin;
    }
}

__host__ void cudaUpdateMetropolis4() {
	double rr[length * length];
	for(int ii = 0; ii < length * length; ii++) rr[ii] = gen.randF();
	gpuErrchk(cudaMemcpyToSymbol(d_rr, rr, sizeof(rr)));
	int dead_spin = gen.randL(0, length * length);
    
	squareMetropolisStepGPU4<<<grid, thread_block, shared_size * sizeof(int)>>>(0, dead_spin, length);
	gpuErrchk(cudaDeviceSynchronize());
	squareMetropolisStepGPU4<<<grid, thread_block, shared_size * sizeof(int)>>>(1, dead_spin, length);
	gpuErrchk(cudaDeviceSynchronize());
}

__host__ void cudaMeasureEnergyMagnetization4() {
	calculateEnergyMagnetizatonGPU4<<<1,1>>>(length);
	
	gpuErrchk(cudaMemcpyFromSymbol(&energy, d_energy, sizeof(double)));
	gpuErrchk(cudaMemcpyFromSymbol(&magnetization, d_magnetization, sizeof(double)));
	
	double normalization = 1. / (length * length);
	energy *= normalization;
	magnetization *= normalization;
}

// 5) old ma parallello atomicAdd

//__global__ void calculateEnergyMagnetizatonGPU5() {
//	int index = threadIdx.x + blockIdx.x * blockDim.x;
//	int i = index / d_length, j = index % d_length; // sarebbe bello dividere in quadrati ed evitarsi i %
	
//	if (index == 0) {
//		d_energy = 0.;
//		d_magnetization = 0.;
//	}
//	__syncthreads();
	
//	int curr_spin = d_spin[index2d(i,j,d_length)];
//	int neighbours = d_spin[index2d(i+1,j,d_length)] + d_spin[index2d(i,j+1,d_length)];
//	myAtomicAdd(&d_energy, - curr_spin * neighbours);
//	myAtomicAdd(&d_magnetization, curr_spin);
	
	/* .................. da fare in measureEnergyMagnetization() ..................
	double normalization = 1. / (d_length * d_length);
	d_magnetization *= normalization;
	d_energy *= JACC * normalization;
	d_energy += - d_extrafield * magnetization;
	*/
//}

// 6) old ma parallelo reduction

//__global__ void calculateEnergyMagnetizatonGPU6() {	
//	int index = threadIdx.x + blockIdx.x * blockDim.x;
//	int i = index / d_length, j = index % d_length; // sarebbe bello dividere in quadrati ed evitarsi i %
	
//	if (index == 0) {
//		d_energy = 0.;
//		d_magnetization = 0.;
//	}
	
//	int curr_spin = d_spin[index2d(i,j,d_length)];
//	int neighbours = d_spin[index2d(i+1,j,d_length)] + d_spin[index2d(i,j+1,d_length)];
//	d_delta_energy[index2d(i,j,d_length)] = - curr_spin * neighbours;
//	d_delta_magnetization[index2d(i,j,d_length)] = curr_spin;
	
	/* .................. da fare in measureEnergyMagnetization() ..................
	double normalization = 1. / (d_length * d_length);
	d_magnetization *= normalization;
	d_energy *= JACC * normalization;
	d_energy += - d_extrafield * magnetization;
	*/
//}

