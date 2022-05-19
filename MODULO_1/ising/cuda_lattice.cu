#include "cuda_lattice.cuh"
#include <cuda.h>
#include <stdio.h>

__global__ void calculateEnergyMagnetizatonGPU(CudaLattice2D lattice) {
	*lattice.d_energy = 0.;
	*lattice.d_magnetization = 0.;
	for (int i = 0; i < lattice.length * lattice.length; i++) {
		for (int j = 0; j < lattice.links_per_spin; j++) {
			*lattice.d_energy += - JACC / 2. * lattice.d_spin[i] * lattice.d_spin[lattice.d_links[lattice.links_per_spin * i + j]] / (1. * lattice.length * lattice.length);
		}
		*lattice.d_energy += - lattice.d_spin[i] * lattice.extrafield / (1. * lattice.length * lattice.length);
		*lattice.d_magnetization += lattice.d_spin[i] / (1. * lattice.length * lattice.length);
	}
}

__global__ void squareMetropolisStepGPU(CudaLattice2D lattice, int bw, int dead_spin) {
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	if(index >= lattice.nbw) return;
	int x = 2 * index + ((2 * index / lattice.length + bw) % 2);
	if(x == dead_spin) return; // avoid loops!
	int newspin = - lattice.d_spin[x];
	double denergy = -2. * newspin * lattice.extrafield;
	for (int j = 0; j < lattice.links_per_spin; j++) {
		denergy += -2. * JACC * newspin * lattice.d_spin[lattice.d_links[lattice.links_per_spin * x + j]];
	}
	double r = exp(-lattice.beta * denergy);
	if (lattice.d_rr[2 * index + bw] < r) lattice.d_spin[x] = newspin;
}

__host__ CudaLattice2D::CudaLattice2D() {
	;
}

__host__ CudaLattice2D::~CudaLattice2D() {
	;
}

__host__ void CudaLattice2D::cudaInitFromLattice(Lattice2D* h_lattice) {
	assert(h_lattice->length % 2 == 0 && h_lattice->links_per_spin == 4);
	length = h_lattice->length;
	links_per_spin = h_lattice->links_per_spin;
	beta = h_lattice->beta;
	extrafield = h_lattice->extrafield;
	nbw = length * length / 2;
	nthreads = nbw > 256 ? 256 : nbw;
	nblocks = (nbw - 1) / nthreads + 1;
	gen = h_lattice->gen;
	
	int device_number;
	gpuErrchk(cudaGetDevice(&device_number));
	std::cout << "device:" << device_number << std::endl;
	
	size_t totsize;
	totsize = length * length * (1 + links_per_spin) * sizeof(int) + (2 * nbw + 2) * sizeof(double);
	gpuErrchk(cudaDeviceSetLimit(cudaLimitStackSize, totsize));
	
	gpuErrchk(cudaMalloc((void **)&d_spin, length * length * sizeof(int)));
	gpuErrchk(cudaMalloc((void **)&d_links, length * length * links_per_spin * sizeof(int)));
	gpuErrchk(cudaMalloc((void **)&d_rr, 2 * nbw * sizeof(double)));
	gpuErrchk(cudaMalloc((void **)&d_energy, sizeof(double)));
	gpuErrchk(cudaMalloc((void **)&d_magnetization, sizeof(double)));
	gpuErrchk(cudaMemcpy(d_spin, &h_lattice->spin, length * length * sizeof(int), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_links, &h_lattice->links, length * length * links_per_spin * sizeof(int), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_energy, &h_lattice->energy, sizeof(double), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_magnetization, &h_lattice->magnetization, sizeof(double), cudaMemcpyHostToDevice));
}

__host__ void CudaLattice2D::cudaRetrieveLattice(Lattice2D* h_lattice) {
	h_lattice->length = length;
	h_lattice->links_per_spin = links_per_spin;
	h_lattice->beta = beta;
	h_lattice->extrafield = extrafield;
	h_lattice->gen = gen;
	
	gpuErrchk(cudaMemcpy(&h_lattice->spin, d_spin, length * length * sizeof(int), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(&h_lattice->links, d_links, length * length * links_per_spin * sizeof(int), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(&h_lattice->energy, d_energy, sizeof(double), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(&h_lattice->magnetization, d_magnetization, sizeof(double), cudaMemcpyDeviceToHost));
}

__host__ void CudaLattice2D::cudaDestroyLattice() {
	gpuErrchk(cudaFree(d_spin));
	gpuErrchk(cudaFree(d_links));
	gpuErrchk(cudaFree(d_rr));
	gpuErrchk(cudaFree(d_energy));
	gpuErrchk(cudaFree(d_magnetization));
	cudaDeviceReset();
}

__host__ void CudaLattice2D::cudaMeasureEnergyMagnetization(double* h_energy, double* h_magnetization) {
	calculateEnergyMagnetizatonGPU<<<1,1>>>(*this);
	gpuErrchk(cudaDeviceSynchronize());
	gpuErrchk(cudaMemcpy(h_energy, d_energy, sizeof(double), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(h_magnetization, d_magnetization, sizeof(double), cudaMemcpyDeviceToHost));
}

__host__ void CudaLattice2D::cudaUpdateMetropolis() {
	double rr[2 * nbw];
	for(int ii = 0; ii < 2 * nbw; ii++) rr[ii] = gen.randF();
	gpuErrchk(cudaMemcpy(d_rr, rr, 2 * nbw * sizeof(double), cudaMemcpyHostToDevice));
	int dead_spin = gen.randL(0, 2 * nbw);
	
	squareMetropolisStepGPU<<<nblocks,nthreads>>>(*this, 0, dead_spin);
	gpuErrchk(cudaDeviceSynchronize());
	squareMetropolisStepGPU<<<nblocks,nthreads>>>(*this, 1, dead_spin);
	gpuErrchk(cudaDeviceSynchronize());
}