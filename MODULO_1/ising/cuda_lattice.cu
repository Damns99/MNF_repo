#include "lattice.h"

#include <cuda.h>

__global__ void squareMetropolisStepGPU(Lattice2D* lattice, double* rr, int size, int bw, int dead_spin) {
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	if(index >= size) return;
	int x = 2 * index + ((2 * index / lattice->length + bw) % 2);
	if(x == dead_spin) return; // avoid loops!
	int newspin = - lattice->spin[x];
	double denergy = -2. * newspin * lattice->extrafield;
	for (int j = 0; j < lattice->links_per_spin; j++) {
		denergy += -2. * JACC * newspin * lattice->spin[lattice->links[lattice->links_per_spin * x + j]];
	}
	double r = exp(-lattice->beta * denergy);
	if (rr[index] < r) lattice->spin[x] = newspin;
}

void cudaSquareUpdateMetropolis(Lattice2D* lattice) { // accettanza?
	assert(lattice->length % 2 == 0);
	int nbw = lattice->length * lattice->length / 2;
	int sizedbl = nbw * sizeof(double);
	double rrb[nbw], rrw[nbw];
	for(int ii = 0; ii < nbw; ii++) {
		rrb[ii] = lattice->gen.randF();
		rrw[ii] = lattice->gen.randF();
	}
	int dead_spin = lattice->gen.randL(0, lattice->length * lattice->length);
	
	double* d_rrb;
	double* d_rrw;
	Lattice2D* d_lattice;
	cudaMalloc((void **)&d_rrb, sizedbl);
	cudaMalloc((void **)&d_rrw, sizedbl);
	cudaMalloc((void **)&d_lattice, sizeof(Lattice2D));
	
	cudaMemcpy(d_rrb, rrb, sizedbl, cudaMemcpyHostToDevice);
	cudaMemcpy(d_rrw, rrw, sizedbl, cudaMemcpyHostToDevice);
	cudaMemcpy(d_lattice, lattice, sizeof(Lattice2D), cudaMemcpyHostToDevice);
	
	int nthreads = nbw > 256 ? 256 : nbw;
	int nblocks = (nbw - 1) / nthreads + 1;
	squareMetropolisStepGPU<<<nblocks,nthreads>>>(d_lattice, d_rrb, nbw, 0, dead_spin);
	cudaDeviceSynchronize();
	squareMetropolisStepGPU<<<nblocks,nthreads>>>(d_lattice, d_rrw, nbw, 1, dead_spin);
	cudaDeviceSynchronize();
	
	cudaMemcpy(lattice, d_lattice, sizeof(Lattice2D), cudaMemcpyDeviceToHost);
	lattice->calculateEnergyMagnetizaton();
	
	cudaFree(d_rrb);
	cudaFree(d_rrw);
	cudaFree(d_lattice);
}