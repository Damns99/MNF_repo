#include "cuda_lattice.cuh"
#include <cuda.h>
#include <stdio.h>

__host__ void cudaInitFromLattice() {
	grid = dim3(iDivUp(length, 2 * BLOCK_SIDE));
	
	gpuErrchk(cudaMemcpyToSymbol(d_length, &length, sizeof(int)));
	gpuErrchk(cudaMemcpyToSymbol(d_beta, &beta, sizeof(double)));
	gpuErrchk(cudaMemcpyToSymbol(d_p_length, &p_length, sizeof(int)));
	gpuErrchk(cudaMemcpyToSymbol(d_nparticles, &nparticles, sizeof(int)));
	gpuErrchk(cudaMemcpyToSymbol(d_y, y, sizeof(y)));
	gpuErrchk(cudaMemcpyToSymbol(d_links, links, sizeof(links)));
	
	/* gpuErrchk(cudaMemcpyToSymbol(d_rules, rules, sizeof(rules)));
	gpuErrchk(cudaMemcpyToSymbol(d_repetitions, repetitions, sizeof(repetitions)));
	gpuErrchk(cudaMemcpyToSymbol(d_nrules, &nrules, sizeof(int))); */
	
	gpuErrchk(cudaMemcpyToSymbol(d_obs1, obs1, sizeof(obs1)));
	gpuErrchk(cudaMemcpyToSymbol(d_obs2, obs2, sizeof(obs2)));
	
	gpuErrchk(cudaMemcpyToSymbol(d_pars, pars, sizeof(pars)));
    
	double* temp_delta_obs1;
	gpuErrchk(cudaGetSymbolAddress((void **)&temp_delta_obs1, d_delta_obs1));
	gpuErrchk(cudaMemset(temp_delta_obs1, 0, sizeof(y)));
	double* temp_delta_obs2;
	gpuErrchk(cudaGetSymbolAddress((void **)&temp_delta_obs2, d_delta_obs2));
	gpuErrchk(cudaMemset(temp_delta_obs2, 0, sizeof(y)));
	
	func_ds = harm_pot;
}

__host__ void cudaRetrieveLattice() {
	gpuErrchk(cudaMemcpyFromSymbol(y, d_y, sizeof(y)));
}

__host__ void cudaDestroyLattice() {
	cudaDeviceReset();
}

// 2) delta + GPU
__global__ void metropolisStepGPU(int bw, int dead_site, int length_, Function_ds func) {	
	// thread position
	int j = blockIdx.x * blockDim.x + threadIdx.x;
	// spin position
	int jj = 2 * j + (bw % 2);
	if(jj >= length_) return;
	
	if(cuda_index1d(jj,length_) == dead_spin) return; // avoid loops!
	
	double eta = d_beta / d_p_length;
	double delta = 2. * sqrt(eta);
	
	double y0 = d_y[cuda_index1d(jj,length_)]; //centre
	double y1 = d_y[d_links[2*jj+0]]; //left
	double y2 = d_y[d_links[2*jj+1]]; //right
	double yp = (d_rr[cuda_index1d(jj,length_) * 2. - 1.) * delta + y0;
	
	double outputs[2];
	double ds = func(y0, y1, y2, yp, d_pars, outputs);
	double r = exp(-ds);
	if (d_rr[cuda_index1d(jj,length_)+length_] < r) {
        d_y[cuda_index1d(jj,length_)] = yp;
		d_delta_obs1[cuda_index1d(jj,length_)] += outputs[0];
		d_delta_obs2[cuda_index1d(jj,length_)] += outputs[1];
    }
}

__host__ void cudaUpdateMetropolis() {
	double rr[2*length];
	for(int ii = 0; ii < 2*length; ii++) rr[ii] = gen.randF();
	gpuErrchk(cudaMemcpyToSymbol(d_rr, rr, sizeof(rr)));
	int dead_site = gen.randL(0, length);
    
	metropolisStepGPU<<<grid, thread_block>>>(0, dead_site, length, func_ds);
	gpuErrchk(cudaDeviceSynchronize());
	metropolisStepGPU<<<grid, thread_block>>>(1, dead_site, length, func_ds);
	gpuErrchk(cudaDeviceSynchronize());
}

__host__ void calculateObsGPU() {	
	double h_delta_obs1[length], h_delta_obs2[length];
	gpuErrchk(cudaMemcpyFromSymbol(h_delta_obs1, d_delta_obs1, length * sizeof(double)));
	gpuErrchk(cudaMemcpyFromSymbol(h_delta_obs2, d_delta_obs2, length * sizeof(double)));
	
	double normalization = 1. / p_length;
	for(int i = 0; i < nparticles; i++) {
		double sum1 = 0., sum2 = 0.;
		for(int j = 0; j < p_length; j++) {
			sum1 += h_delta_obs1[j + p_length * i];
			sum2 += h_delta_obs2[j + p_length * i];
		}
		obs1[i] += sum1 * normalization;
		obs2[i] += sum2 * normalization;
	}
	
	double* temp_delta_energy;
    gpuErrchk(cudaGetSymbolAddress((void **)&temp_delta_energy, d_delta_energy));
    gpuErrchk(cudaMemset(temp_delta_energy, 0, sizeof(y)));
    double* temp_delta_magnetization;
    gpuErrchk(cudaGetSymbolAddress((void **)&temp_delta_magnetization, d_delta_magnetization));
    gpuErrchk(cudaMemset(temp_delta_magnetization, 0, sizeof(y)));
}