#include <iostream>
#include <cuda.h>
#include <stdio.h>

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, std::string file, int line, bool abort=false)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file.c_str(), line);
      if (abort) exit(code);
   }
}

class myClass {
	public:
		int* d_value;
		size_t freem = 1, totm = 1;
		
		__host__ myClass() {;}
		__host__ ~myClass() {;}
		__host__ void initClass(int* init);
		__host__ void destroyClass();
		__host__ void runClass();
		__host__ void retrieveClass(int* h_value);
};

__global__ void setOne(myClass c) {
	*c.d_value = 1;
}

__host__ void myClass::initClass(int* init) {
	gpuErrchk(cudaMemGetInfo(&freem, &totm));
	std::cout << "free: " << freem << " total: " << totm << std::endl;
	gpuErrchk(cudaMalloc((void **)&d_value, sizeof(int)));
	gpuErrchk(cudaMemcpy(d_value, init, sizeof(int), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemGetInfo(&freem, &totm));
	std::cout << "free: " << freem << " total: " << totm << std::endl;
}
__host__ void myClass::destroyClass() {
	gpuErrchk(cudaMemGetInfo(&freem, &totm));
	std::cout << "free: " << freem << " total: " << totm << std::endl;
	gpuErrchk(cudaFree(d_value));
	gpuErrchk(cudaMemGetInfo(&freem, &totm));
	std::cout << "free: " << freem << " total: " << totm << std::endl;
}
__host__ void myClass::runClass() {
	setOne<<<1,1>>>(*this);
	gpuErrchk(cudaDeviceSynchronize());
}
__host__ void myClass::retrieveClass(int* h_value) {
	gpuErrchk(cudaMemGetInfo(&freem, &totm));
	std::cout << "free: " << freem << " total: " << totm << std::endl;
	gpuErrchk(cudaMemcpy(h_value, d_value, sizeof(int), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemGetInfo(&freem, &totm));
	std::cout << "free: " << freem << " total: " << totm << std::endl;
}

int main() {
	int a = 101;
	myClass c;
	c.initClass(&a);
	c.runClass();
	c.retrieveClass(&a);
	c.destroyClass();
}