#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <filesystem>
namespace fs = std::filesystem;

#include "lattice.h"
#include "cuda_lattice.cuh"

#include <TError.h>

std::string outfoldername = "test";

constexpr int length_        = 32;
constexpr int init_mode_     = HOT;
constexpr double extrafield_ = 0.;
constexpr double beta_       = 0.3;

std::string makeFixedLength(const int i, const int l)
{
    std::stringstream s;
    s << std::setfill('0') << std::setw(l) << i;
    return s.str();
}

void calculateEnergyMagnetizaton_GPU() {
	energy = 0.;
	magnetization = 0.;
	
	for (int i = 0; i < length; i++) for (int j = 0; j < length; j++) {
		int current = spin[index2d(i,j,length)];
		int neighbours = spin[index2d(i+1,j,length)] + spin[index2d(i,j+1,length)];
		energy += - current * neighbours;
		magnetization += current;
	}
	
	double normalization = 1. / (length * length);
	magnetization *= normalization;
	energy *= JACC * normalization;
	energy += - extrafield * magnetization;
}

void timeTest(int nmeas, int ncycles, std::ofstream* measfile) {
	int lognmeas = log10(nmeas) + 1;
	
    std::chrono::time_point<std::chrono::steady_clock> start, end0, end1, end2, end3, end4, end5, end6;
	double time0 = 0., time1 = 0., time2 = 0., time3 = 0., time4 = 0., time5 = 0., time6 = 0.;
	
	for(int ii = 0; ii < nmeas; ii++) {		
        start = std::chrono::steady_clock::now();
		
		for(int jj = 0; jj < ncycles; jj++) cudaUpdateMetropolis1();
		cudaMeasureEnergyMagnetization1();
		
		end0 = std::chrono::steady_clock::now();
		time0 += 0.001 * (end0 - start).count() / nmeas / ncycles;
		
		*measfile << '1' << '\t' << energy << std::endl;
		
		start = std::chrono::steady_clock::now();
		
		for(int jj = 0; jj < ncycles; jj++) cudaUpdateMetropolis2();
		
		end1 = std::chrono::steady_clock::now();
		time1 += 0.001 * (end1 - start).count() / nmeas / ncycles;
		
		start = std::chrono::steady_clock::now();
		
		cudaMeasureEnergyMagnetization2();
		
		end2 = std::chrono::steady_clock::now();
		time2 += 0.001 * (end2 - start).count() / nmeas;
		
		*measfile << '2' << '\t' << energy << std::endl;
		
		start = std::chrono::steady_clock::now();
		
		for(int jj = 0; jj < ncycles; jj++) cudaUpdateMetropolis3();
		
		end3 = std::chrono::steady_clock::now();
		time3 += 0.001 * (end3 - start).count() / nmeas / ncycles;
		
		start = std::chrono::steady_clock::now();
		
		cudaMeasureEnergyMagnetization3();
		
		end4 = std::chrono::steady_clock::now();
		time4 += 0.001 * (end4 - start).count() / nmeas;
		
		*measfile << '3' << '\t' << energy << std::endl;
		
		start = std::chrono::steady_clock::now();
		
		for(int jj = 0; jj < ncycles; jj++) cudaUpdateMetropolis4();
		
		end5 = std::chrono::steady_clock::now();
		time5 += 0.001 * (end5 - start).count() / nmeas / ncycles;
		
		start = std::chrono::steady_clock::now();
		
		cudaMeasureEnergyMagnetization4();
		
		end6 = std::chrono::steady_clock::now();
		time6 += 0.001 * (end6 - start).count() / nmeas;
		
		*measfile << '4' << '\t' << energy << std::endl;
	}
	
	std::cout << std::endl;
	std::cout << "nmeas = " << nmeas << "\tncycles = " << ncycles << std::endl;
	std::cout << "1) updated + measured in: " << time0 << " us per cycle" << std::endl;
	std::cout << "2) updated + measured in: " << time1 + time2 / ncycles << " us per cycle" << std::endl;
	std::cout << "2) updated in:            " << time1 << " us per cycle" << std::endl;
	std::cout << "2) measured in:           " << time2 << " us" << std::endl;
	std::cout << "3) updated + measured in: " << time3 + time4 / ncycles << " us per cycle" << std::endl;
	std::cout << "3) updated in:            " << time3 << " us per cycle" << std::endl;
	std::cout << "3) measured in:           " << time4 << " us" << std::endl;
	std::cout << "4) updated + measured in: " << time5 + time6 / ncycles << " us per cycle" << std::endl;
	std::cout << "4) updated in:            " << time5 << " us per cycle" << std::endl;
	std::cout << "4) measured in:           " << time6 << " us" << std::endl;
	
	std::cout << std::endl;
	std::cout << "1) total time: " << 0.001 * time0 * nmeas * ncycles << " ms" << std::endl;
	std::cout << "2) total time: " << 0.001 * (time1 * nmeas * ncycles + time2 * nmeas) << " ms" << std::endl;
	std::cout << "3) total time: " << 0.001 * (time3 * nmeas * ncycles + time4 * nmeas) << " ms" << std::endl;
	std::cout << "4) total time: " << 0.001 * (time5 * nmeas * ncycles + time6 * nmeas) << " ms" << std::endl;
	std::cout << std::endl;
}

int main(int argc, char* argv[]) {
	gErrorIgnoreLevel = kWarning;
	
    long seed = -std::time(NULL);
	
	createLattice(length_, beta_, extrafield_, seed, init_mode_, "");
	cudaInitFromLattice();
	
	fs::current_path(fs::current_path() / "measures");
	fs::create_directory(outfoldername);
	fs::current_path(fs::current_path() / outfoldername);
	
	std::ofstream measfile;
	measfile.open("cuda_metro_ising_test_energy.txt", std::fstream::out);
	measfile << "#algorithm \tenergy" << std::endl;
	
	int ntotal = 1024;
	for(int nmeas = ntotal; nmeas >= 1; nmeas/=2) timeTest(nmeas, ntotal / nmeas, &measfile);
	
	measfile.close();
	
	cudaDestroyLattice();
	
	std::cout << std::endl;
	return 0;
}