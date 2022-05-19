#include "box_muller.h"

boxMuller::BoxMuller::BoxMuller(long seed) {
	long tmp = seed < 0 ? seed : -seed;
	r.reset(tmp);
}

void boxMuller::BoxMuller::reset(long seed) {
	r.reset(seed);
}

void boxMuller::BoxMuller::rand(float mean, float stdev, float& r1, float& r2) {
	float x1 = r.randF(), x2 = r.randF();
	float phi = 2. * M_PI * x1, rho2 = -log(1. - x2);
	r1 = sqrt(2. * stdev * stdev * rho2) * sin(phi) + mean;
	r2 = sqrt(2. * stdev * stdev * rho2) * cos(phi) + mean;
}

std::vector<float> boxMuller::BoxMuller::randToVec(float mean, float stdev, int n) {
	std::vector<float> res;
	res.reserve(n);
	float r1, r2;
	if(n % 2 != 0) {
		rand(mean, stdev, r1, r2);
		res.push_back(r1);
		n--;
	}
	for(int i = 0; i < n / 2; i++) {
		rand(mean, stdev, r1, r2);
		res.push_back(r1);
		res.push_back(r2);
	}
	return res;
}

void boxMuller::BoxMuller::randToFile(float mean, float stdev, int n, const std::string filepath) {
	std::ofstream outfile;
	outfile.open(filepath, std::fstream::out);
	float r1, r2;
	if(n % 2 != 0) {
		rand(mean, stdev, r1, r2);
		outfile << r1 << std::endl;
		n--;
	}
	for(int i = 0; i < n / 2; i++) {
		rand(mean, stdev, r1, r2);
		outfile << r1 << std::endl << r2 << std::endl;
	}
}