#include <iostream>
#include <vector>
#include "vec_hist.h"
#include "ran2.h"

#define N 1000000
#define B 5
#define MEAN N/B

int main() {
	std::cout << "Testing ran2: ";
	
	long seed = -42;
	ran2::RandomGenerator gen(seed);
	std::vector<long> vec;
	vec.reserve(N);
	
	for(long ii = 0; ii < N; ii++) vec.push_back(gen.randL(10, 20));
	auto hist = vecHist::makeHist<long>(vec, B, 10, 20);
	
	bool result = true;
	for(int ii = 0; ii < B; ii++) {
		long aa = hist[ii] - MEAN;
		if(aa * aa > 9 * MEAN) {
			std::cout << aa << "\t" << hist[ii] << std::endl;
			result = false;
		}
	}
	
	gen.toFile("ran2.txt");
	ran2::RandomGenerator gen2(seed);
	gen2.fromFile("ran2.txt");
	if(gen.randF() != gen2.randF()) result = false;
	
	std::string res = result ? "passed" : "failed";
	std::cout << res << std::endl << std::endl;
}