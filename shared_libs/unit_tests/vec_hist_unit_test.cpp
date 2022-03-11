#include <iostream>
#include "vec_hist.h"

#define N 100

int main() {
	std::cout << "Testing vecHist: ";
	
	std::vector<float> input;
	input.reserve(N);
	for(unsigned i = 0; i < N; i++) input.push_back(i * i / 2.);
	
	auto output = vecHist::makeHist<float>(input, 10, 0, N * N / 2);
	
	bool result = true;
	std::vector<unsigned> pred_output {32, 13, 10, 9, 7, 7, 6, 6, 5, 5};
	for(unsigned i = 0; i < 10; i++) if(output[i] != pred_output[i]) result = false;
	
	std::string res = result ? "passed" : "failed";
	std::cout << res << std::endl << std::endl;
}