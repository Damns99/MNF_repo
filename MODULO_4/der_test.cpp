#include <iostream>
#include <string>

#include "bound_cond_vecs.h"
#include "derivators.h"

int main() {
	std::cout << "Testing derivators:" << std::endl;
	
	int N = 10;
	bound_cond_vecs::BoundCondVec<double> v1(N);
	for(int ii = 0; ii < N; ii++) v1[ii] = ii;
	
	bound_cond_vecs::BoundCondVec<double> dv1 = derivators::fwd_derive(v1, 1.);
	for(int ii = 0; ii < N; ii++) std::cout << v1[ii] << " ";
	std::cout << std::endl;
	for(int ii = 0; ii < N; ii++) std::cout << dv1[ii] << " ";
	std::cout << std::endl;
	
	bool result = true;
	for(int ii = 0; ii < N - 1; ii++) if(dv1[ii] != 1.) result = false;
	
	std::string res = result ? "passed" : "failed";
	std::cout << res << std::endl << std::endl;
}