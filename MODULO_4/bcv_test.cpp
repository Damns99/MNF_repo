#include <iostream>
#include <string>

#include "bound_cond_vecs.h"

int main() {
	std::cout << "Testing bound_cond_vecs:" << std::endl;
	
	int N = 7;
	double v[N] = {};
	for(int i = 0; i < N; i++) {
		v[i] = 3.14 * ((13 * i) % N);
		std::cout << v[i] << '\t';
	}
	std::cout << std::endl;
	
	bound_cond_vecs::BoundCondVec v1(N, v, PERIODIC_BC);
	bound_cond_vecs::BoundCondVec v2(N, v, ABSORBANT_BC);
	bound_cond_vecs::BoundCondVec v3(N, v, REFLECTIVE_BC);
	
	for(int i = -2 * N; i < 3 * N; i++) std::cout << i << '\t' << v1[i] << '\t' << v2[i] << '\t' << v3[i] << std::endl;
	
	bool result = true;
	if(v1[-3] != v1[2 * N - 3]) result = false; 
	else if(int(v2[-3]) != 0) result = false; 
	else if(v3[-3] != v3[2]) result = false; 
	
	std::string res = result ? "passed" : "failed";
	std::cout << res << std::endl << std::endl;
}