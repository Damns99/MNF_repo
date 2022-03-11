#include <iostream>
#include "text_io.h"

#define N 4

int main() {
	std::cout << "Testing textIo: ";
	
	int x[N] = {1, 2, 3, 4};
	std::string y[N] = {"a", "bb", "ccc", "dddd"};
	float z[N] = {1.23, 2.34, 3.45, 4.56};
	textIo::textOut("text_io.txt", '\t', '#', " text_io.txt", N, false, x, y, z);
	
	int a[N];
	std::string b[N];
	float c[N];
	textIo::textIn("text_io.txt", '\t', '#', a, b, c);
	
	bool result = true;
	if(result) for(unsigned int i = 0; i < N; i++) if(a[i] != x[i]) result = false;
	if(result) for(unsigned int i = 0; i < N; i++) if(b[i] != y[i]) result = false;
	if(result) for(unsigned int i = 0; i < N; i++) if(c[i] != z[i]) result = false;
	
	std::string res = result ? "passed" : "failed";
	std::cout << res << std::endl << std::endl;
}