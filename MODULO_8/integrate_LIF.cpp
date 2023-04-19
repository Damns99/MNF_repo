#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <math.h>
#include <chrono>
#include <ctime>
#include <filesystem>
namespace fs = std::filesystem;

// single neuron Leaky Integrate and Fire model with external input
// {
//    Cm * dV/dt = -g * (V - Vrest) + I(t)
//    if V >= Vth => V -> V_reset
// }

std::vector<double> linspace(double x0, double x1, double n = 100) {
	std::vector<double> res(n);
	double dx = (x1 - x0) / n;
	for(int ii = 0; ii < n; ii++) res[ii] = x0 + dx * ii;
	return res;
}

std::vector<double> logspace(double e0, double e1, double n = 100) {
	std::vector<double> res(n);
	double de = (e1 - e0) / n;
	for(int ii = 0; ii < n; ii++) res[ii] = pow(10, e0 + de * ii);
	return res;
}

// Forward Euler method
// dV(n)/dt = (V(n+1) - V(n)) / h
std::vector<double> fwdEuler1(double V0, double h, int N, std::vector<double> I, double params[5]) {
	double Cm = params[0], g = params[1], Vrest = params[2], Vth = params[3], Vreset = params[4];
	std::vector<double> V(N+1);
	V[0] = V0;
	for(int n = 0; n < N; n++) {
		V[n+1] = V[n] * (1 - g/Cm*h) + h/Cm * I[n];
	}
	return V;
}

std::pair<std::vector<double>,std::vector<double>> fwdEuler2(double V0, double t0, double h, int N, std::vector<double> I, double params[5]) {
	double Cm = params[0], g = params[1], Vrest = params[2], Vth = params[3], Vreset = params[4];
	std::vector<double> V(N+1), t(N+1);
	V[0] = V0; t[0] = t0;
	for(int n = 0; n < N; n++) {
		V[n+1] = V[n] * (1 - g/Cm*h) + h/Cm * I[n];
		t[n+1] = t0 + (n+1) * h;
	}
	return std::make_pair(V,t);
}

int main() {
	std::chrono::time_point<std::chrono::steady_clock> start, end;
	start = std::chrono::steady_clock::now();
	end = std::chrono::steady_clock::now();
	std::cout << "updated + measured in: " << (end - start).count() / 1000000. << " s" << std::endl;
}