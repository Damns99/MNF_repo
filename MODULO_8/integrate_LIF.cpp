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
std::vector<double> fwdEuler1(double V0, double h, int N, std::vector<double>& I, double params[5]) {
	double Cm = params[0], g = params[1], Vrest = params[2], Vth = params[3], Vreset = params[4];
	std::vector<double> V(N+1);
	V[0] = V0;
	for(int n = 0; n < N; n++) {
		V[n+1] = V[n] * (1 - g/Cm*h) + h/Cm * I[n];
	}
	return V;
}

std::pair<std::vector<double>,std::vector<double>> fwdEuler2(double V0, double t0, double h, int N, std::vector<double>& I, double params[5]) {
	double Cm = params[0], g = params[1], Vrest = params[2], Vth = params[3], Vreset = params[4];
	std::vector<double> V(N+1), t(N+1);
	V[0] = V0; t[0] = t0;
	for(int n = 0; n < N; n++) {
		V[n+1] = V[n] * (1 - g/Cm*h) + h/Cm * I[n];
		t[n+1] = t0 + (n+1) * h;
	}
	return std::make_pair(V,t);
}

std::vector<double> fwdEuler3(double V0, double h, int N, std::vector<double>& I, double params[5]) {
	double Cm = params[0], g = params[1], Vrest = params[2], Vth = params[3], Vreset = params[4];
	std::vector<double> V;
	V.push_back(V0);
	for(int n = 0; n < N; n++) {
		V.push_back(V[n] * (1 - g/Cm*h) + h/Cm * I[n]);
	}
	return V;
}

void fwdEuler4(double V0, double h, int N, double I[], double params[5], double* V) {
	double Cm = params[0], g = params[1], Vrest = params[2], Vth = params[3], Vreset = params[4];
	V[0] = V0;
	for(int n = 0; n < N; n++) {
		V[n+1] = V[n] * (1 - g/Cm*h) + h/Cm * I[n];
	}
}

void fwdEuler5(double V0, double t0, double h, int N, double I[], double params[5], double* V, double* t) {
	double Cm = params[0], g = params[1], Vrest = params[2], Vth = params[3], Vreset = params[4];
	V[0] = V0; t[0] = t0;
	for(int n = 0; n < N; n++) {
		V[n+1] = V[n] * (1 - g/Cm*h) + h/Cm * I[n];
		t[n+1] = t0 + (n+1) * h;
	}
}

int main() {
	double params[5] = {
		1., // Cm [ms/ohm]
		1., // g [1/ohm]
		-70., // Vrest [mV]
		40., // Vth [mV]
		-65.  // Vreset [mV]
	};
	double h = 0.1; // ms
	int N = 1000;
	std::vector<double> I = square_wave(); // da definireeeeeeeeeeeeeeee // mA
	double t0 = 0.; // ms
	double V0 = -70.; // mV
	
	{
	std::chrono::time_point<std::chrono::steady_clock> start, end;
	start = std::chrono::steady_clock::now();
	std::vector<double> V = fwdEuler1(V0, h, N, I, params);
	std::vector<double> t = linspace(t0, t0 + N*h, N);
	end = std::chrono::steady_clock::now();
	std::cout << "fwdEuler1: " << (end - start).count() / 1000000. << " s" << std::endl;
	}
	
	{
	std::chrono::time_point<std::chrono::steady_clock> start, end;
	start = std::chrono::steady_clock::now();
	std::pair<std::vector<double>,std::vector<double>> res = fwdEuler2(V0, t0, h, N, I, params);
	std::vector<double> V = res.first, t = res.second;
	end = std::chrono::steady_clock::now();
	std::cout << "fwdEuler2: " << (end - start).count() / 1000000. << " s" << std::endl;
	}
}