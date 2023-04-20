#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <math.h>
#include <chrono>
#include <ctime>
#include <filesystem>
namespace fs = std::filesystem;

#include <TGraph.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TFrame.h>
#include <TLegend.h>

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

std::vector<double> square_wave(int n, int period, double amplitude, int phase, double offset) {
	std::vector<double> res(n);
	for(int ii = 0; ii < n; ii++) {
		if((ii + phase + period) % period < period/2) res[ii] = offset + amplitude;
		else res[ii] = offset - amplitude;
	}
	return res;
}

// Forward Euler method
// dV(n)/dt = (V(n+1) - V(n)) / h
std::vector<double> fwdEuler(double V0, double h, int N, std::vector<double>& I, double params[5]) {
	double Cm = params[0], g = params[1], Vrest = params[2], Vth = params[3], Vreset = params[4];
	std::vector<double> V(N+1);
	V[0] = V0;
	for(int n = 0; n < N; n++) {
		V[n+1] = V[n] * (1 - g/Cm*h) + h/Cm * (I[n] + g * Vrest);
		if(V[n+1] >= Vth) V[n+1] = Vreset;
	}
	return V;
}

// Backward Euler method
// dV(n)/dt = (V(n+1) - V(n)) / h
std::vector<double> bwdEuler(double V0, double h, int N, std::vector<double>& I, double params[5]) {
	double Cm = params[0], g = params[1], Vrest = params[2], Vth = params[3], Vreset = params[4];
	std::vector<double> V(N+1);
	V[0] = V0;
	for(int n = 1; n <= N; n++) {
		V[n] = (V[n-1] + h/Cm * (I[n] + g * Vrest)) / (1 + g/Cm*h);
		if(V[n] >= Vth) V[n] = Vreset;
	}
	return V;
}

int main() {
	double params[5] = {
		500e-9, // Cm [ms/ohm]
		25e-9, // g [1/ohm]
		-70, // Vrest [mV]
		-52, // Vth [mV]
		-59  // Vreset [mV]
	};
	double h = 1e-1; // [ms]
	int N = 100;
	double deltaI = 13.75e-9*300;
	std::vector<double> I = square_wave(N, 100, deltaI/2, -10, deltaI/2); // [mA]
	double t0 = 0; // [ms]
	double V0 = -70; // [mV]
	
	int nrep = 100;
	TCanvas* canvas = new TCanvas("canvas0", "Canvas0", 600, 400);
	TMultiGraph* multigraph = new TMultiGraph();
	TLegend* legend = new TLegend(0.75, 0.75, 0.85, 0.85);
	TGraph* graph[5];
	
	// fwdEuler
	std::chrono::time_point<std::chrono::steady_clock> start1, end1;
	start1 = std::chrono::steady_clock::now();
	std::vector<double> V1, t1;
	for(int jj = 0; jj < nrep; jj++) {
		V1 = fwdEuler(V0, h, N, I, params);
		t1 = linspace(t0, t0 + (N+1)*h, N+1);
	}
	end1 = std::chrono::steady_clock::now();
	// for(int i = 0; i < N+1; i++) std::cout << V1[i] << ";" << t1[i] << std::endl;
	std::cout << "fwdEuler: " << std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1).count() / 1000. << " ms" << std::endl;
	graph[0] = new TGraph(N+1, t1.data(), V1.data());
	graph[0]->SetMarkerColor(1);
	multigraph->Add(graph[0]);
	legend->AddEntry(graph[0], "fwdEuler", "p");
	
	// bwdEuler
	std::chrono::time_point<std::chrono::steady_clock> start2, end2;
	start2 = std::chrono::steady_clock::now();
	std::vector<double> V2, t2;
	for(int jj = 0; jj < nrep; jj++) {
		V2 = bwdEuler(V0, h, N, I, params);
		t2 = linspace(t0, t0 + (N+1)*h, N+1);
	}
	end2 = std::chrono::steady_clock::now();
	std::cout << "bwdEuler: " << std::chrono::duration_cast<std::chrono::microseconds>(end2 - start2).count() / 1000. << " ms" << std::endl;
	graph[1] = new TGraph(N+1, t2.data(), V2.data());
	graph[1]->SetMarkerColor(2);
	multigraph->Add(graph[1]);
	legend->AddEntry(graph[1], "bwdEuler", "p");
	
	/* auto iGraph = new TGraph(N+1, t1.data(), I.data());
	iGraph->SetMarkerColor(7);
	multigraph->Add(iGraph);
	legend->AddEntry(iGraph, "I(t)", "p"); */
	
	canvas->cd();
	canvas->SetGrid();
	multigraph->Draw("AP");
	legend->Draw();
	canvas->SaveAs("Euler_comparison.pdf");
	
}