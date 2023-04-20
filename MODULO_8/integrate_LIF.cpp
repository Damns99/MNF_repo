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
std::vector<double> fwdEuler1(double V0, double h, int N, std::vector<double>& I, double params[5]) {
	double Cm = params[0], g = params[1], Vrest = params[2], Vth = params[3], Vreset = params[4];
	std::vector<double> V(N+1);
	V[0] = V0;
	for(int n = 0; n < N; n++) {
		V[n+1] = V[n] * (1 - g/Cm*h) + h/Cm * (I[n] + g * Vrest);
		if(V[n+1] >= Vth) V[n+1] = Vreset;
	}
	return V;
}

std::pair<std::vector<double>,std::vector<double>> fwdEuler2(double V0, double t0, double h, int N, std::vector<double>& I, double params[5]) {
	double Cm = params[0], g = params[1], Vrest = params[2], Vth = params[3], Vreset = params[4];
	std::vector<double> V(N+1), t(N+1);
	V[0] = V0; t[0] = t0;
	for(int n = 0; n < N; n++) {
		V[n+1] = V[n] * (1 - g/Cm*h) + h/Cm * (I[n] + g * Vrest);
		t[n+1] = t0 + (n+1) * h;
		if(V[n+1] >= Vth) V[n+1] = Vreset;
	}
	return std::make_pair(V,t);
}

std::vector<double> fwdEuler3(double V0, double h, int N, std::vector<double>& I, double params[5]) {
	double Cm = params[0], g = params[1], Vrest = params[2], Vth = params[3], Vreset = params[4];
	std::vector<double> V;
	V.push_back(V0);
	for(int n = 0; n < N; n++) {
		V.push_back(V[n] * (1 - g/Cm*h) + h/Cm * (I[n] + g * Vrest));
		if(V[n+1] >= Vth) V[n+1] = Vreset;
	}
	return V;
}

void fwdEuler4(double V0, double h, int N, double* I, double params[5], double* V) {
	double Cm = params[0], g = params[1], Vrest = params[2], Vth = params[3], Vreset = params[4];
	V[0] = V0;
	for(int n = 0; n < N; n++) {
		V[n+1] = V[n] * (1 - g/Cm*h) + h/Cm * (I[n] + g * Vrest);
		if(V[n+1] >= Vth) V[n+1] = Vreset;
	}
}

void fwdEuler5(double V0, double t0, double h, int N, double* I, double params[5], double* V, double* t) {
	double Cm = params[0], g = params[1], Vrest = params[2], Vth = params[3], Vreset = params[4];
	V[0] = V0; t[0] = t0;
	for(int n = 0; n < N; n++) {
		V[n+1] = V[n] * (1 - g/Cm*h) + h/Cm * (I[n] + g * Vrest);
		t[n+1] = t0 + (n+1) * h;
		if(V[n+1] >= Vth) V[n+1] = Vreset;
	}
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
	
	// fwdEuler1
	std::chrono::time_point<std::chrono::steady_clock> start1, end1;
	start1 = std::chrono::steady_clock::now();
	std::vector<double> V1, t1;
	for(int jj = 0; jj < nrep; jj++) {
		V1 = fwdEuler1(V0, h, N, I, params);
		t1 = linspace(t0, t0 + (N+1)*h, N+1);
	}
	end1 = std::chrono::steady_clock::now();
	// for(int i = 0; i < N+1; i++) std::cout << V1[i] << ";" << t1[i] << std::endl;
	std::cout << "fwdEuler1: " << std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1).count() / 1000. << " ms" << std::endl;
	graph[0] = new TGraph(N+1, t1.data(), V1.data());
	graph[0]->SetMarkerColor(1);
	multigraph->Add(graph[0]);
	legend->AddEntry(graph[0], "fwdEuler1", "p");
	
	// fwdEuler2
	std::chrono::time_point<std::chrono::steady_clock> start2, end2;
	start2 = std::chrono::steady_clock::now();
	std::vector<double> V2, t2;
	std::pair<std::vector<double>,std::vector<double>> res2;
	for(int jj = 0; jj < nrep; jj++) {
		res2 = fwdEuler2(V0, t0, h, N, I, params);
		V2 = res2.first;
		t2 = res2.second;
	}
	end2 = std::chrono::steady_clock::now();
	std::cout << "fwdEuler2: " << std::chrono::duration_cast<std::chrono::microseconds>(end2 - start2).count() / 1000. << " ms" << std::endl;
	graph[1] = new TGraph(N+1, t2.data(), V2.data());
	graph[1]->SetMarkerColor(2);
	multigraph->Add(graph[1]);
	legend->AddEntry(graph[1], "fwdEuler2", "p");
	
	// fwdEuler3
	std::chrono::time_point<std::chrono::steady_clock> start3, end3;
	start3 = std::chrono::steady_clock::now();
	std::vector<double> V3, t3;
	for(int jj = 0; jj < nrep; jj++) {
		V3 = fwdEuler3(V0, h, N, I, params);
		t3 = linspace(t0, t0 + (N+1)*h, N+1);
	}
	end3 = std::chrono::steady_clock::now();
	std::cout << "fwdEuler3: " << std::chrono::duration_cast<std::chrono::microseconds>(end3 - start3).count() / 1000. << " ms" << std::endl;
	graph[2] = new TGraph(N+1, t3.data(), V3.data());
	graph[2]->SetMarkerColor(3);
	multigraph->Add(graph[2]);
	legend->AddEntry(graph[2], "fwdEuler3", "p");
	
	// fwdEuler4
	std::chrono::time_point<std::chrono::steady_clock> start4, end4;
	start4 = std::chrono::steady_clock::now();
	double V4[N+1], t4[N+1];
	for(int jj = 0; jj < nrep; jj++) {
		fwdEuler4(V0, h, N, I.data(), params, V4);
		for(int i = 0; i <= N; i++) t4[i] = t0 + i * h;
	}
	end4 = std::chrono::steady_clock::now();
	std::cout << "fwdEuler4: " << std::chrono::duration_cast<std::chrono::microseconds>(end4 - start4).count() / 1000. << " ms" << std::endl;
	graph[3] = new TGraph(N+1, t4, V4);
	graph[3]->SetMarkerColor(4);
	multigraph->Add(graph[3]);
	legend->AddEntry(graph[3], "fwdEuler4", "p");
	
	// fwdEuler5
	std::chrono::time_point<std::chrono::steady_clock> start5, end5;
	start5 = std::chrono::steady_clock::now();
	double V5[N+1], t5[N+1];
	for(int jj = 0; jj < nrep; jj++) {
		fwdEuler5(V0, t0, h, N, I.data(), params, V5, t5);
	}
	end5 = std::chrono::steady_clock::now();
	std::cout << "fwdEuler5: " << std::chrono::duration_cast<std::chrono::microseconds>(end5 - start5).count() / 1000. << " ms" << std::endl;
	graph[4] = new TGraph(N+1, t5, V5);
	graph[4]->SetMarkerColor(5);
	multigraph->Add(graph[4]);
	legend->AddEntry(graph[4], "fwdEuler5", "p");
	
	/* auto iGraph = new TGraph(N+1, t1.data(), I.data());
	iGraph->SetMarkerColor(7);
	multigraph->Add(iGraph);
	legend->AddEntry(iGraph, "I(t)", "p"); */
	
	canvas->cd();
	canvas->SetGrid();
	multigraph->Draw("AP");
	legend->Draw();
	canvas->SaveAs("fwdEuler_comparison.pdf");
	
}