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

#include "integrate_LIF.h"

int main() {
	double params[6] = {
		500e-9, // Cm [ms/ohm]
		25e-9, // g [1/ohm]
		-70, // Vrest [mV]
		-52, // Vth [mV]
		-59,  // Vreset [mV]
		1. // taurefr [ms]
	};
	double h = 1e-3; // [ms]
	int N = 100000;
	double deltaI = 13.75e-9*300;
	std::vector<double> I = int_lif::currents::square_wave(N+1, 10000, deltaI/2, -10, deltaI/2); // [mA]
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
		V1 = int_lif::fwdEuler(V0, h, N, I, params);
		t1 = int_lif::utils::linspace(t0, t0 + (N+1)*h, N+1);
	}
	end1 = std::chrono::steady_clock::now();
	// for(int i = 0; i < N+1; i++) std::cout << I[i] << "->" << V1[i] << ";" << t1[i] << std::endl;
	std::cout << "fwdEuler: " << std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1).count() / 1000. << " ms" << std::endl;
	graph[0] = new TGraph(N+1, t1.data(), V1.data());
	graph[0]->SetMarkerColor(1);
	graph[0]->SetMarkerStyle(8);
	graph[0]->SetMarkerSize(0.5);
	multigraph->Add(graph[0]);
	legend->AddEntry(graph[0], "fwdEuler", "p");
	
	// bwdEuler
	std::chrono::time_point<std::chrono::steady_clock> start2, end2;
	start2 = std::chrono::steady_clock::now();
	std::vector<double> V2, t2;
	for(int jj = 0; jj < nrep; jj++) {
		V2 = int_lif::bwdEuler(V0, h, N, I, params);
		t2 = int_lif::utils::linspace(t0, t0 + (N+1)*h, N+1);
	}
	end2 = std::chrono::steady_clock::now();
	std::cout << "bwdEuler: " << std::chrono::duration_cast<std::chrono::microseconds>(end2 - start2).count() / 1000. << " ms" << std::endl;
	graph[1] = new TGraph(N+1, t2.data(), V2.data());
	graph[1]->SetMarkerColor(2);
	graph[1]->SetMarkerStyle(8);
	graph[1]->SetMarkerSize(0.5);
	multigraph->Add(graph[1]);
	legend->AddEntry(graph[1], "bwdEuler", "p");
	
	/* auto iGraph = new TGraph(N+1, t1.data(), I.data());
	iGraph->SetMarkerColor(7);
	multigraph->Add(iGraph);
	legend->AddEntry(iGraph, "I(t)", "p"); */
	
	canvas->cd();
	canvas->SetGrid();
	multigraph->Draw("AP");
	multigraph->SetTitle(";t [ms];V [mV]");
	legend->Draw();
	canvas->SaveAs("Euler_comparison.pdf");
	
}