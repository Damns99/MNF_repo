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
		1e-10, // Cm [ms/ohm]
		1e-8, // g [1/ohm]
		-70e-3, // Vrest [V]
		-50e-3, // Vth [V]
		-70e-3,  // Vreset [V]
		1e-3 // taurefr [s]
	};
	double h = 0.5e-3; // [s]
	int N = int(1./h);
	double aI = 0.3e-9; // [A]
	std::vector<double> I = int_lif::currents::pulse_train(N+2, int(0.1/h), int(0.01/h), aI, int(0.05/h), 25);
	double t0 = 0; // [s]
	double V0 = -70e-3; // [V]
	
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
		t1 = int_lif::utils::linspace(t0, t0 + N*h, N);
	}
	end1 = std::chrono::steady_clock::now();
	// for(int i = 0; i < N+1; i++) std::cout << I[i] << "->" << V1[i] << ";" << t1[i] << std::endl;
	std::cout << "fwdEuler: " << std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1).count() / 1000. << " ms" << std::endl;
	graph[0] = new TGraph(N, t1.data(), V1.data());
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
		t2 = int_lif::utils::linspace(t0, t0 + N*h, N);
	}
	end2 = std::chrono::steady_clock::now();
	std::cout << "bwdEuler: " << std::chrono::duration_cast<std::chrono::microseconds>(end2 - start2).count() / 1000. << " ms" << std::endl;
	graph[1] = new TGraph(N, t2.data(), V2.data());
	graph[1]->SetMarkerColor(2);
	graph[1]->SetMarkerStyle(8);
	graph[1]->SetMarkerSize(0.5);
	multigraph->Add(graph[1]);
	legend->AddEntry(graph[1], "bwdEuler", "p");
	
	// Heun
	std::chrono::time_point<std::chrono::steady_clock> start3, end3;
	start3 = std::chrono::steady_clock::now();
	std::vector<double> V3, t3;
	for(int jj = 0; jj < nrep; jj++) {
		V3 = int_lif::Heun(V0, h, N, I, params);
		t3 = int_lif::utils::linspace(t0, t0 + N*h, N);
	}
	end3 = std::chrono::steady_clock::now();
	std::cout << "Heun: " << std::chrono::duration_cast<std::chrono::microseconds>(end3 - start3).count() / 1000. << " ms" << std::endl;
	graph[2] = new TGraph(N, t3.data(), V3.data());
	graph[2]->SetMarkerColor(3);
	graph[2]->SetMarkerStyle(8);
	graph[2]->SetMarkerSize(0.5);
	multigraph->Add(graph[2]);
	legend->AddEntry(graph[2], "Heun", "p");
	
	// RK4
	std::vector<double> IRK4 = int_lif::currents::pulse_train(2*N+3, 2*int(0.1/h), 2*int(0.01/h), aI, 2*int(0.05/h), 25);
	std::chrono::time_point<std::chrono::steady_clock> start4, end4;
	start4 = std::chrono::steady_clock::now();
	std::vector<double> V4, t4;
	for(int jj = 0; jj < nrep; jj++) {
		V4 = int_lif::RK4(V0, h, N, IRK4, params);
		t4 = int_lif::utils::linspace(t0, t0 + N*h, N);
	}
	end4 = std::chrono::steady_clock::now();
	std::cout << "RK4: " << std::chrono::duration_cast<std::chrono::microseconds>(end4 - start4).count() / 1000. << " ms" << std::endl;
	graph[3] = new TGraph(N, t4.data(), V4.data());
	graph[3]->SetMarkerColor(4);
	graph[3]->SetMarkerStyle(8);
	graph[3]->SetMarkerSize(0.5);
	multigraph->Add(graph[3]);
	legend->AddEntry(graph[3], "RK4", "p");
	
	/* auto iGraph = new TGraph(N+1, t1.data(), I.data());
	iGraph->SetMarkerColor(7);
	multigraph->Add(iGraph);
	legend->AddEntry(iGraph, "I(t)", "p"); */
	
	std::cout << std::endl;
	
	std::cout << "fwdEuler - RK4 :" << std::endl;
	std::cout << "mse = " << int_lif::utils::mse(V1,V4) << std::endl;
	std::cout << "mae = " << int_lif::utils::mae(V1,V4) << std::endl;
	std::cout << std::endl;
	
	std::cout << "bwdEuler - RK4 :" << std::endl;
	std::cout << "mse = " << int_lif::utils::mse(V2,V4) << std::endl;
	std::cout << "mae = " << int_lif::utils::mae(V2,V4) << std::endl;
	std::cout << std::endl;
	
	std::cout << "Heun - RK4 :" << std::endl;
	std::cout << "mse = " << int_lif::utils::mse(V3,V4) << std::endl;
	std::cout << "mae = " << int_lif::utils::mae(V3,V4) << std::endl;
	std::cout << std::endl;
	
	canvas->cd();
	canvas->SetGrid();
	multigraph->Draw("AP");
	multigraph->SetTitle(";t [s];V [V]");
	legend->Draw();
	canvas->SaveAs("Euler_comparison.pdf");
	
}