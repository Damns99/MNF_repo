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

#include "mystyle.h"
#include "integrate_LIF.h"

void addToMultigraph(TMultiGraph* multigraph, TLegend* legend, std::vector<double>& x, std::vector<double>& y, int n, int color, const char* name, const char* labeltype = "p", const char* drawoption = "") {
	TGraph* graph = new TGraph(n, x.data(), y.data());
	graph->SetMarkerColor(color);
	graph->SetLineColor(color);
	multigraph->Add(graph, drawoption);
	legend->AddEntry(graph, name, labeltype);
}

int main() {
	double params[6] = {
		1e-10, // Cm [ms/ohm]
		1e-8, // g [1/ohm]
		-70e-3, // Vrest [V]
		-50e-3, // Vth [V]
		-70e-3,  // Vreset [V]
		1e-3 // taurefr [s]
	};
	double h = 1e-3; // [s]
	double simtime = 1.; // [s]
	int N = int(simtime/h);
	double aI = 0.1e-9; // [A]
	double startI = 0.05, durI = 0.05, intI = 0.02;
	int numI = 25;
	std::vector<double> I = int_lif::currents::pulse_train(N, int(startI/h), int(durI/h), aI, int(intI/h), numI);
	double t0 = 0; // [s]
	double V0 = -70e-3; // [V]
	
	int nrep = 100;
	
	myStyle();
	
	TCanvas* canvas = new TCanvas("canvas0", "Canvas0", 600, 400);
	TMultiGraph* multigraph = new TMultiGraph();
	TLegend* legend = new TLegend(0.75, 0.75, 0.85, 0.85);
	
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
	
	// RK4
	std::vector<double> IRK4 = int_lif::currents::pulse_train(2*N, 2*int(startI/h), 2*int(durI/h), aI, 2*int(intI/h), numI);
	std::chrono::time_point<std::chrono::steady_clock> start4, end4;
	start4 = std::chrono::steady_clock::now();
	std::vector<double> V4, t4;
	for(int jj = 0; jj < nrep; jj++) {
		V4 = int_lif::RK4(V0, h, N, IRK4, params);
		t4 = int_lif::utils::linspace(t0, t0 + N*h, N);
	}
	end4 = std::chrono::steady_clock::now();
	std::cout << "RK4: " << std::chrono::duration_cast<std::chrono::microseconds>(end4 - start4).count() / 1000. << " ms" << std::endl;
	
	std::cout << std::endl;
	
	// Benchmark
	double h_bench = 1e-6;
	int N_bench = int(simtime/h_bench);
	std::vector<double> I_bench = int_lif::currents::pulse_train(2*N_bench, 2*int(startI/h_bench), 2*int(durI/h_bench), aI, 2*int(intI/h_bench), numI);
	std::vector<double> V_bench_tmp, t_bench_tmp;
	V_bench_tmp = int_lif::RK4(V0, h_bench, N_bench, I_bench, params);
	t_bench_tmp = int_lif::utils::linspace(t0, t0 + N_bench*h_bench, N_bench);
	
	std::vector<double> V_bench(N);
	for(int i = 0; i < N; i++) V_bench[i] = V_bench_tmp[int(h/h_bench*i)];
	
	std::cout << "fwdEuler - benchmark :" << std::endl;
	std::cout << "mse = " << int_lif::utils::mse(V1,V_bench) << std::endl;
	std::cout << "mae = " << int_lif::utils::mae(V1,V_bench) << std::endl;
	std::cout << std::endl;
	
	std::cout << "bwdEuler - benchmark :" << std::endl;
	std::cout << "mse = " << int_lif::utils::mse(V2,V_bench) << std::endl;
	std::cout << "mae = " << int_lif::utils::mae(V2,V_bench) << std::endl;
	std::cout << std::endl;
	
	std::cout << "Heun - benchmark :" << std::endl;
	std::cout << "mse = " << int_lif::utils::mse(V3,V_bench) << std::endl;
	std::cout << "mae = " << int_lif::utils::mae(V3,V_bench) << std::endl;
	std::cout << std::endl;
	
	std::cout << "RK4 - benchmark :" << std::endl;
	std::cout << "mse = " << int_lif::utils::mse(V4,V_bench) << std::endl;
	std::cout << "mae = " << int_lif::utils::mae(V4,V_bench) << std::endl;
	std::cout << std::endl;
	
	addToMultigraph(multigraph, legend, t_bench_tmp, V_bench_tmp, N_bench, kBlack, "benchmark", "l", "L");
	addToMultigraph(multigraph, legend, t1, V1, N, 1, "fwdEuler", "p", "");
	addToMultigraph(multigraph, legend, t2, V2, N, 2, "bwdEuler", "p", "");
	addToMultigraph(multigraph, legend, t3, V3, N, 3, "Heun", "p", "");
	addToMultigraph(multigraph, legend, t4, V4, N, 4, "RK4", "p", "");
	addToMultigraph(multigraph, legend, t4, V_bench, N, 5, "V_bench", "p", "");
	
	canvas->cd();
	canvas->SetGrid();
	multigraph->Draw("AP");
	multigraph->SetTitle("LIF integrators comparison;t [s];V [V]");
	legend->Draw();
	canvas->SaveAs("integration_comparison.pdf");
	
	TCanvas* canvas2 = new TCanvas("canvas2", "Canvas2", 600, 400);
	TMultiGraph* multigraph2 = new TMultiGraph();
	TLegend* legend2 = new TLegend(0.75, 0.75, 0.85, 0.85);
	TGraph* graph2[8];
	
	std::vector<double> diff1(N,0.), diff2(N,0.), diff3(N,0.), diff4(N,0.);
	for(int i = 0; i < N; i++) {
		if(abs(V_bench[i] - V0) != 0) {
			diff1[i] = abs(V1[i] - V_bench[i]) / abs(V_bench[i] - V0);
			diff2[i] = abs(V2[i] - V_bench[i]) / abs(V_bench[i] - V0);
			diff3[i] = abs(V3[i] - V_bench[i]) / abs(V_bench[i] - V0);
			diff4[i] = abs(V4[i] - V_bench[i]) / abs(V_bench[i] - V0);
		}
	}
	
	addToMultigraph(multigraph2, legend2, t1, diff1, N, 1, "fwdEuler", "p", "");
	addToMultigraph(multigraph2, legend2, t2, diff2, N, 2, "bwdEuler", "p", "");
	addToMultigraph(multigraph2, legend2, t3, diff3, N, 3, "Heun", "p", "");
	addToMultigraph(multigraph2, legend2, t4, diff4, N, 4, "RK4", "p", "");
	
	canvas2->cd();
	canvas2->SetGrid();
	multigraph2->SetMinimum(1e-4);
	multigraph2->Draw("AP");
	multigraph2->SetTitle("LIF integrators deviation from reference;t [s];abs of relative deviation []");
	legend2->Draw();
	canvas2->SetLogy();
	canvas2->SaveAs("deviation_comparison.pdf");
	
}