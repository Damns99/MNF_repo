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
	// Voltages in [Vrest = -70e-3 V], Currents in [g*Vrest = -7e-10 A], Times in [Cm/g = 1e-2 s]
	double params[2] = {
		0.7, // yth
		1.,  // yreset
	};
	double h = 1e-1;
	double simtime = 1e2;
	int N = int(simtime/h);
	double az = -0.16;
	// double startz = 5, durz = 50, intz = 20;
	// int numz = 25;
	// std::vector<double> z = int_lif::currents::pulse_train(N, int(startz/h), int(durz/h), az, int(intz/h), numz);
	double periodz = 40, phasez = 3.923, offsetz = az;
	std::vector<double> z = int_lif::currents::sine_wave(N, int(periodz/h), az, int(phasez/h), offsetz);
	double x0 = 0.;
	double y0 = 1.;
	
	int nrep = 100;
	
	myStyle();
	
	TCanvas* canvas = new TCanvas("canvas0", "Canvas0", 600, 400);
	TMultiGraph* multigraph = new TMultiGraph();
	TLegend* legend = new TLegend(0.75, 0.75, 0.85, 0.85);
	
	// fwdEuler
	std::chrono::time_point<std::chrono::steady_clock> start1, end1;
	start1 = std::chrono::steady_clock::now();
	std::vector<double> y1, x1, spkt1;
	for(int jj = 0; jj < nrep; jj++) {
		y1 = int_lif::fwdEuler(y0, h, N, z, params);
		x1 = int_lif::utils::linspace(x0, x0 + N*h, N);
	}
	end1 = std::chrono::steady_clock::now();
	std::cout << "fwdEuler: " << std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1).count() / 1000. << " ms" << std::endl;
	y1 = int_lif::fwdEuler(y0, h, N, z, params, &spkt1);
	std::cout << "total spikes: " << spkt1.size() << std::endl;
	std::cout << "spike times [Cm/g]";
	for(auto& sp: spkt1) std::cout << " " << sp;
	std::cout << std::endl;
	
	// bwdEuler
	std::chrono::time_point<std::chrono::steady_clock> start2, end2;
	start2 = std::chrono::steady_clock::now();
	std::vector<double> y2, x2, spkt2;
	for(int jj = 0; jj < nrep; jj++) {
		y2 = int_lif::bwdEuler(y0, h, N, z, params);
		x2 = int_lif::utils::linspace(x0, x0 + N*h, N);
	}
	end2 = std::chrono::steady_clock::now();
	std::cout << "bwdEuler: " << std::chrono::duration_cast<std::chrono::microseconds>(end2 - start2).count() / 1000. << " ms" << std::endl;
	y2 = int_lif::bwdEuler(y0, h, N, z, params, &spkt2);
	std::cout << "total spikes: " << spkt2.size() << std::endl;
	std::cout << "spike times [Cm/g]";
	for(auto& sp: spkt2) std::cout << " " << sp;
	std::cout << std::endl;
	
	// Heun
	std::chrono::time_point<std::chrono::steady_clock> start3, end3;
	start3 = std::chrono::steady_clock::now();
	std::vector<double> y3, x3, spkt3;
	for(int jj = 0; jj < nrep; jj++) {
		y3 = int_lif::Heun(y0, h, N, z, params);
		x3 = int_lif::utils::linspace(x0, x0 + N*h, N);
	}
	end3 = std::chrono::steady_clock::now();
	std::cout << "Heun: " << std::chrono::duration_cast<std::chrono::microseconds>(end3 - start3).count() / 1000. << " ms" << std::endl;
	y3 = int_lif::Heun(y0, h, N, z, params, &spkt3);
	std::cout << "total spikes: " << spkt3.size() << std::endl;
	std::cout << "spike times [Cm/g]";
	for(auto& sp: spkt3) std::cout << " " << sp;
	std::cout << std::endl;
	
	// RK4
	std::vector<double> zRK4 = int_lif::currents::sine_wave(2*N, 2*int(periodz/h), az, 2*int(phasez/h), offsetz);
	std::chrono::time_point<std::chrono::steady_clock> start4, end4;
	start4 = std::chrono::steady_clock::now();
	std::vector<double> y4, x4, spkt4;
	for(int jj = 0; jj < nrep; jj++) {
		y4 = int_lif::RK4(y0, h, N, zRK4, params);
		x4 = int_lif::utils::linspace(x0, x0 + N*h, N);
	}
	end4 = std::chrono::steady_clock::now();
	std::cout << "RK4: " << std::chrono::duration_cast<std::chrono::microseconds>(end4 - start4).count() / 1000. << " ms" << std::endl;
	y4 = int_lif::RK4(y0, h, N, zRK4, params, &spkt4);
	std::cout << "total spikes: " << spkt4.size() << std::endl;
	std::cout << "spike times [Cm/g]";
	for(auto& sp: spkt4) std::cout << " " << sp;
	std::cout << std::endl;
	
	std::cout << std::endl;
	
	// Benchmark
	double h_bench = 1e-4;
	int N_bench = int(simtime/h_bench);
	std::vector<double> z_bench = int_lif::currents::sine_wave(2*N_bench, 2*int(periodz/h_bench), az, 2*int(phasez/h_bench), offsetz);
	std::vector<double> y_bench_tmp, x_bench_tmp, spkt_bench;
	y_bench_tmp = int_lif::RK4(y0, h_bench, N_bench, z_bench, params, &spkt_bench);
	x_bench_tmp = int_lif::utils::linspace(x0, x0 + N_bench*h_bench, N_bench);
	std::vector<double> x_z_bench = int_lif::utils::linspace(x0, x0 + N_bench*h_bench, 2*N_bench);
	std::cout << "Benchmark" << std::endl;
	std::cout << "total spikes: " << spkt_bench.size() << std::endl;
	std::cout << "spike times [Cm/g]";
	for(auto& sp: spkt_bench) std::cout << " " << sp;
	std::cout << std::endl << std::endl;
	
	std::vector<double> y_bench(N);
	for(int i = 0; i < N; i++) y_bench[i] = y_bench_tmp[int(h/h_bench*i)];
	/* for(int i = 0; i < N; i++) {
		if(i*h < startz) y_bench[i] = y0;
		else if(i*h < startz+durz) y_bench[i] = y0+az*(1-exp(-(i*h-startz)));
		else if(i*h < startz+durz+intz) y_bench[i] = y0+az*(exp(-(i*h-startz-durz)));
		else y_bench[i] = y0+az*(1-exp(-(i*h-startz-durz-intz)));
	} */
	
	std::cout << "fwdEuler - benchmark :" << std::endl;
	std::cout << "mse = " << int_lif::utils::mse(y1,y_bench) << std::endl;
	std::cout << "mae = " << int_lif::utils::mae(y1,y_bench) << std::endl;
	std::cout << std::endl;
	
	std::cout << "bwdEuler - benchmark :" << std::endl;
	std::cout << "mse = " << int_lif::utils::mse(y2,y_bench) << std::endl;
	std::cout << "mae = " << int_lif::utils::mae(y2,y_bench) << std::endl;
	std::cout << std::endl;
	
	std::cout << "Heun - benchmark :" << std::endl;
	std::cout << "mse = " << int_lif::utils::mse(y3,y_bench) << std::endl;
	std::cout << "mae = " << int_lif::utils::mae(y3,y_bench) << std::endl;
	std::cout << std::endl;
	
	std::cout << "RK4 - benchmark :" << std::endl;
	std::cout << "mse = " << int_lif::utils::mse(y4,y_bench) << std::endl;
	std::cout << "mae = " << int_lif::utils::mae(y4,y_bench) << std::endl;
	std::cout << std::endl;
	
	addToMultigraph(multigraph, legend, x_bench_tmp, y_bench_tmp, N_bench, kBlack, "benchmark", "l", "L");
	// addToMultigraph(multigraph, legend, x_z_bench, z_bench, 2*N_bench, kBlack+2, "input", "l", "L");
	addToMultigraph(multigraph, legend, x1, y1, N, 1, "fwdEuler", "p", "");
	addToMultigraph(multigraph, legend, x2, y2, N, 2, "bwdEuler", "p", "");
	addToMultigraph(multigraph, legend, x3, y3, N, 3, "Heun", "p", "");
	addToMultigraph(multigraph, legend, x4, y4, N, 4, "RK4", "p", "");
	addToMultigraph(multigraph, legend, x4, y_bench, N, 5, "y_bench", "p", "");
	
	canvas->cd();
	canvas->SetGrid();
	multigraph->Draw("AP");
	multigraph->SetTitle("LIF integrators comparison;t [Cm/g];V [Vrest]");
	legend->Draw();
	canvas->SaveAs("integration_comparison.pdf");
	
	TCanvas* canvas2 = new TCanvas("canvas2", "Canvas2", 600, 400);
	TMultiGraph* multigraph2 = new TMultiGraph();
	TLegend* legend2 = new TLegend(0.75, 0.75, 0.85, 0.85);
	TGraph* graph2[8];
	
	std::vector<double> diff1(N,0.), diff2(N,0.), diff3(N,0.), diff4(N,0.);
	for(int i = 0; i < N; i++) {
		if(abs(y_bench[i] - y0) != 0) {
			diff1[i] = abs(y1[i] - y_bench[i]) / abs(y_bench[i] - y0);
			diff2[i] = abs(y2[i] - y_bench[i]) / abs(y_bench[i] - y0);
			diff3[i] = abs(y3[i] - y_bench[i]) / abs(y_bench[i] - y0);
			diff4[i] = abs(y4[i] - y_bench[i]) / abs(y_bench[i] - y0);
		}
	}
	
	addToMultigraph(multigraph2, legend2, x1, diff1, N, 1, "fwdEuler", "p", "");
	addToMultigraph(multigraph2, legend2, x2, diff2, N, 2, "bwdEuler", "p", "");
	addToMultigraph(multigraph2, legend2, x3, diff3, N, 3, "Heun", "p", "");
	addToMultigraph(multigraph2, legend2, x4, diff4, N, 4, "RK4", "p", "");
	
	canvas2->cd();
	canvas2->SetGrid();
	multigraph2->SetMinimum(1e-6);
	multigraph2->Draw("AP");
	multigraph2->SetTitle("LIF integrators deviation from reference;t [Cm/g];abs of relative deviation []");
	legend2->Draw();
	canvas2->SetLogy();
	canvas2->SaveAs("deviation_comparison.pdf");
	
}