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
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TFrame.h>
#include <TLegend.h>

#include "mystyle.h"
#include "integrate_LIF.h"

void addToMultigraph(TMultiGraph* multigraph, TLegend* legend, std::vector<double>& x, std::vector<double>& y, std::vector<double>& dy, int n, int color, const char* name, const char* labeltype = "p", const char* drawoption = "") {
	TGraphErrors* graph = new TGraphErrors(n, x.data(), y.data(), nullptr, dy.data());
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
	double simtime = 1e2;
	double az = -0.16;
	// double startz = 5, durz = 50, intz = 20;
	// int numz = 25;
	// std::vector<double> z = int_lif::currents::pulse_train(N, int(startz/h), int(durz/h), az, int(intz/h), numz);
	double periodz = 40, phasez = 3.923, offsetz = az;
	double x0 = 0.;
	double y0 = 1.;
	
	int nrep = 100;
	
	myStyle();
	
	TCanvas* canvas = new TCanvas("canvas0", "Canvas0", 600, 400);
	TMultiGraph* multigraph = new TMultiGraph();
	TLegend* legend = new TLegend(0.75, 0.75, 0.85, 0.85);
	
	std::vector<double> hvec = int_lif::utils::logspace(-4, 0, 400);
	
	// fwdEuler
	std::vector<double> timevec1, dtimevec1;
	int percent1 = 0, ii1 = -1;
	for(double h : hvec) {
		int_lif::utils::printPercent(ii1++, percent1, hvec.size(), "fwdEuler: ");
		int N = int(simtime/h);
		std::vector<double> z = int_lif::currents::sine_wave(N, int(periodz/h), az, int(phasez/h), offsetz);
		std::vector<double> time(nrep);
		for(int jj = 0; jj < nrep; jj++) {
			std::vector<double> y, x;
			std::chrono::time_point<std::chrono::steady_clock> start, end;
			start = std::chrono::steady_clock::now();
			y = int_lif::fwdEuler(y0, h, N, z, params);
			end = std::chrono::steady_clock::now();
			time[jj] = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1.e9; // [s]
			x = int_lif::utils::linspace(x0, x0 + N*h, N);
		}
		timevec1.push_back(int_lif::utils::mean(time));
		dtimevec1.push_back(int_lif::utils::stdev(time)/sqrt(nrep));
	}
	std::cout << std::endl;
	
	// bwdEuler
	std::vector<double> timevec2, dtimevec2;
	int percent2 = 0, ii2 = -1;
	for(double h : hvec) {
		int_lif::utils::printPercent(ii2++, percent2, hvec.size(), "bwdEuler: ");
		int N = int(simtime/h);
		std::vector<double> z = int_lif::currents::sine_wave(N, int(periodz/h), az, int(phasez/h), offsetz);
		std::vector<double> time(nrep);
		for(int jj = 0; jj < nrep; jj++) {
			std::vector<double> y, x;
			std::chrono::time_point<std::chrono::steady_clock> start, end;
			start = std::chrono::steady_clock::now();
			y = int_lif::bwdEuler(y0, h, N, z, params);
			x = int_lif::utils::linspace(x0, x0 + N*h, N);
			end = std::chrono::steady_clock::now();
			time[jj] = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1.e9; // [s]
		}
		timevec2.push_back(int_lif::utils::mean(time));
		dtimevec2.push_back(int_lif::utils::stdev(time)/sqrt(nrep));
	}
	std::cout << std::endl;
	
	// Heun
	std::vector<double> timevec3, dtimevec3;
	int percent3 = 0, ii3 = -1;
	for(double h : hvec) {
		int_lif::utils::printPercent(ii3++, percent3, hvec.size(), "Heun: ");
		int N = int(simtime/h);
		std::vector<double> z = int_lif::currents::sine_wave(N, int(periodz/h), az, int(phasez/h), offsetz);
		std::vector<double> time(nrep);
		for(int jj = 0; jj < nrep; jj++) {
			std::vector<double> y, x;
			std::chrono::time_point<std::chrono::steady_clock> start, end;
			start = std::chrono::steady_clock::now();
			y = int_lif::Heun(y0, h, N, z, params);
			x = int_lif::utils::linspace(x0, x0 + N*h, N);
			end = std::chrono::steady_clock::now();
			time[jj] = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1.e9; // [s]
		}
		timevec3.push_back(int_lif::utils::mean(time));
		dtimevec3.push_back(int_lif::utils::stdev(time)/sqrt(nrep));
	}
	std::cout << std::endl;
	
	// RK4
	std::vector<double> timevec4, dtimevec4;
	int percent4 = 0, ii4 = -1;
	for(double h : hvec) {
		int_lif::utils::printPercent(ii4++, percent4, hvec.size(), "RK4: ");
		int N = int(simtime/h);
		std::vector<double> zRK4 = int_lif::currents::sine_wave(2*N, 2*int(periodz/h), az, 2*int(phasez/h), offsetz);
		std::vector<double> time(nrep);
		for(int jj = 0; jj < nrep; jj++) {
			std::vector<double> y, x;
			std::chrono::time_point<std::chrono::steady_clock> start, end;
			start = std::chrono::steady_clock::now();
			y = int_lif::RK4(y0, h, N, zRK4, params);
			x = int_lif::utils::linspace(x0, x0 + N*h, N);
			end = std::chrono::steady_clock::now();
			time[jj] = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1.e9; // [s]
		}
		timevec4.push_back(int_lif::utils::mean(time));
		dtimevec4.push_back(int_lif::utils::stdev(time)/sqrt(nrep));
	}
	std::cout << std::endl;
	
	addToMultigraph(multigraph, legend, hvec, timevec1, dtimevec1, hvec.size(), 1, "fwdEuler", "p", "");
	addToMultigraph(multigraph, legend, hvec, timevec2, dtimevec2, hvec.size(), 2, "bwdEuler", "p", "");
	addToMultigraph(multigraph, legend, hvec, timevec3, dtimevec3, hvec.size(), 3, "Heun", "p", "");
	addToMultigraph(multigraph, legend, hvec, timevec4, dtimevec4, hvec.size(), 4, "RK4", "p", "");
	
	canvas->cd();
	canvas->SetGrid();
	multigraph->Draw("AP");
	multigraph->SetTitle("timings comparison;h [s];time [s]");
	legend->Draw();
	canvas->SaveAs("timings_comparison.pdf");
	multigraph->GetXaxis()->SetLimits(0.5e-4, 2e0);
	canvas->SetLogx();
	canvas->SaveAs("timings_comparison_lx.pdf");
	multigraph->SetMaximum(pow(10, ceil(log10(multigraph->GetHistogram()->GetMaximum()))));
	multigraph->SetMinimum(pow(10, -7));
	multigraph->Draw("APX"); // no error bars in logy
	legend->Draw();
	canvas->SetLogy();
	canvas->SaveAs("timings_comparison_lxly.pdf");
	
}