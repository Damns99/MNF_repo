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
	double params[6] = {
		1e-10, // Cm [s/ohm]
		1e-8, // g [1/ohm]
		-70e-3, // Vrest [V]
		1., // Vth [V]
		-70e-3,  // Vreset [V]
		1e-3 // taurefr [s]
	};
	double aI = 1e-9; // [A]
	double startI = 0.1, durI = 0.008, intI = 0.05;
	int numI = 25;
	double t0 = 0; // [s]
	double V0 = -70e-3; // [V]
	
	double simtime = 1.; // [s]
	int nrep = 100;
	
	myStyle();
	
	TCanvas* canvas = new TCanvas("canvas0", "Canvas0", 600, 400);
	TMultiGraph* multigraph = new TMultiGraph();
	TLegend* legend = new TLegend(0.75, 0.75, 0.85, 0.85);
	
	std::vector<double> hvec = int_lif::utils::logspace(-4, -3, 100);
	
	// fwdEuler
	std::vector<double> timevec1, dtimevec1;
	int percent1 = 0, ii1 = -1;
	for(double h : hvec) {
		int_lif::utils::printPercent(ii1++, percent1, hvec.size(), "fwdEuler: ");
		int N = int(simtime/h);
		std::vector<double> I = int_lif::currents::pulse_train(N, int(startI/h), int(durI/h), aI, int(intI/h), numI);
		std::vector<double> time(nrep);
		for(int jj = 0; jj < nrep; jj++) {
			std::vector<double> V, t;
			std::chrono::time_point<std::chrono::steady_clock> start, end;
			start = std::chrono::steady_clock::now();
			V = int_lif::fwdEuler(V0, h, N, I, params);
			t = int_lif::utils::linspace(t0, t0 + N*h, N);
			end = std::chrono::steady_clock::now();
			time[jj] = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1.e9; // [s]
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
		std::vector<double> I = int_lif::currents::pulse_train(N, int(startI/h), int(durI/h), aI, int(intI/h), numI);
		std::vector<double> time(nrep);
		for(int jj = 0; jj < nrep; jj++) {
			std::vector<double> V, t;
			std::chrono::time_point<std::chrono::steady_clock> start, end;
			start = std::chrono::steady_clock::now();
			V = int_lif::bwdEuler(V0, h, N, I, params);
			t = int_lif::utils::linspace(t0, t0 + N*h, N);
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
		std::vector<double> I = int_lif::currents::pulse_train(N, int(startI/h), int(durI/h), aI, int(intI/h), numI);
		std::vector<double> time(nrep);
		for(int jj = 0; jj < nrep; jj++) {
			std::vector<double> V, t;
			std::chrono::time_point<std::chrono::steady_clock> start, end;
			start = std::chrono::steady_clock::now();
			V = int_lif::Heun(V0, h, N, I, params);
			t = int_lif::utils::linspace(t0, t0 + N*h, N);
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
		std::vector<double> I = int_lif::currents::pulse_train(2*N, 2*int(startI/h), 2*int(durI/h), aI, 2*int(intI/h), numI);
		std::vector<double> time(nrep);
		for(int jj = 0; jj < nrep; jj++) {
			std::vector<double> V, t;
			std::chrono::time_point<std::chrono::steady_clock> start, end;
			start = std::chrono::steady_clock::now();
			V = int_lif::RK4(V0, h, N, I, params);
			t = int_lif::utils::linspace(t0, t0 + N*h, N);
			end = std::chrono::steady_clock::now();
			time[jj] = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1.e9; // [s]
		}
		timevec4.push_back(int_lif::utils::mean(time));
		dtimevec4.push_back(int_lif::utils::stdev(time)/sqrt(nrep));
	}
	std::cout << std::endl;
	
	// empty
	std::vector<double> timevec5, dtimevec5;
	int percent5 = 0, ii5 = -1;
	for(double h : hvec) {
		int_lif::utils::printPercent(ii5++, percent5, hvec.size(), "empty: ");
		int N = int(simtime/h);
		std::vector<double> I = int_lif::currents::pulse_train(N, int(startI/h), int(durI/h), aI, int(intI/h), numI);
		std::vector<double> time(nrep);
		for(int jj = 0; jj < nrep; jj++) {
			std::vector<double> V, t;
			std::chrono::time_point<std::chrono::steady_clock> start, end;
			start = std::chrono::steady_clock::now();
			V = int_lif::empty(V0, h, N, I, params);
			t = int_lif::utils::linspace(t0, t0 + N*h, N);
			end = std::chrono::steady_clock::now();
			time[jj] = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1.e9; // [s]
		}
		timevec5.push_back(int_lif::utils::mean(time));
		dtimevec5.push_back(int_lif::utils::stdev(time)/sqrt(nrep));
	}
	std::cout << std::endl;
	
	// empty2
	std::vector<double> timevec6, dtimevec6;
	int percent6 = 0, ii6 = -1;
	for(double h : hvec) {
		int_lif::utils::printPercent(ii6++, percent6, hvec.size(), "empty2: ");
		int N = int(simtime/h);
		std::vector<double> I = int_lif::currents::pulse_train(N, int(startI/h), int(durI/h), aI, int(intI/h), numI);
		std::vector<double> time(nrep);
		for(int jj = 0; jj < nrep; jj++) {
			std::vector<double> V, t;
			std::chrono::time_point<std::chrono::steady_clock> start, end;
			start = std::chrono::steady_clock::now();
			V = int_lif::empty2(V0, h, N, I, params);
			t = int_lif::utils::linspace(t0, t0 + N*h, N);
			end = std::chrono::steady_clock::now();
			time[jj] = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1.e9; // [s]
		}
		timevec6.push_back(int_lif::utils::mean(time));
		dtimevec6.push_back(int_lif::utils::stdev(time)/sqrt(nrep));
	}
	std::cout << std::endl;
	
	addToMultigraph(multigraph, legend, hvec, timevec1, dtimevec1, hvec.size(), 1, "fwdEuler", "p", "");
	addToMultigraph(multigraph, legend, hvec, timevec2, dtimevec2, hvec.size(), 2, "bwdEuler", "p", "");
	addToMultigraph(multigraph, legend, hvec, timevec3, dtimevec3, hvec.size(), 3, "Heun", "p", "");
	addToMultigraph(multigraph, legend, hvec, timevec4, dtimevec4, hvec.size(), 4, "RK4", "p", "");
	addToMultigraph(multigraph, legend, hvec, timevec5, dtimevec5, hvec.size(), 5, "empty", "p", "");
	addToMultigraph(multigraph, legend, hvec, timevec6, dtimevec6, hvec.size(), 6, "empty2", "p", "");
	
	canvas->cd();
	canvas->SetGrid();
	multigraph->Draw("AP");
	multigraph->SetTitle("timings comparison;h [s];time [s]");
	legend->Draw();
	canvas->SaveAs("timings_comparison.pdf");
	multigraph->GetXaxis()->SetLimits(0.5e-6, 2e-3);
	canvas->SetLogx();
	canvas->SaveAs("timings_comparison_lx.pdf");
	multigraph->SetMaximum(pow(10, ceil(log10(multigraph->GetHistogram()->GetMaximum()))));
	multigraph->SetMinimum(pow(10, -7));
	multigraph->Draw("APX"); // no error bars in logy
	legend->Draw();
	canvas->SetLogy();
	canvas->SaveAs("timings_comparison_lxly.pdf");
	
}