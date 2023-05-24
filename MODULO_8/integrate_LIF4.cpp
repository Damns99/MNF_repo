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

#include "integrate_LIF.h"

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
	double t0 = 0; // [s]
	double V0 = -70e-3; // [V]
	
	double simtime = 1.; // [s]
	int nrep = 100;
	
	TCanvas* canvas = new TCanvas("canvas0", "Canvas0", 600, 400);
	TMultiGraph* multigraph = new TMultiGraph();
	TLegend* legend = new TLegend(0.75, 0.75, 0.85, 0.85);
	TGraphErrors* graph[5];
	
	std::vector<double> hvec = int_lif::utils::logspace(-6, -3, 100);
	
	// fwdEuler
	std::vector<double> timevec1, dtimevec1;
	int percent1 = 0, ii1 = -1;
	for(double h : hvec) {
		int_lif::utils::printPercent(ii1++, percent1, hvec.size(), "fwdEuler: ");
		int N = int(simtime/h);
		std::vector<double> I = int_lif::currents::pulse_train(N, int(0.1/h), int(0.01/h), aI, int(0.05/h), 25);
		std::vector<double> time1(nrep);
		for(int jj = 0; jj < nrep; jj++) {
			std::vector<double> V1, t1;
			std::chrono::time_point<std::chrono::steady_clock> start1, end1;
			start1 = std::chrono::steady_clock::now();
			V1 = int_lif::fwdEuler(V0, h, N, I, params);
			t1 = int_lif::utils::linspace(t0, t0 + N*h, N);
			end1 = std::chrono::steady_clock::now();
			time1[jj] = std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1).count() / 1e6; // [s]
		}
		// std::cout << "JB = " << int_lif::utils::JarqueBera(time1) << std::endl;
		timevec1.push_back(int_lif::utils::mean(time1));
		dtimevec1.push_back(int_lif::utils::stdev(time1)/sqrt(nrep));
	}
	graph[0] = new TGraphErrors(hvec.size(), hvec.data(), timevec1.data(), nullptr, dtimevec1.data());
	graph[0]->SetMarkerColor(1);
	graph[0]->SetMarkerStyle(8);
	graph[0]->SetMarkerSize(0.25);
	graph[0]->SetLineColor(1);
	multigraph->Add(graph[0]);
	legend->AddEntry(graph[0], "fwdEuler", "p");
	std::cout << std::endl;
	/* for(auto it1 = timevec1.begin(), it2 = dtimevec1.begin(); it1 != timevec1.end() && it2 != dtimevec1.end(); it1++, it2++) {
		std::cout << "[" << (*it1)-(*it2) << " " << (*it1)+(*it2) << "]" << std::endl;
	}
	std::cout << std::endl; */
	
	// bwdEuler
	std::vector<double> timevec2, dtimevec2;
	int percent2 = 0, ii2 = -1;
	for(double h : hvec) {
		int_lif::utils::printPercent(ii2++, percent2, hvec.size(), "bwdEuler: ");
		int N = int(simtime/h);
		std::vector<double> I = int_lif::currents::pulse_train(N, int(0.1/h), int(0.01/h), aI, int(0.05/h), 25);
		std::vector<double> time2(nrep);
		for(int jj = 0; jj < nrep; jj++) {
			std::vector<double> V2, t2;
			std::chrono::time_point<std::chrono::steady_clock> start2, end2;
			start2 = std::chrono::steady_clock::now();
			V2 = int_lif::bwdEuler(V0, h, N, I, params);
			t2 = int_lif::utils::linspace(t0, t0 + N*h, N);
			end2 = std::chrono::steady_clock::now();
			time2[jj] = std::chrono::duration_cast<std::chrono::microseconds>(end2 - start2).count() / 1e6; // [s]
		}
		timevec2.push_back(int_lif::utils::mean(time2));
		dtimevec2.push_back(int_lif::utils::stdev(time2)/sqrt(nrep));
	}
	graph[1] = new TGraphErrors(hvec.size(), hvec.data(), timevec2.data(), nullptr, dtimevec2.data());
	graph[1]->SetMarkerColor(2);
	graph[1]->SetMarkerStyle(8);
	graph[1]->SetMarkerSize(0.25);
	graph[1]->SetLineColor(2);
	multigraph->Add(graph[1]);
	legend->AddEntry(graph[1], "bwdEuler", "p");
	std::cout << std::endl;
	/* for(auto it1 = timevec2.begin(), it2 = dtimevec2.begin(); it1 != timevec2.end() && it2 != dtimevec2.end(); it1++, it2++) {
		std::cout << "[" << (*it1)-(*it2) << " " << (*it1)+(*it2) << "]" << std::endl;
	}
	std::cout << std::endl; */
	
	// Heun
	std::vector<double> timevec3, dtimevec3;
	int percent3 = 0, ii3 = -1;
	for(double h : hvec) {
		int_lif::utils::printPercent(ii3++, percent3, hvec.size(), "Heun: ");
		int N = int(simtime/h);
		std::vector<double> I = int_lif::currents::pulse_train(N, int(0.1/h), int(0.01/h), aI, int(0.05/h), 15);
		std::vector<double> time3(nrep);
		for(int jj = 0; jj < nrep; jj++) {
			std::vector<double> V3, t3;
			std::chrono::time_point<std::chrono::steady_clock> start3, end3;
			start3 = std::chrono::steady_clock::now();
			V3 = int_lif::Heun(V0, h, N, I, params);
			t3 = int_lif::utils::linspace(t0, t0 + N*h, N);
			end3 = std::chrono::steady_clock::now();
			time3[jj] = std::chrono::duration_cast<std::chrono::microseconds>(end3 - start3).count() / 1e6; // [s]
		}
		timevec3.push_back(int_lif::utils::mean(time3));
		dtimevec3.push_back(int_lif::utils::stdev(time3)/sqrt(nrep));
	}
	graph[2] = new TGraphErrors(hvec.size(), hvec.data(), timevec3.data(), nullptr, dtimevec3.data());
	graph[2]->SetMarkerColor(3);
	graph[2]->SetMarkerStyle(8);
	graph[2]->SetMarkerSize(0.25);
	graph[2]->SetLineColor(3);
	multigraph->Add(graph[2]);
	legend->AddEntry(graph[2], "Heun", "p");
	std::cout << std::endl;
	/* for(auto it1 = timevec3.begin(), it2 = dtimevec3.begin(); it1 != timevec3.end() && it2 != dtimevec3.end(); it1++, it2++) {
		std::cout << "[" << (*it1)-(*it2) << " " << (*it1)+(*it2) << "]" << std::endl;
	}
	std::cout << std::endl; */
	
	// RK4
	std::vector<double> timevec4, dtimevec4;
	int percent4 = 0, ii4 = -1;
	for(double h : hvec) {
		int_lif::utils::printPercent(ii4++, percent4, hvec.size(), "RK4: ");
		int N = int(simtime/h);
		std::vector<double> I = int_lif::currents::pulse_train(2*N, 2*int(0.1/h), 2*int(0.01/h), aI, 2*int(0.05/h), 15);
		std::vector<double> time4(nrep);
		for(int jj = 0; jj < nrep; jj++) {
			std::vector<double> V4, t4;
			std::chrono::time_point<std::chrono::steady_clock> start4, end4;
			start4 = std::chrono::steady_clock::now();
			V4 = int_lif::RK4(V0, h, N, I, params);
			t4 = int_lif::utils::linspace(t0, t0 + N*h, N);
			end4 = std::chrono::steady_clock::now();
			time4[jj] = std::chrono::duration_cast<std::chrono::microseconds>(end4 - start4).count() / 1e6; // [s]
		}
		timevec4.push_back(int_lif::utils::mean(time4));
		dtimevec4.push_back(int_lif::utils::stdev(time4)/sqrt(nrep));
	}
	graph[3] = new TGraphErrors(hvec.size(), hvec.data(), timevec4.data(), nullptr, dtimevec4.data());
	graph[3]->SetMarkerColor(4);
	graph[3]->SetMarkerStyle(8);
	graph[3]->SetMarkerSize(0.25);
	graph[3]->SetLineColor(4);
	multigraph->Add(graph[3]);
	legend->AddEntry(graph[3], "RK4", "p");
	std::cout << std::endl;
	/* for(auto it1 = timevec4.begin(), it2 = dtimevec4.begin(); it1 != timevec4.end() && it2 != dtimevec4.end(); it1++, it2++) {
		std::cout << "[" << (*it1)-(*it2) << " " << (*it1)+(*it2) << "]" << std::endl;
	}
	std::cout << std::endl; */
	
	canvas->cd();
	canvas->SetGrid();
	multigraph->Draw("AP");
	multigraph->SetTitle(";h [s];time [s]");
	legend->Draw();
	canvas->SaveAs("timings_comparison.pdf");
	multigraph->GetXaxis()->SetLimits(0.5e-6, 2e-3);
	canvas->SetLogx();
	canvas->SaveAs("timings_comparison_lx.pdf");
	multigraph->SetMaximum(pow(10, ceil(log10(timevec4[0]))));
	multigraph->SetMinimum(pow(10, floor(log10(timevec1[hvec.size()-1]))));
	multigraph->Draw("APX"); // no error bars in logy
	legend->Draw();
	canvas->SetLogy();
	canvas->SaveAs("timings_comparison_lxly.pdf");
	
}