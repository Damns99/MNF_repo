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
	
	std::vector<double> hvec = int_lif::utils::logspace(-6, -2, 400);
	std::vector<double> msevec1, msevec2, msevec3;
	std::vector<double> maevec1, maevec2, maevec3;
	
	int percent1 = 0, ii1 = -1;
	for(double h : hvec) {
		int_lif::utils::printPercent(ii1++, percent1, hvec.size(), "");
		int N = int(simtime/h);
		std::vector<double> I = int_lif::currents::pulse_train(N, int(0.1/h), int(0.01/h), aI, int(0.05/h), 25);
		std::vector<double> IRK4 = int_lif::currents::pulse_train(2*N, 2*int(0.1/h), 2*int(0.01/h), aI, 2*int(0.05/h), 25);
		std::vector<double> V1, V2, V3, V4;
		V1 = int_lif::fwdEuler(V0, h, N, I, params);
		V2 = int_lif::bwdEuler(V0, h, N, I, params);
		V3 = int_lif::Heun(V0, h, N, I, params);
		V4 = int_lif::RK4(V0, h, N, IRK4, params);
		msevec1.push_back(int_lif::utils::mse(V1,V4));
		msevec2.push_back(int_lif::utils::mse(V2,V4));
		msevec3.push_back(int_lif::utils::mse(V3,V4));
		maevec1.push_back(int_lif::utils::mae(V1,V4));
		maevec2.push_back(int_lif::utils::mae(V2,V4));
		maevec3.push_back(int_lif::utils::mae(V3,V4));
	}
	std::cout << std::endl;
	
	TCanvas* canvas1 = new TCanvas("canvas1", "Canvas1", 600, 400);
	TMultiGraph* multigraph1 = new TMultiGraph();
	TLegend* legend1 = new TLegend(0.15, 0.75, 0.25, 0.85);
	TGraph* graph1[5];
	
	graph1[0] = new TGraph(hvec.size(), hvec.data(), msevec1.data());
	graph1[0]->SetMarkerColor(1);
	graph1[0]->SetMarkerStyle(8);
	graph1[0]->SetMarkerSize(0.25);
	graph1[0]->SetLineColor(1);
	multigraph1->Add(graph1[0]);
	legend1->AddEntry(graph1[0], "fwdEuler", "p");
	
	graph1[1] = new TGraph(hvec.size(), hvec.data(), msevec2.data());
	graph1[1]->SetMarkerColor(2);
	graph1[1]->SetMarkerStyle(8);
	graph1[1]->SetMarkerSize(0.25);
	graph1[1]->SetLineColor(2);
	multigraph1->Add(graph1[1]);
	legend1->AddEntry(graph1[1], "bwdEuler", "p");
	
	graph1[2] = new TGraph(hvec.size(), hvec.data(), msevec3.data());
	graph1[2]->SetMarkerColor(3);
	graph1[2]->SetMarkerStyle(8);
	graph1[2]->SetMarkerSize(0.25);
	graph1[2]->SetLineColor(3);
	multigraph1->Add(graph1[2]);
	legend1->AddEntry(graph1[2], "Heun", "p");
	
	canvas1->cd();
	canvas1->SetGrid();
	multigraph1->Draw("AP");
	multigraph1->SetTitle(";h [s];mse [V^2]");
	legend1->Draw();
	canvas1->SaveAs("mse_comparison.pdf");
	multigraph1->GetXaxis()->SetLimits(0.5e-6, 2e-2);
	canvas1->SetLogx();
	canvas1->SaveAs("mse_comparison_lx.pdf");
	multigraph1->SetMaximum(pow(10, ceil(log10(msevec1[hvec.size()-1]))));
	multigraph1->SetMinimum(pow(10, floor(log10(msevec3[0]))));
	canvas1->SetLogy();
	canvas1->SaveAs("mse_comparison_lxly.pdf");
	
	TCanvas* canvas2 = new TCanvas("canvas2", "Canvas2", 600, 400);
	TMultiGraph* multigraph2 = new TMultiGraph();
	TLegend* legend2 = new TLegend(0.15, 0.75, 0.25, 0.85);
	TGraph* graph2[5];
	
	graph2[0] = new TGraph(hvec.size(), hvec.data(), maevec1.data());
	graph2[0]->SetMarkerColor(1);
	graph2[0]->SetMarkerStyle(8);
	graph2[0]->SetMarkerSize(0.25);
	graph2[0]->SetLineColor(1);
	multigraph2->Add(graph2[0]);
	legend2->AddEntry(graph2[0], "fwdEuler", "p");
	
	graph2[1] = new TGraph(hvec.size(), hvec.data(), maevec2.data());
	graph2[1]->SetMarkerColor(2);
	graph2[1]->SetMarkerStyle(8);
	graph2[1]->SetMarkerSize(0.25);
	graph2[1]->SetLineColor(2);
	multigraph2->Add(graph2[1]);
	legend2->AddEntry(graph2[1], "bwdEuler", "p");
	
	graph2[2] = new TGraph(hvec.size(), hvec.data(), maevec3.data());
	graph2[2]->SetMarkerColor(3);
	graph2[2]->SetMarkerStyle(8);
	graph2[2]->SetMarkerSize(0.25);
	graph2[2]->SetLineColor(3);
	multigraph2->Add(graph2[2]);
	legend2->AddEntry(graph2[2], "Heun", "p");
	
	canvas2->cd();
	canvas2->SetGrid();
	multigraph2->Draw("AP");
	multigraph2->SetTitle(";h [s];mae [V]");
	legend2->Draw();
	canvas2->SaveAs("mae_comparison.pdf");
	multigraph2->GetXaxis()->SetLimits(0.5e-6, 2e-2);
	canvas2->SetLogx();
	canvas2->SaveAs("mae_comparison_lx.pdf");
	multigraph2->SetMaximum(pow(10, ceil(log10(maevec1[hvec.size()-1]))));
	multigraph2->SetMinimum(pow(10, floor(log10(maevec3[0]))));
	canvas2->SetLogy();
	canvas2->SaveAs("mae_comparison_lxly.pdf");
	
}