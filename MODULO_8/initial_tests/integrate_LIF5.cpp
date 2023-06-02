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
#include <TH1F.h>

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
		1e-10, // Cm [s/ohm]
		1e-8, // g [1/ohm]
		-70e-3, // Vrest [V]
		1., // Vth [V]
		-70e-3,  // Vreset [V]
		1e-3 // taurefr [s]
	};
	double aI = 0.1e-9; // [A]
	double startI = 0.05, durI = 0.5, intI = 0.2;
	int numI = 25;
	double t0 = 0; // [s]
	double V0 = -70e-3; // [V]
	
	double simtime = 1.; // [s]
	
	std::vector<double> hvec = int_lif::utils::logspace(-7, -2, 400);
	std::vector<double> msevec1, msevec2, msevec3, msevec4;
	std::vector<double> maevec1, maevec2, maevec3, maevec4;
	
	// Benchmark
	double h_bench = 1e-7;
	int N_bench = int(simtime/h_bench);
	std::vector<double> I_bench = int_lif::currents::pulse_train(2*N_bench, 2*int(startI/h_bench), 2*int(durI/h_bench), aI, 2*int(intI/h_bench), numI);
	std::vector<double> V_bench_tmp(N_bench), t_bench_tmp(N_bench);
	// V_bench_tmp = int_lif::RK4(V0, h_bench, N_bench, I_bench, params);
	for(int i = 0; i < N_bench; i++) {
		if(i*h_bench < startI) V_bench_tmp[i] = -70e-3;
		else if(i*h_bench < startI+durI) V_bench_tmp[i] = -70e-3+10e-3*(1-exp(-(i*h_bench-startI)/1e-2));
		else if(i*h_bench < startI+durI+intI) V_bench_tmp[i] = -70e-3+10e-3*(exp(-(i*h_bench-startI-durI)/1e-2));
		else V_bench_tmp[i] = -70e-3+10e-3*(1-exp(-(i*h_bench-startI-durI-intI)/1e-2));
	}
	
	int percent1 = 0, ii1 = -1;
	for(double h : hvec) {
		int_lif::utils::printPercent(ii1++, percent1, hvec.size(), "");
		int N = int(simtime/h);
		std::vector<double> I = int_lif::currents::pulse_train(N, int(startI/h), int(durI/h), aI, int(intI/h), numI);
		std::vector<double> IRK4 = int_lif::currents::pulse_train(2*N, 2*int(startI/h), 2*int(durI/h), aI, 2*int(intI/h), numI);
		std::vector<double> V1, V2, V3, V4;
		V1 = int_lif::fwdEuler(V0, h, N, I, params);
		V2 = int_lif::bwdEuler(V0, h, N, I, params);
		V3 = int_lif::Heun(V0, h, N, I, params);
		V4 = int_lif::RK4(V0, h, N, IRK4, params);
		std::vector<double> V_bench(N);
		// for(int i = 0; i < N; i++) V_bench[i] = V_bench_tmp[int(h/h_bench*i)];
		for(int i = 0; i < N; i++) {
			if(i*h < startI) V_bench[i] = -70e-3;
			else if(i*h < startI+durI) V_bench[i] = -70e-3+10e-3*(1-exp(-(i*h-startI)/1e-2));
			else if(i*h < startI+durI+intI) V_bench[i] = -70e-3+10e-3*(exp(-(i*h-startI-durI)/1e-2));
			else V_bench[i] = -70e-3+10e-3*(1-exp(-(i*h-startI-durI-intI)/1e-2));
		}
		msevec1.push_back(int_lif::utils::mse(V1,V_bench));
		msevec2.push_back(int_lif::utils::mse(V2,V_bench));
		msevec3.push_back(int_lif::utils::mse(V3,V_bench));
		msevec4.push_back(int_lif::utils::mse(V4,V_bench));
		maevec1.push_back(int_lif::utils::mae(V1,V_bench));
		maevec2.push_back(int_lif::utils::mae(V2,V_bench));
		maevec3.push_back(int_lif::utils::mae(V3,V_bench));
		maevec4.push_back(int_lif::utils::mae(V4,V_bench));
	}
	std::cout << std::endl;
	
	myStyle();
	
	TCanvas* canvas1 = new TCanvas("canvas1", "Canvas1", 600, 400);
	TMultiGraph* multigraph1 = new TMultiGraph();
	TLegend* legend1 = new TLegend(0.15, 0.75, 0.25, 0.85);
	
	addToMultigraph(multigraph1, legend1, hvec, msevec1, hvec.size(), 1, "fwdEuler", "p", "");
	addToMultigraph(multigraph1, legend1, hvec, msevec2, hvec.size(), 2, "bwdEuler", "p", "");
	addToMultigraph(multigraph1, legend1, hvec, msevec3, hvec.size(), 3, "Heun", "p", "");
	addToMultigraph(multigraph1, legend1, hvec, msevec4, hvec.size(), 4, "RK4", "p", "");
	
	canvas1->cd();
	canvas1->SetGrid();
	multigraph1->Draw("AP");
	multigraph1->SetTitle("MSE comparison;h [s];mse [V^2]");
	legend1->Draw();
	canvas1->SaveAs("mse_comparison.pdf");
	multigraph1->GetXaxis()->SetLimits(0.5e-8, 2e-2);
	canvas1->SetLogx();
	canvas1->SaveAs("mse_comparison_lx.pdf");
	multigraph1->SetMaximum(pow(10, ceil(log10(multigraph1->GetHistogram()->GetMaximum()))));
	multigraph1->SetMinimum(pow(10, -15));
	canvas1->SetLogy();
	canvas1->SaveAs("mse_comparison_lxly.pdf");
	
	TCanvas* canvas2 = new TCanvas("canvas2", "Canvas2", 600, 400);
	TMultiGraph* multigraph2 = new TMultiGraph();
	TLegend* legend2 = new TLegend(0.15, 0.75, 0.25, 0.85);
	
	addToMultigraph(multigraph2, legend2, hvec, maevec1, hvec.size(), 1, "fwdEuler", "p", "");
	addToMultigraph(multigraph2, legend2, hvec, maevec2, hvec.size(), 2, "bwdEuler", "p", "");
	addToMultigraph(multigraph2, legend2, hvec, maevec3, hvec.size(), 3, "Heun", "p", "");
	addToMultigraph(multigraph2, legend2, hvec, maevec4, hvec.size(), 4, "RK4", "p", "");
	
	canvas2->cd();
	canvas2->SetGrid();
	multigraph2->Draw("AP");
	multigraph2->SetTitle("MAE comparison;h [s];mae [V]");
	legend2->Draw();
	canvas2->SaveAs("mae_comparison.pdf");
	multigraph2->GetXaxis()->SetLimits(0.5e-8, 2e-2);
	canvas2->SetLogx();
	canvas2->SaveAs("mae_comparison_lx.pdf");
	multigraph2->SetMaximum(pow(10, ceil(log10(multigraph2->GetHistogram()->GetMaximum()))));
	multigraph2->SetMinimum(pow(10, -8));
	canvas2->SetLogy();
	canvas2->SaveAs("mae_comparison_lxly.pdf");
	
	std::vector<double> locerror1, locerror2, locerror3, locerror4;
	int Npoints = 100;
	int percent101 = 0, ii101 = -1;
	for(double h : hvec) {
		int_lif::utils::printPercent(ii101++, percent101, hvec.size(), "");
		locerror1.push_back(int_lif::fwdEulerLocError(V_bench_tmp, h_bench, h, Npoints, I_bench, params));
		locerror2.push_back(int_lif::bwdEulerLocError(V_bench_tmp, h_bench, h, Npoints, I_bench, params));
		locerror3.push_back(int_lif::HeunLocError(V_bench_tmp, h_bench, h, Npoints, I_bench, params));
		locerror4.push_back(int_lif::RK4LocError(V_bench_tmp, h_bench, h, Npoints, I_bench, params));
	}
	std::cout << std::endl;
	
	TCanvas* canvas3 = new TCanvas("canvas3", "Canvas3", 600, 400);
	TMultiGraph* multigraph3 = new TMultiGraph();
	TLegend* legend3 = new TLegend(0.15, 0.75, 0.25, 0.85);
	
	addToMultigraph(multigraph3, legend3, hvec, locerror1, hvec.size(), 1, "fwdEuler", "p", "");
	addToMultigraph(multigraph3, legend3, hvec, locerror2, hvec.size(), 2, "bwdEuler", "p", "");
	addToMultigraph(multigraph3, legend3, hvec, locerror3, hvec.size(), 3, "Heun", "p", "");
	addToMultigraph(multigraph3, legend3, hvec, locerror4, hvec.size(), 4, "RK4", "p", "");
	
	canvas3->cd();
	canvas3->SetGrid();
	multigraph3->Draw("AP");
	multigraph3->SetTitle("locerror comparison;h [s];locerror [V]");
	legend3->Draw();
	canvas3->SaveAs("locerror_comparison.pdf");
	multigraph3->GetXaxis()->SetLimits(0.5e-8, 2e-2);
	canvas3->SetLogx();
	canvas3->SaveAs("locerror_comparison_lx.pdf");
	multigraph3->SetMaximum(pow(10, ceil(log10(multigraph3->GetHistogram()->GetMaximum()))));
	multigraph3->SetMinimum(pow(10, -8));
	canvas3->SetLogy();
	canvas3->SaveAs("locerror_comparison_lxly.pdf");
	
}