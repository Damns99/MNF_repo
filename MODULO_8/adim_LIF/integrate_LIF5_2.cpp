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
#include "text_io.h"

void addToMultigraph(TMultiGraph* multigraph, TLegend* legend, std::vector<double>& x, std::vector<double>& y, int n, int color, const char* name, const char* labeltype = "p", const char* drawoption = "") {
	TGraph* graph = new TGraph(n, x.data(), y.data());
	graph->SetMarkerColor(color);
	graph->SetLineColor(color);
	multigraph->Add(graph, drawoption);
	legend->AddEntry(graph, name, labeltype);
}

void addToMultigraph(TMultiGraph* multigraph, TLegend* legend, std::vector<int>& x, std::vector<double>& y, int n, int color, const char* name, const char* labeltype = "p", const char* drawoption = "") {
	std::vector<double> xx(n);
	for(int ii = 0; ii < n; ii++) xx[ii] = 1.9e0/x[ii];
	TGraph* graph = new TGraph(n, xx.data(), y.data());
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
	double simtime = 1.9e0;
	double az = -1.;
	// double startz = 5, durz = 50, intz = 20;
	// int numz = 25;
	// std::vector<double> z = int_lif::currents::pulse_train(N, int(startz/h), int(durz/h), az, int(intz/h), numz);
	// double periodz = 40, phasez = 3.923, offsetz = -0.20;
	double periodz = 1., phasez = 0., offsetz = 0.;
	double x0 = 0.;
	double y0 = 1.;
	
	int nNvec = 400;
	std::vector<int> Nvec; // occhio che sono double ma sono interi
	for(int i = 0; i < nNvec; i++) {
		int newN = floor(pow(10, 2.+(3.*i)/nNvec));
		if(i == 0 || newN != Nvec[Nvec.size()-1]) {
			Nvec.push_back(newN);
			std::cout << Nvec[Nvec.size()-1] << " ";
		}
	}
	std::vector<double> msevec1, msevec2, msevec3, msevec4;
	std::vector<double> maevec1, maevec2, maevec3, maevec4;
	std::vector<double> finerrvec1, finerrvec2, finerrvec3, finerrvec4;
	
	// Benchmark
	int N_bench = 1e5;
	double h_bench = simtime/N_bench;
	std::vector<double> z_bench = int_lif::currents::sine_wave(2*N_bench, 2*int(periodz/h_bench), az, 2*int(phasez/h_bench), offsetz);
	std::vector<double> y_bench_tmp, x_bench_tmp, spkt_bench;
	y_bench_tmp = int_lif::RK4(y0, h_bench, int(N_bench), z_bench, params, &spkt_bench);
	x_bench_tmp = int_lif::utils::linspace(x0, x0 + simtime, N_bench);
	std::vector<double> x_z_bench = int_lif::utils::linspace(x0, x0 + simtime, 2*N_bench);
	
	int percent1 = 0, ii1 = -1;
	for(int N : Nvec) {
		int_lif::utils::printPercent(ii1++, percent1, Nvec.size(), "");
		double h = simtime/N;
		std::vector<double> z = int_lif::currents::sine_wave(N, int(periodz/h), az, int(phasez/h), offsetz);
		std::vector<double> zRK4 = int_lif::currents::sine_wave(2*N, 2*int(periodz/h), az, 2*int(phasez/h), offsetz);
		std::vector<double> y1, y2, y3, y4;
		y1 = int_lif::fwdEuler(y0, h, N, z, params);
		y2 = int_lif::bwdEuler(y0, h, N, z, params);
		y3 = int_lif::Heun(y0, h, N, z, params);
		y4 = int_lif::RK4(y0, h, N, zRK4, params);
		std::vector<double> y_bench(N);
		for(int i = 0; i < N; i++) y_bench[i] = y_bench_tmp[int((N_bench-1.)/(N-1.)*i)];
		/* for(int i = 0; i < N; i++) {
			if(i*h < startz) y_bench[i] = y0;
			else if(i*h < startz+durz) y_bench[i] = y0+az*(1-exp(-(i*h-startz)));
			else if(i*h < startz+durz+intz) y_bench[i] = y0+az*(exp(-(i*h-startz-durz)));
			else y_bench[i] = y0+az*(1-exp(-(i*h-startz-durz-intz)));
		} */
		msevec1.push_back(int_lif::utils::mse(y1,y_bench));
		msevec2.push_back(int_lif::utils::mse(y2,y_bench));
		msevec3.push_back(int_lif::utils::mse(y3,y_bench));
		msevec4.push_back(int_lif::utils::mse(y4,y_bench));
		maevec1.push_back(int_lif::utils::mae(y1,y_bench));
		maevec2.push_back(int_lif::utils::mae(y2,y_bench));
		maevec3.push_back(int_lif::utils::mae(y3,y_bench));
		maevec4.push_back(int_lif::utils::mae(y4,y_bench));
		finerrvec1.push_back(abs(y1[N-1] - y_bench[N-1]));
		finerrvec2.push_back(abs(y2[N-1] - y_bench[N-1]));
		finerrvec3.push_back(abs(y3[N-1] - y_bench[N-1]));
		finerrvec4.push_back(abs(y4[N-1] - y_bench[N-1]));
	}
	std::cout << std::endl;
	
	myStyle();
	
	TCanvas* canvas1 = new TCanvas("canvas1", "Canvas1", 600, 400);
	TMultiGraph* multigraph1 = new TMultiGraph();
	TLegend* legend1 = new TLegend(0.15, 0.75, 0.25, 0.85);
	
	addToMultigraph(multigraph1, legend1, Nvec, msevec1, Nvec.size(), 1, "fwdEuler", "p", "");
	addToMultigraph(multigraph1, legend1, Nvec, msevec2, Nvec.size(), 2, "bwdEuler", "p", "");
	addToMultigraph(multigraph1, legend1, Nvec, msevec3, Nvec.size(), 3, "Heun", "p", "");
	addToMultigraph(multigraph1, legend1, Nvec, msevec4, Nvec.size(), 4, "RK4", "p", "");
	
	canvas1->cd();
	canvas1->SetGrid();
	multigraph1->Draw("AP");
	multigraph1->SetTitle("MSE comparison;h [];mse [V^2]");
	legend1->Draw();
	canvas1->SaveAs("mse_comparison_2.pdf");
	multigraph1->GetXaxis()->SetLimits(0.5e-5, 2e0);
	canvas1->SetLogx();
	canvas1->SaveAs("mse_comparison_lx_2.pdf");
	multigraph1->SetMaximum(pow(10, ceil(log10(multigraph1->GetHistogram()->GetMaximum()))));
	multigraph1->SetMinimum(pow(10, -12));
	canvas1->SetLogy();
	canvas1->SaveAs("mse_comparison_lxly_2.pdf");
	
	TCanvas* canvas2 = new TCanvas("canvas2", "Canvas2", 600, 400);
	TMultiGraph* multigraph2 = new TMultiGraph();
	TLegend* legend2 = new TLegend(0.15, 0.75, 0.25, 0.85);
	
	addToMultigraph(multigraph2, legend2, Nvec, maevec1, Nvec.size(), 1, "fwdEuler", "p", "");
	addToMultigraph(multigraph2, legend2, Nvec, maevec2, Nvec.size(), 2, "bwdEuler", "p", "");
	addToMultigraph(multigraph2, legend2, Nvec, maevec3, Nvec.size(), 3, "Heun", "p", "");
	addToMultigraph(multigraph2, legend2, Nvec, maevec4, Nvec.size(), 4, "RK4", "p", "");
	
	canvas2->cd();
	canvas2->SetGrid();
	multigraph2->Draw("AP");
	multigraph2->SetTitle("MAE comparison;h [];mae [V]");
	legend2->Draw();
	canvas2->SaveAs("mae_comparison_2.pdf");
	multigraph2->GetXaxis()->SetLimits(0.5e-5, 2e0);
	canvas2->SetLogx();
	canvas2->SaveAs("mae_comparison_lx_2.pdf");
	multigraph2->SetMaximum(pow(10, ceil(log10(multigraph2->GetHistogram()->GetMaximum()))));
	multigraph2->SetMinimum(pow(10, -6));
	canvas2->SetLogy();
	canvas2->SaveAs("mae_comparison_lxly_2.pdf");
	
	/* std::vector<double> locerror1, locerror2, locerror3, locerror4;
	int Npoints = 100;
	int percent101 = 0, ii101 = -1;
	for(double N : Nvec) {
		int_lif::utils::printPercent(ii101++, percent101, Nvec.size(), "");
		locerror1.push_back(int_lif::fwdEulerLocError(y_bench_tmp, h_bench, h, Npoints, z_bench, params));
		locerror2.push_back(int_lif::bwdEulerLocError(y_bench_tmp, h_bench, h, Npoints, z_bench, params));
		locerror3.push_back(int_lif::HeunLocError(y_bench_tmp, h_bench, h, Npoints, z_bench, params));
		locerror4.push_back(int_lif::RK4LocError(y_bench_tmp, h_bench, h, Npoints, z_bench, params));
	}
	std::cout << std::endl;
	
	TCanvas* canvas3 = new TCanvas("canvas3", "Canvas3", 600, 400);
	TMultiGraph* multigraph3 = new TMultiGraph();
	TLegend* legend3 = new TLegend(0.15, 0.75, 0.25, 0.85);
	
	addToMultigraph(multigraph3, legend3, Nvec, locerror1, Nvec.size(), 1, "fwdEuler", "p", "");
	addToMultigraph(multigraph3, legend3, Nvec, locerror2, Nvec.size(), 2, "bwdEuler", "p", "");
	addToMultigraph(multigraph3, legend3, Nvec, locerror3, Nvec.size(), 3, "Heun", "p", "");
	addToMultigraph(multigraph3, legend3, Nvec, locerror4, Nvec.size(), 4, "RK4", "p", "");
	
	canvas3->cd();
	canvas3->SetGrid();
	multigraph3->Draw("AP");
	multigraph3->SetTitle("locerror comparison;N [];locerror [V]");
	legend3->Draw();
	canvas3->SaveAs("locerror_comparison_2.pdf");
	multigraph3->GetXaxis()->SetLimits(0.5e0, 2e5);
	canvas3->SetLogx();
	canvas3->SaveAs("locerror_comparison_lx_2.pdf");
	multigraph3->SetMaximum(pow(10, ceil(log10(multigraph3->GetHistogram()->GetMaximum()))));
	multigraph3->SetMinimum(pow(10, -5));
	canvas3->SetLogy();
	canvas3->SaveAs("locerror_comparison_lxly_2.pdf"); */
	
	TCanvas* canvas4 = new TCanvas("canvas4", "Canvas4", 600, 400);
	TMultiGraph* multigraph4 = new TMultiGraph();
	TLegend* legend4 = new TLegend(0.15, 0.75, 0.25, 0.85);
	
	addToMultigraph(multigraph4, legend4, Nvec, finerrvec1, Nvec.size(), 1, "fwdEuler", "p", "");
	addToMultigraph(multigraph4, legend4, Nvec, finerrvec2, Nvec.size(), 2, "bwdEuler", "p", "");
	addToMultigraph(multigraph4, legend4, Nvec, finerrvec3, Nvec.size(), 3, "Heun", "p", "");
	addToMultigraph(multigraph4, legend4, Nvec, finerrvec4, Nvec.size(), 4, "RK4", "p", "");
	
	canvas4->cd();
	canvas4->SetGrid();
	multigraph4->Draw("AP");
	multigraph4->SetTitle("final error comparison;h [];fin. error [V]");
	legend4->Draw();
	canvas4->SaveAs("finerr_comparison_2.pdf");
	multigraph4->GetXaxis()->SetLimits(0.5e-5, 2e0);
	canvas4->SetLogx();
	canvas4->SaveAs("finerr_comparison_lx_2.pdf");
	multigraph4->SetMaximum(pow(10, ceil(log10(multigraph4->GetHistogram()->GetMaximum()))));
	multigraph4->SetMinimum(pow(10, -6));
	canvas4->SetLogy();
	canvas4->SaveAs("finerr_comparison_lxly_2.pdf");
	
	std::vector<double> newNvec;
	for(auto& a : Nvec) newNvec.push_back(a);
	fs::current_path(fs::current_path() / "measures");
	textIo::textOut("LIF5_2_mae.txt", '\t', '#', "h\tfwdEuler\tbwdEuler\tHeun\tRK4", Nvec.size(), false, newNvec, maevec1, maevec2, maevec3, maevec4);
	textIo::textOut("LIF5_2_mse.txt", '\t', '#', "h\tfwdEuler\tbwdEuler\tHeun\tRK4", Nvec.size(), false, newNvec, msevec1, msevec2, msevec3, msevec4);
	// textIo::textOut("LIF5_2_locerror.txt", '\t', '#', "h\tfwdEuler\tbwdEuler\tHeun\tRK4", Nvec.size(), false, Nvec, locerror1, locerror2, locerror3, locerror4);
	textIo::textOut("LIF5_2_finerr.txt", '\t', '#', "h\tfwdEuler\tbwdEuler\tHeun\tRK4", Nvec.size(), false, newNvec, finerrvec1, finerrvec2, finerrvec3, finerrvec4);
}