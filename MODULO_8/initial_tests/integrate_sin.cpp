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

std::vector<double> sin_fwdEuler(double V0, double h, int N, std::vector<double>& I) {
	std::vector<double> V(N);
	V[0] = V0;
	double new_V = V0;
	for(int n = 1; n < N; n++) {
		new_V = new_V + h * I[n-1];
		V[n] = new_V;
	}
	return V;
}

double sin_fwdEulerLocError(std::vector<double>& V_ref, double h_ref, double h, int Npoints, std::vector<double>& I_ref) {
	int N_ref = V_ref.size();
	int methodSteps = 1;
	int init_ii = ceil(methodSteps * h/h_ref);
	std::vector<double> locerrorvec;
	for(int ii = init_ii; ii < N_ref; ii += int((N_ref-init_ii)/Npoints)) {
		double new_V = V_ref[int(ii-h/h_ref)] + h * I_ref[int(ii-h/h_ref)];
		locerrorvec.push_back(abs(new_V - V_ref[ii]));
	}
	return int_lif::utils::mean(locerrorvec);
}

std::vector<double> sin_bwdEuler(double V0, double h, int N, std::vector<double>& I) {
	std::vector<double> V(N);
	V[0] = V0;
	double new_V = V0;
	for(int n = 1; n < N; n++) {
		new_V = new_V + h * I[n];
		V[n] = new_V;
	}
	return V;
}

double sin_bwdEulerLocError(std::vector<double>& V_ref, double h_ref, double h, int Npoints, std::vector<double>& I_ref) {
	int N_ref = V_ref.size();
	int methodSteps = 1;
	int init_ii = ceil(methodSteps * h/h_ref);
	std::vector<double> locerrorvec;
	for(int ii = init_ii; ii < N_ref; ii += int((N_ref-init_ii)/Npoints)) {
		double new_V = V_ref[int(ii-h/h_ref)] + h * I_ref[int(ii)];
		locerrorvec.push_back(abs(new_V - V_ref[ii]));
	}
	return int_lif::utils::mean(locerrorvec);
}

std::vector<double> sin_Heun(double V0, double h, int N, std::vector<double>& I) {
	std::vector<double> V(N);
	V[0] = V0;
	double new_V = V0;
	for(int n = 1; n < N; n++) {
		new_V = new_V + h/2. * (I[n-1] + I[n]);
		V[n] = new_V;
	}
	return V;
}

double sin_HeunLocError(std::vector<double>& V_ref, double h_ref, double h, int Npoints, std::vector<double>& I_ref) {
	int N_ref = V_ref.size();
	int methodSteps = 2;
	int init_ii = ceil(methodSteps * h/h_ref);
	std::vector<double> locerrorvec;
	for(int ii = init_ii; ii < N_ref; ii += int((N_ref-init_ii)/Npoints)) {
		double new_V = V_ref[int(ii-h/h_ref)] + h/2. * (I_ref[int(ii-h/h_ref)] + I_ref[int(ii)]);
		locerrorvec.push_back(abs(new_V - V_ref[ii]));
	}
	return int_lif::utils::mean(locerrorvec);
}

std::vector<double> sin_RK4(double V0, double h, int N, std::vector<double>& I) {
	std::vector<double> V(N);
	V[0] = V0;
	double new_V = V0;
	for(int n = 1; n < N; n++) {
		new_V = new_V + h/6. * (I[2*n-2] + 4*I[2*n-1] + I[2*n]);
		V[n] = new_V;
	}
	return V;
}

double sin_RK4LocError(std::vector<double>& V_ref, double h_ref, double h, int Npoints, std::vector<double>& I_ref) {
	int N_ref = V_ref.size();
	int methodSteps = 2;
	int init_ii = ceil(methodSteps * h/h_ref);
	std::vector<double> locerrorvec;
	for(int ii = init_ii; ii < N_ref; ii += int((N_ref-init_ii)/Npoints)) {
		double new_V = V_ref[int(ii-h/h_ref)] + h/6. * (I_ref[int(ii-h/h_ref)] + 4*I_ref[int(ii-h/h_ref/2)] + I_ref[int(ii)]);
		locerrorvec.push_back(abs(new_V - V_ref[ii]));
	}
	return int_lif::utils::mean(locerrorvec);
}

// 0th derivative discontinuous
/* double func1(int i, int N, double h) {
	if(i<N/4) return 1.;
	else if(i<N/2) return -1.;
	else if(i<3*N/4) return 1.;
	else return -1.;
}

double func2(int i, int N, double h) {
	if(i<N/4) return i*h;
	else if(i<N/2) return (N/2-i)*h;
	else if(i<3*N/4) return (i-N/2)*h;
	else return (N-i)*h;
} */

// 1st derivative discontinuous
/* double func1(int i, int N, double h) {
	if(i<N/2) return i*h;
	else return (N-i)*h;
}

double func2(int i, int N, double h) {
	if(i<N/2) return pow(i*h, 2)/2;
	else return pow(N*h, 2)/4 - pow((N-i)*h, 2)/2;
} */

// 2nd derivative discontinuous
double func1(int i, int N, double h) {
	if(i<N/2) return pow(i*h, 2)/2;
	else return pow(N*h, 2)/4 - pow((N-i)*h, 2)/2;
}

double func2(int i, int N, double h) {
	if(i<N/2) return pow(i*h, 3)/6;
	else return pow(N*h, 3)/48 + 1./48 * (pow(N*h, 3) - 12*pow(N*h, 2)*i*h + 24*h*N*pow(i*h, 2) - 8*pow(i*h, 3));
}

int main() {
	
	double t0 = 0; // [s]
	double V0 = 0; // [V]
	
	double simtime = 1.; // [s]
	
	std::vector<double> hvec = int_lif::utils::logspace(-7, -2, 400);
	std::vector<double> msevec1, msevec2, msevec3, msevec4;
	std::vector<double> maevec1, maevec2, maevec3, maevec4;
	
	// Benchmark
	double h_bench = 1e-7;
	int N_bench = int(simtime/h_bench);
	std::vector<double> I_bench(N_bench);
	//  for(int i = 0; i < N_bench; i++) I_bench[i] = cos(i*h_bench/simtime);
	//   for(int i = 0; i < N_bench; i++) I_bench[i] = i < N_bench/2 ? 1. : -1.;
	for(int i = 0; i < N_bench; i++) I_bench[i] = func1(i, N_bench, h_bench);
	std::vector<double> V_bench_tmp(N_bench), t_bench_tmp(N_bench);
	// V_bench_tmp = sin_RK4(V0, h_bench, N_bench, I_bench);
	//  for(int i = 0; i < N_bench; i++) V_bench_tmp[i] = sin(i*h_bench/simtime);
	//   for(int i = 0; i < N_bench; i++) V_bench_tmp[i] = i < N_bench/2 ? i*h_bench : (N_bench-i)*h_bench;
	for(int i = 0; i < N_bench; i++) V_bench_tmp[i] = func2(i, N_bench, h_bench);
	t_bench_tmp = int_lif::utils::linspace(t0, t0 + N_bench*h_bench, N_bench);
	
	int percent1 = 0, ii1 = -1;
	for(double h : hvec) {
		int_lif::utils::printPercent(ii1++, percent1, hvec.size(), "");
		int N = int(simtime/h);
		std::vector<double> I(N);
		//  for(int i = 0; i < N; i++) I[i] = cos(i*h/simtime);
		//    for(int i = 0; i < N; i++) I[i] = i < N/2 ? 1. : -1.;
		for(int i = 0; i < N; i++) I[i] = func1(i, N, h);
		std::vector<double> IRK4(2*N);
		//  for(int i = 0; i < 2*N; i++) IRK4[i] = cos(i*h/simtime/2.);
		//    for(int i = 0; i < 2*N; i++) IRK4[i] = i < N ? 1. : -1.;
		for(int i = 0; i < 2*N; i++) IRK4[i] = func1(i, 2*N, h/2);
		std::vector<double> V1, V2, V3, V4;
		V1 = sin_fwdEuler(V0, h, N, I);
		V2 = sin_bwdEuler(V0, h, N, I);
		V3 = sin_Heun(V0, h, N, I);
		V4 = sin_RK4(V0, h, N, IRK4);
		std::vector<double> V_bench(N);
		// for(int i = 0; i < N; i++) V_bench[i] = V_bench_tmp[int(h/h_bench*i)];
		//  for(int i = 0; i < N; i++) V_bench[i] = sin(i*h/simtime);
		//    for(int i = 0; i < N; i++) V_bench[i] = i < N/2 ? i*h : (N-i)*h;
		for(int i = 0; i < N; i++) V_bench[i] = func2(i, N, h);
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
	canvas1->SaveAs("sin_mse_comparison.pdf");
	multigraph1->GetXaxis()->SetLimits(0.5e-8, 2e-2);
	canvas1->SetLogx();
	canvas1->SaveAs("sin_mse_comparison_lx.pdf");
	multigraph1->SetMaximum(pow(10, ceil(log10(multigraph1->GetHistogram()->GetMaximum()))));
	multigraph1->SetMinimum(pow(10, -30));
	canvas1->SetLogy();
	canvas1->SaveAs("sin_mse_comparison_lxly.pdf");
	
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
	canvas2->SaveAs("sin_mae_comparison.pdf");
	multigraph2->GetXaxis()->SetLimits(0.5e-8, 2e-2);
	canvas2->SetLogx();
	canvas2->SaveAs("sin_mae_comparison_lx.pdf");
	multigraph2->SetMaximum(pow(10, ceil(log10(multigraph2->GetHistogram()->GetMaximum()))));
	multigraph2->SetMinimum(pow(10, -15));
	canvas2->SetLogy();
	canvas2->SaveAs("sin_mae_comparison_lxly.pdf");
	
	std::vector<double> locerror1, locerror2, locerror3, locerror4;
	int Npoints = 100;
	int percent101 = 0, ii101 = -1;
	for(double h : hvec) {
		int_lif::utils::printPercent(ii101++, percent101, hvec.size(), "");
		locerror1.push_back(sin_fwdEulerLocError(V_bench_tmp, h_bench, h, Npoints, I_bench));
		locerror2.push_back(sin_bwdEulerLocError(V_bench_tmp, h_bench, h, Npoints, I_bench));
		locerror3.push_back(sin_HeunLocError(V_bench_tmp, h_bench, h, Npoints, I_bench));
		locerror4.push_back(sin_RK4LocError(V_bench_tmp, h_bench, h, Npoints, I_bench));
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
	canvas3->SaveAs("sin_locerror_comparison.pdf");
	multigraph3->GetXaxis()->SetLimits(0.5e-8, 2e-2);
	canvas3->SetLogx();
	canvas3->SaveAs("sin_locerror_comparison_lx.pdf");
	multigraph3->SetMaximum(pow(10, ceil(log10(multigraph3->GetHistogram()->GetMaximum()))));
	multigraph3->SetMinimum(pow(10, -10));
	canvas3->SetLogy();
	canvas3->SaveAs("sin_locerror_comparison_lxly.pdf");
	
	double h = 1e-4;
	int N = int(simtime/h);
	std::vector<double> t = int_lif::utils::linspace(t0, t0 + N*h, N);
	std::vector<double> V1, V2, V3, V4;
	std::vector<double> I(N), IRK4(2*N);
	//  for(int i = 0; i < N; i++) I[i] = cos(i*h/simtime);
	//    for(int i = 0; i < N; i++) I[i] = i < N/2 ? 1. : -1.;
	for(int i = 0; i < N; i++) I[i] = func1(i, N, h);
	for(int i = 0; i < 2*N; i++) IRK4[i] = func1(i, 2*N, h/2);
	
	V1 = sin_fwdEuler(V0, h, N, I);
	V2 = sin_bwdEuler(V0, h, N, I);
	V3 = sin_Heun(V0, h, N, I);
	V4 = sin_RK4(V0, h, N, IRK4);
	
	TCanvas* canvas4 = new TCanvas("canvas4", "Canvas4", 600, 400);
	TMultiGraph* multigraph4 = new TMultiGraph();
	TLegend* legend4 = new TLegend(0.15, 0.75, 0.25, 0.85);
	
	addToMultigraph(multigraph4, legend4, t_bench_tmp, V_bench_tmp, N_bench, kBlack, "benchmark", "l", "L");
	addToMultigraph(multigraph4, legend4, t, V1, N, 1, "fwdEuler", "p", "");
	addToMultigraph(multigraph4, legend4, t, V2, N, 2, "bwdEuler", "p", "");
	addToMultigraph(multigraph4, legend4, t, V3, N, 3, "Heun", "p", "");
	addToMultigraph(multigraph4, legend4, t, V4, N, 4, "RK4", "p", "");
	
	canvas4->cd();
	canvas4->SetGrid();
	multigraph4->Draw("AP");
	multigraph4->SetTitle("sin_integrator_comparison;h [s];locerror [V]");
	legend4->Draw();
	canvas4->SaveAs("sin_integrator_comparison.pdf");
	
}