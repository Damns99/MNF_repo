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

void addToMultigraph(TMultiGraph* multigraph, TLegend* legend, std::vector<int>& x, std::vector<double>& y, int n, double simtime, int color, const char* name, const char* labeltype = "p", const char* drawoption = "") {
	std::vector<double> xx(n);
	for(int ii = 0; ii < n; ii++) xx[ii] = simtime/(x[ii]-1);
	TGraph* graph = new TGraph(n, xx.data(), y.data());
	graph->SetMarkerColor(color);
	graph->SetLineColor(color);
	multigraph->Add(graph, drawoption);
	legend->AddEntry(graph, name, labeltype);
}

/* double exact_sol(double x, double offsetz, double periodz, double az) {
	return 2. - ((1. - offsetz) * (1 - exp(-x)) + exp(-x) / (4*M_PI*M_PI/periodz/periodz + 1) * az * (-exp(x)*sin(2*M_PI/periodz*x) + exp(x)*2*M_PI/periodz*cos(2*M_PI/periodz*x) - 2*M_PI/periodz) + 1.*exp(-x)); 
} */
/* double exact_sol(double x, double startz, double durz, double az) {
	if(x<startz) return 1.;
	//else if(x<=startz+durz) return ((az+1)*exp(-(x-startz))*(exp((x-startz))+1) + 1.*exp(-(x-startz)));
	else if(x<startz+durz) return 1. + az - az*exp(-(x-startz));
	else return (exact_sol((startz+durz-1e-10),startz,durz,az)-1)*exp(-(x-startz-durz)) + 1;
} */

double func1(double x, double simtime, double az) {
	if(x<simtime/4) return az*x;
	else if(x<simtime*3/4) return az*(simtime/2-x);
	else return az*(-simtime+x);
}
double exact_sol(double x, double simtime, double az) {
	if(x<simtime/4) return (az-1)*exp(-x)+az*(x-1)+1+1*exp(-x);
	else if(x<simtime*3/4) return exact_sol(simtime/4-1e-10,simtime,az)*exp(simtime/4-x) - (az*simtime/4+az+1)*exp(simtime/4-x) + az*(2*simtime/4-x+1) +1;
	else return exact_sol(simtime*3/4-1e-10,simtime,az)*exp(simtime*3/4-x) + (az*simtime/4+az-1)*exp(simtime*3/4-x) + az*(-simtime+x-1) +1;
}

int main() {
	// Voltages in [Vrest = -70e-3 V], Currents in [g*Vrest = -7e-10 A], Times in [Cm/g = 1e-2 s]
	double params[2] = {
		0.7, // yth
		1.,  // yreset
	};
	double simtime = 1e1;
	// Currents
	double az = -0.1;
	// pulse_train
	double startz = 20, durz = 50, intz = 200;
	int numz = 25;
	// sine_wave
	double periodz = 1., phasez = 0., offsetz = -0.5;
	double x0 = 0.;
	double y0 = 1.;
	
	int nNvec = 400;
	double Nvec0 = 2., Nvec1 = 7.;
	std::vector<int> Nvec; // occhio che sono double ma sono interi
	for(int i = 0; i < nNvec; i++) {
		int newN = floor(pow(10, Nvec0+((Nvec1-Nvec0)*i)/nNvec));
		if(i == 0 || newN != Nvec[Nvec.size()-1]) {
			Nvec.push_back(newN);
			//std::cout << Nvec[Nvec.size()-1] << " ";
		}
	}
	std::vector<double> msevec1, msevec2, msevec3, msevec4;
	std::vector<double> maevec1, maevec2, maevec3, maevec4;
	std::vector<double> finerrvec1, finerrvec2, finerrvec3, finerrvec4;
	std::vector<double> vcfvec1, vcfvec2, vcfvec3, vcfvec4;
	std::vector<double> scfvec1, scfvec2, scfvec3, scfvec4;
	double vcfdelta = 0.05, scfdelta = 0.05;
	
	// Benchmark
	int N_bench = 1e8;
	double h_bench = simtime/(N_bench-1);
	std::vector<double> z_bench = int_lif::currents::sine_wave(2*N_bench-1, h_bench/2, periodz, az, phasez, offsetz);
	std::vector<double> y_bench_tmp, x_bench_tmp, spkt_bench;
	y_bench_tmp = int_lif::RK4(y0, h_bench, N_bench, z_bench, params, &spkt_bench);
	x_bench_tmp = int_lif::utils::linspace(x0, x0 + simtime, N_bench);
	std::vector<double> x_z_bench = int_lif::utils::linspace(x0, x0 + simtime, 2*N_bench-1);
	
	int percent1 = 0, ii1 = -1;
	for(int N : Nvec) {
		int_lif::utils::printPercent(ii1++, percent1, Nvec.size(), "");
		double h = simtime/(N-1);
		//std::vector<double> z = int_lif::currents::sine_wave(N, h, periodz, az, phasez, offsetz);
		//std::vector<double> zRK4 = int_lif::currents::sine_wave(2*N-1, h/2, periodz, az, phasez, offsetz);
	//std::vector<double> z = int_lif::currents::pulse_train(N, int(startz/h), int(durz/h), az, int(intz/h), numz);
	//std::vector<double> zRK4 = int_lif::currents::pulse_train(2*N-1, int(startz/(h/2)), int(durz/(h/2)), az, int(intz/(h/2)), numz);
	std::vector<double> z, zRK4;
	for(int i = 0; i < N; i++) z.push_back(func1(h*i,simtime,az));
	for(int i = 0; i < 2*N-1; i++) zRK4.push_back(func1(h/2*i,simtime,az));
		std::vector<double> y1, y2, y3, y4;
		std::vector<double> spkt1, spkt2, spkt3, spkt4;
		y1 = int_lif::fwdEuler(y0, h, N, z, params, &spkt1);
		y2 = int_lif::bwdEuler(y0, h, N, z, params, &spkt2);
		y3 = int_lif::Heun(y0, h, N, z, params, &spkt3);
		y4 = int_lif::RK4(y0, h, N, zRK4, params, &spkt4);
		std::vector<double> y_bench(N);
		std::vector<double> x = int_lif::utils::linspace(x0, x0 + simtime, N);
		// for(int i = 0; i < N; i++) y_bench[i] = y_bench_tmp[int((N_bench-1.)/(N-1.)*i)];
		// for(int i = 0; i < N; i++) y_bench[i] = exact_sol(x[i], offsetz, periodz, az);
		// for(int i = 0; i < N; i++) y_bench[i] = exact_sol(x[i], startz, durz, az);
		for(int i = 0; i < N; i++) y_bench[i] = exact_sol(x[i], simtime,az);
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
		vcfvec1.push_back(1-int_lif::utils::vcf(y1,y_bench,vcfdelta));
		vcfvec2.push_back(1-int_lif::utils::vcf(y2,y_bench,vcfdelta));
		vcfvec3.push_back(1-int_lif::utils::vcf(y3,y_bench,vcfdelta));
		vcfvec4.push_back(1-int_lif::utils::vcf(y4,y_bench,vcfdelta));
		scfvec1.push_back(int_lif::utils::scf(spkt1,spkt_bench,scfdelta,simtime));
		scfvec2.push_back(int_lif::utils::scf(spkt2,spkt_bench,scfdelta,simtime));
		scfvec3.push_back(int_lif::utils::scf(spkt3,spkt_bench,scfdelta,simtime));
		scfvec4.push_back(int_lif::utils::scf(spkt4,spkt_bench,scfdelta,simtime));
	}
	std::cout << std::endl;
	
	myStyle();
	
	TCanvas* canvas1 = new TCanvas("canvas1", "Canvas1", 600, 400);
	TMultiGraph* multigraph1 = new TMultiGraph();
	TLegend* legend1 = new TLegend(0.15, 0.75, 0.25, 0.85);
	
	addToMultigraph(multigraph1, legend1, Nvec, msevec1, Nvec.size(), simtime, 1, "fwdEuler", "p", "");
	addToMultigraph(multigraph1, legend1, Nvec, msevec2, Nvec.size(), simtime, 2, "bwdEuler", "p", "");
	addToMultigraph(multigraph1, legend1, Nvec, msevec3, Nvec.size(), simtime, 3, "Heun", "p", "");
	addToMultigraph(multigraph1, legend1, Nvec, msevec4, Nvec.size(), simtime, 4, "RK4", "p", "");
	
	canvas1->cd();
	canvas1->SetGrid();
	multigraph1->Draw("AP");
	multigraph1->SetTitle("MSE comparison;h [];mse [V^2]");
	legend1->Draw();
	canvas1->SaveAs("mse_comparison_2.pdf");
	multigraph1->GetXaxis()->SetLimits(0.5e-6, 2e-1);
	canvas1->SetLogx();
	canvas1->SaveAs("mse_comparison_lx_2.pdf");
	multigraph1->SetMaximum(pow(10, ceil(log10(multigraph1->GetHistogram()->GetMaximum()))));
	multigraph1->SetMinimum(pow(10, -28));
	canvas1->SetLogy();
	canvas1->SaveAs("mse_comparison_lxly_2.pdf");
	
	TCanvas* canvas2 = new TCanvas("canvas2", "Canvas2", 600, 400);
	TMultiGraph* multigraph2 = new TMultiGraph();
	TLegend* legend2 = new TLegend(0.15, 0.75, 0.25, 0.85);
	
	addToMultigraph(multigraph2, legend2, Nvec, maevec1, Nvec.size(), simtime, 1, "fwdEuler", "p", "");
	addToMultigraph(multigraph2, legend2, Nvec, maevec2, Nvec.size(), simtime, 2, "bwdEuler", "p", "");
	addToMultigraph(multigraph2, legend2, Nvec, maevec3, Nvec.size(), simtime, 3, "Heun", "p", "");
	addToMultigraph(multigraph2, legend2, Nvec, maevec4, Nvec.size(), simtime, 4, "RK4", "p", "");
	
	canvas2->cd();
	canvas2->SetGrid();
	multigraph2->Draw("AP");
	multigraph2->SetTitle("MAE comparison;h [];mae [V]");
	legend2->Draw();
	canvas2->SaveAs("mae_comparison_2.pdf");
	multigraph2->GetXaxis()->SetLimits(0.5e-6, 2e-1);
	canvas2->SetLogx();
	canvas2->SaveAs("mae_comparison_lx_2.pdf");
	multigraph2->SetMaximum(pow(10, ceil(log10(multigraph2->GetHistogram()->GetMaximum()))));
	multigraph2->SetMinimum(pow(10, -15));
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
	
	addToMultigraph(multigraph4, legend4, Nvec, finerrvec1, Nvec.size(), simtime, 1, "fwdEuler", "p", "");
	addToMultigraph(multigraph4, legend4, Nvec, finerrvec2, Nvec.size(), simtime, 2, "bwdEuler", "p", "");
	addToMultigraph(multigraph4, legend4, Nvec, finerrvec3, Nvec.size(), simtime, 3, "Heun", "p", "");
	addToMultigraph(multigraph4, legend4, Nvec, finerrvec4, Nvec.size(), simtime, 4, "RK4", "p", "");
	
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
	
	TCanvas* canvas5 = new TCanvas("canvas5", "Canvas5", 600, 400);
	TMultiGraph* multigraph5 = new TMultiGraph();
	TLegend* legend5 = new TLegend(0.15, 0.75, 0.25, 0.85);
	
	addToMultigraph(multigraph5, legend5, Nvec, vcfvec1, Nvec.size(), simtime, 1, "fwdEuler", "p", "");
	addToMultigraph(multigraph5, legend5, Nvec, vcfvec2, Nvec.size(), simtime, 2, "bwdEuler", "p", "");
	addToMultigraph(multigraph5, legend5, Nvec, vcfvec3, Nvec.size(), simtime, 3, "Heun", "p", "");
	addToMultigraph(multigraph5, legend5, Nvec, vcfvec4, Nvec.size(), simtime, 4, "RK4", "p", "");
	
	canvas5->cd();
	canvas5->SetGrid();
	multigraph5->Draw("AP");
	multigraph5->SetTitle("VCF comparison;h [];vcf []");
	legend5->Draw();
	canvas5->SaveAs("vcf_comparison_2.pdf");
	multigraph5->GetXaxis()->SetLimits(0.5e-5, 2e0);
	canvas5->SetLogx();
	canvas5->SaveAs("vcf_comparison_lx_2.pdf");
	multigraph5->SetMinimum(pow(10, -5));
	canvas5->SetLogy();
	canvas5->SaveAs("vcf_comparison_lxly_2.pdf");
	
	TCanvas* canvas6 = new TCanvas("canvas6", "Canvas6", 600, 400);
	TMultiGraph* multigraph6 = new TMultiGraph();
	TLegend* legend6 = new TLegend(0.15, 0.75, 0.25, 0.85);
	
	addToMultigraph(multigraph6, legend6, Nvec, scfvec1, Nvec.size(), simtime, 1, "fwdEuler", "p", "");
	addToMultigraph(multigraph6, legend6, Nvec, scfvec2, Nvec.size(), simtime, 2, "bwdEuler", "p", "");
	addToMultigraph(multigraph6, legend6, Nvec, scfvec3, Nvec.size(), simtime, 3, "Heun", "p", "");
	addToMultigraph(multigraph6, legend6, Nvec, scfvec4, Nvec.size(), simtime, 4, "RK4", "p", "");
	
	canvas6->cd();
	canvas6->SetGrid();
	multigraph6->Draw("AP");
	multigraph6->SetTitle("SCF comparison;h [];scf []");
	multigraph6->SetMinimum(-1);
	multigraph6->SetMaximum(1);
	canvas6->SetLogx();
	legend6->Draw();
	canvas6->SaveAs("scf_comparison_2.pdf");
	
	std::vector<double> newNvec;
	for(auto& a : Nvec) newNvec.push_back(a);
	fs::current_path(fs::current_path() / "measures");
	textIo::textOut("LIF5_2_mae.txt", '\t', '#', "h\tfwdEuler\tbwdEuler\tHeun\tRK4", Nvec.size(), false, newNvec, maevec1, maevec2, maevec3, maevec4);
	textIo::textOut("LIF5_2_mse.txt", '\t', '#', "h\tfwdEuler\tbwdEuler\tHeun\tRK4", Nvec.size(), false, newNvec, msevec1, msevec2, msevec3, msevec4);
	// textIo::textOut("LIF5_2_locerror.txt", '\t', '#', "h\tfwdEuler\tbwdEuler\tHeun\tRK4", Nvec.size(), false, Nvec, locerror1, locerror2, locerror3, locerror4);
	textIo::textOut("LIF5_2_finerr.txt", '\t', '#', "h\tfwdEuler\tbwdEuler\tHeun\tRK4", Nvec.size(), false, newNvec, finerrvec1, finerrvec2, finerrvec3, finerrvec4);
}