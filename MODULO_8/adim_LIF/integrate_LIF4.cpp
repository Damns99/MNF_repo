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
#include <TF1.h>
#include <TFitResult.h>

#include "mystyle.h"
#include "integrate_LIF.h"
#include "text_io.h"

TGraphErrors* addToMultigraph(TMultiGraph* multigraph, TLegend* legend, std::vector<double>& x, std::vector<double>& y, std::vector<double>& dy, int n, int color, const char* name, const char* labeltype = "p", const char* drawoption = "") {
	TGraphErrors* graph = new TGraphErrors(n, x.data(), y.data(), nullptr, dy.data());
	graph->SetMarkerColor(color);
	graph->SetLineColor(color);
	multigraph->Add(graph, drawoption);
	legend->AddEntry(graph, name, labeltype);
	return graph;
}

/* TF1* fitOneOverH(const char* namefunc, double xleft, double xright, TGraphErrors* graph, int color) {
	TF1 *f = new TF1(namefunc, "[0]*x^[1]", xleft, xright);
	f->SetParameters(1.,-1.);
	auto fitres = graph->Fit(f, "SQ");
	std::cout << "Fit results for y = p0*x^p1 :" << std::endl;
	std::cout << "p0 = " << fitres->Parameter(0) << " +- " << fitres->ParError(0) << std::endl;
	std::cout << "p1 = " << fitres->Parameter(1) << " +- " << fitres->ParError(1) << std::endl;
	std::cout << "covariance = " << fitres->CovMatrix(0,1) << std::endl;
	std::cout << "chi2/ndof = " << fitres->Chi2() / fitres->Ndf() << std::endl;
	f = graph->GetFunction(namefunc);
	f->SetLineColor(color);
	f->SetLineWidth(1);
	std::cout << std::endl;
	return f;
} */

TF1* fitOneOverH(const char* namefunc, double xleft, double xright, TGraphErrors* graph, int color) {
	TF1* ftmp = new TF1("ftmp", "[0]+[1]*x", log(xleft), log(xright));
	ftmp->SetParameters(1.,-1.);
	int N = graph->GetN();
	auto graphx = graph->GetX(), graphy = graph->GetY();
	std::vector<double> logx, logy;
	for(int i = 0; i < N; i++) {
		if(graphx[i] > 0 && graphy[i] > 0) {
			logx.push_back(log(graphx[i]));
			logy.push_back(log(graphy[i]));
		}
	}
	TGraph* logloggraph = new TGraph(N, logx.data(), logy.data());
	auto fitres = logloggraph->Fit(ftmp, "SQ");
	std::cout << "Fit results for y = e^p0 * x^p1 :" << std::endl;
	std::cout << "p0 = " << fitres->Parameter(0) << " +- " << fitres->ParError(0) << std::endl;
	std::cout << "p1 = " << fitres->Parameter(1) << " +- " << fitres->ParError(1) << std::endl;
	std::cout << "covariance = " << fitres->CovMatrix(0,1) << std::endl;
	std::cout << "chi2/ndof = " << fitres->Chi2() / fitres->Ndf() << std::endl;
	ftmp = graph->GetFunction(namefunc);
	TF1* f = new TF1(namefunc, "exp([0]) * x^[1]", xleft, xright);
	f->SetParameters(fitres->Parameter(0), fitres->Parameter(1));
	f->SetLineColor(color);
	f->SetLineWidth(1);
	std::cout << std::endl;
	return f;
}

int main() {
	// Voltages in [Vrest = -70e-3 V], Currents in [g*Vrest = -7e-10 A], Times in [Cm/g = 1e-2 s]
	double params[2] = {
		0.7, // yth
		1.,  // yreset
	};
	double simtime = 1e2;
	double az = -0.10;
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
	
	// Heun naive
	std::vector<double> timevec5, dtimevec5;
	int percent5 = 0, ii5 = -1;
	for(double h : hvec) {
		int_lif::utils::printPercent(ii5++, percent5, hvec.size(), "Heun naive: ");
		int N = int(simtime/h);
		std::vector<double> z = int_lif::currents::sine_wave(N, int(periodz/h), az, int(phasez/h), offsetz);
		std::vector<double> time(nrep);
		for(int jj = 0; jj < nrep; jj++) {
			std::vector<double> y, x;
			std::chrono::time_point<std::chrono::steady_clock> start, end;
			start = std::chrono::steady_clock::now();
			y = int_lif::Heun_naive(y0, h, N, z, params);
			x = int_lif::utils::linspace(x0, x0 + N*h, N);
			end = std::chrono::steady_clock::now();
			time[jj] = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1.e9; // [s]
		}
		timevec5.push_back(int_lif::utils::mean(time));
		dtimevec5.push_back(int_lif::utils::stdev(time)/sqrt(nrep));
	}
	std::cout << std::endl;
	
	// RK4 naive
	std::vector<double> timevec6, dtimevec6;
	int percent6 = 0, ii6 = -1;
	for(double h : hvec) {
		int_lif::utils::printPercent(ii6++, percent6, hvec.size(), "RK4_naive: ");
		int N = int(simtime/h);
		std::vector<double> zRK4 = int_lif::currents::sine_wave(2*N, 2*int(periodz/h), az, 2*int(phasez/h), offsetz);
		std::vector<double> time(nrep);
		for(int jj = 0; jj < nrep; jj++) {
			std::vector<double> y, x;
			std::chrono::time_point<std::chrono::steady_clock> start, end;
			start = std::chrono::steady_clock::now();
			y = int_lif::RK4_naive(y0, h, N, zRK4, params);
			x = int_lif::utils::linspace(x0, x0 + N*h, N);
			end = std::chrono::steady_clock::now();
			time[jj] = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1.e9; // [s]
		}
		timevec6.push_back(int_lif::utils::mean(time));
		dtimevec6.push_back(int_lif::utils::stdev(time)/sqrt(nrep));
	}
	std::cout << std::endl;
	
	auto graph1 = addToMultigraph(multigraph, legend, hvec, timevec1, dtimevec1, hvec.size(), 1, "fwdEuler", "p", "");
	auto graph2 = addToMultigraph(multigraph, legend, hvec, timevec2, dtimevec2, hvec.size(), 2, "bwdEuler", "p", "");
	auto graph3 = addToMultigraph(multigraph, legend, hvec, timevec3, dtimevec3, hvec.size(), 3, "Heun", "p", "");
	auto graph4 = addToMultigraph(multigraph, legend, hvec, timevec4, dtimevec4, hvec.size(), 4, "RK4", "p", "");
	auto graph5 = addToMultigraph(multigraph, legend, hvec, timevec5, dtimevec5, hvec.size(), 5, "Heun_naive", "p", "");
	auto graph6 = addToMultigraph(multigraph, legend, hvec, timevec6, dtimevec6, hvec.size(), 6, "RK4_naive", "p", "");
	
	auto f1 = fitOneOverH("f1", 1e-4, 1e0, graph1, 1);
	auto f2 = fitOneOverH("f2", 1e-4, 1e0, graph2, 2);
	auto f3 = fitOneOverH("f3", 1e-4, 1e0, graph3, 3);
	auto f4 = fitOneOverH("f4", 1e-4, 1e0, graph4, 4);
	auto f5 = fitOneOverH("f5", 1e-4, 1e0, graph5, 5);
	auto f6 = fitOneOverH("f6", 1e-4, 1e0, graph6, 6);
	
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
	f1->Draw("SAME");
	f2->Draw("SAME");
	f3->Draw("SAME");
	f4->Draw("SAME");
	f5->Draw("SAME");
	f6->Draw("SAME");
	legend->Draw();
	canvas->SetLogy();
	canvas->SaveAs("timings_comparison_lxly.pdf");
	
	/* TCanvas* canvas2 = new TCanvas("canvas2", "Canvas2", 600, 400);
	TMultiGraph* multigraph2 = new TMultiGraph();
	TLegend* legend2 = new TLegend(0.75, 0.75, 0.85, 0.85);
	
	std::vector<double> diff1(hvec.size(), 0.), diff2(hvec.size(), 0.), diff3(hvec.size(), 0.), diff4(hvec.size(), 0.);
	for(int i = 0; i < hvec.size(); i++) {
		if(f1->Eval(hvec[i]) != 0.) diff1[i] = (f1->Eval(hvec[i]) - timevec1[i]) / f1->Eval(hvec[i]);
		if(f2->Eval(hvec[i]) != 0.) diff2[i] = (f2->Eval(hvec[i]) - timevec2[i]) / f2->Eval(hvec[i]);
		if(f3->Eval(hvec[i]) != 0.) diff3[i] = (f3->Eval(hvec[i]) - timevec3[i]) / f3->Eval(hvec[i]);
		if(f4->Eval(hvec[i]) != 0.) diff4[i] = (f4->Eval(hvec[i]) - timevec4[i]) / f4->Eval(hvec[i]);
	}
	
	TGraph* graphdiff1 = new TGraph(hvec.size(), hvec.data(), diff1.data());
	graphdiff1->SetMarkerColor(1);
	graphdiff1->SetLineColor(1);
	multigraph2->Add(graphdiff1, "");
	legend2->AddEntry(graphdiff1, "fwdEuler", "p");
	
	TGraph* graphdiff2 = new TGraph(hvec.size(), hvec.data(), diff2.data());
	graphdiff2->SetMarkerColor(2);
	graphdiff2->SetLineColor(2);
	multigraph2->Add(graphdiff2, "");
	legend2->AddEntry(graphdiff2, "bwdEuler", "p");
	
	TGraph* graphdiff3 = new TGraph(hvec.size(), hvec.data(), diff3.data());
	graphdiff3->SetMarkerColor(3);
	graphdiff3->SetLineColor(3);
	multigraph2->Add(graphdiff3, "");
	legend2->AddEntry(graphdiff3, "Heun", "p");
	
	TGraph* graphdiff4 = new TGraph(hvec.size(), hvec.data(), diff4.data());
	graphdiff4->SetMarkerColor(4);
	graphdiff4->SetLineColor(4);
	multigraph2->Add(graphdiff4, "");
	legend2->AddEntry(graphdiff4, "RK4", "p");
	
	canvas2->cd();
	canvas2->SetGrid();
	multigraph2->Draw("AP");
	multigraph2->SetTitle("diff fit-timings comparison;h [s];(fit - time) / fit []");
	legend2->Draw();
	multigraph2->GetXaxis()->SetLimits(0.5e-4, 2e0);
	canvas2->SetLogx();
	canvas2->SaveAs("timings_diff_fit_comparison.pdf");
	
	fs::current_path(fs::current_path() / "measures");
	textIo::textOut("LIF4_timings.txt", '\t', '#', "h\tfwdEuler\terr\tbwdEuler\terr\tHeun\terr\tRK4\terr", hvec.size(), false, hvec, timevec1, dtimevec1, timevec2, dtimevec2, timevec3, dtimevec3, timevec4, dtimevec4);
	textIo::textOut("LIF4_difftimes.txt", '\t', '#', "h\tfwdEuler\tbwdEuler\tHeun\tRK4", hvec.size(), false, hvec, diff1, diff2, diff3, diff4); */
}