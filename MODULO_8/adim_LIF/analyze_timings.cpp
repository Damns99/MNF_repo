#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <math.h>
#include <chrono>
#include <ctime>
#include <filesystem>
namespace fs = std::filesystem;
#include <algorithm>

#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TFrame.h>
#include <TLegend.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TH1D.h>

#include "mystyle.h"
#include "integrate_LIF.h"
#include "text_io.h"
#include "cmdline_parser.h"

TGraph* addToMultigraph(TMultiGraph* multigraph, TLegend* legend, std::vector<double>& x, std::vector<double>& y, int n, int color, const char* name, const char* labeltype = "p", const char* drawoption = "") {
	TGraph* graph = new TGraph(n, x.data(), y.data());
	graph->SetMarkerColor(color);
	graph->SetLineColor(color);
	multigraph->Add(graph, drawoption);
	legend->AddEntry(graph, name, labeltype);
	return graph;
}

TF1* fitOneOverH(const char* namefunc, double xleft, double xright, TGraph* graph, int color) {
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

int main(int argc, char* argv[]) {
	myStyle();
	
	int nbins = 10;
	std::vector<double> hvec = int_lif::utils::logspace(-4, 0, 400);
	
	fs::current_path(fs::current_path() / "measures" / "5_timings");
	
	// fwdEuler
	std::vector<double> mean1, stdev1, skew1, kurt1;
	int percent1 = 0, ii1 = -1;
	for(auto& h : hvec) {
		int_lif::utils::printPercent(ii1++, percent1, hvec.size(), "fwdEuler: ");
		std::string filename = std::string("LIF4_all_timings_")+"fwdEuler"+"_h_"+std::to_string(h)+".txt";
		std::vector<double> times;
		textIo::textIn(filename, '\t', '#', &times);
		times.pop_back();
		double mean = int_lif::utils::mean(times);
		double stdev = int_lif::utils::stdev(times);
		double minx = mean - 2*stdev;
		double maxx = mean + 8*stdev;
		//double minx = 1e-8*floor(*std::min_element(times.begin(), times.end())*1e8);
		//double maxx = 1e-8*ceil(*std::max_element(times.begin(), times.end())*1e8);
		TH1D* hist = new TH1D("hist", (filename+";time [s];counts []").c_str(), nbins, minx, maxx);
		for(auto& tt : times) hist->Fill(tt);
		mean1.push_back(hist->GetMean());
		stdev1.push_back(hist->GetStdDev());
		skew1.push_back(hist->GetSkewness());
		kurt1.push_back(hist->GetKurtosis());
		delete hist;
	}
	std::cout << std::endl;
	
	// bwdEuler
	std::vector<double> mean2, stdev2, skew2, kurt2;
	int percent2 = 0, ii2 = -1;
	for(auto& h : hvec) {
		int_lif::utils::printPercent(ii2++, percent2, hvec.size(), "bwdEuler: ");
		std::string filename = std::string("LIF4_all_timings_")+"bwdEuler"+"_h_"+std::to_string(h)+".txt";
		std::vector<double> times;
		textIo::textIn(filename, '\t', '#', &times);
		times.pop_back();
		double mean = int_lif::utils::mean(times);
		double stdev = int_lif::utils::stdev(times);
		double minx = mean - 2*stdev;
		double maxx = mean + 8*stdev;
		//double minx = 1e-8*floor(*std::min_element(times.begin(), times.end())*1e8);
		//double maxx = 1e-8*ceil(*std::max_element(times.begin(), times.end())*1e8);
		TH1D* hist = new TH1D("hist", (filename+";time [s];counts []").c_str(), nbins, minx, maxx);
		for(auto& tt : times) hist->Fill(tt);
		mean2.push_back(hist->GetMean());
		stdev2.push_back(hist->GetStdDev());
		skew2.push_back(hist->GetSkewness());
		kurt2.push_back(hist->GetKurtosis());
		delete hist;
	}
	std::cout << std::endl;
	
	// Heun
	std::vector<double> mean3, stdev3, skew3, kurt3;
	int percent3 = 0, ii3 = -1;
	for(auto& h : hvec) {
		int_lif::utils::printPercent(ii3++, percent3, hvec.size(), "Heun: ");
		std::string filename = std::string("LIF4_all_timings_")+"Heun"+"_h_"+std::to_string(h)+".txt";
		std::vector<double> times;
		textIo::textIn(filename, '\t', '#', &times);
		times.pop_back();
		double mean = int_lif::utils::mean(times);
		double stdev = int_lif::utils::stdev(times);
		double minx = mean - 2*stdev;
		double maxx = mean + 8*stdev;
		//double minx = 1e-8*floor(*std::min_element(times.begin(), times.end())*1e8);
		//double maxx = 1e-8*ceil(*std::max_element(times.begin(), times.end())*1e8);
		TH1D* hist = new TH1D("hist", (filename+";time [s];counts []").c_str(), nbins, minx, maxx);
		for(auto& tt : times) hist->Fill(tt);
		mean3.push_back(hist->GetMean());
		stdev3.push_back(hist->GetStdDev());
		skew3.push_back(hist->GetSkewness());
		kurt3.push_back(hist->GetKurtosis());
		delete hist;
	}
	std::cout << std::endl;
	
	// RK4
	std::vector<double> mean4, stdev4, skew4, kurt4;
	int percent4 = 0, ii4 = -1;
	for(auto& h : hvec) {
		int_lif::utils::printPercent(ii4++, percent4, hvec.size(), "RK4: ");
		std::string filename = std::string("LIF4_all_timings_")+"RK4"+"_h_"+std::to_string(h)+".txt";
		std::vector<double> times;
		textIo::textIn(filename, '\t', '#', &times);
		times.pop_back();
		double mean = int_lif::utils::mean(times);
		double stdev = int_lif::utils::stdev(times);
		double minx = mean - 2*stdev;
		double maxx = mean + 8*stdev;
		//double minx = 1e-8*floor(*std::min_element(times.begin(), times.end())*1e8);
		//double maxx = 1e-8*ceil(*std::max_element(times.begin(), times.end())*1e8);
		TH1D* hist = new TH1D("hist", (filename+";time [s];counts []").c_str(), nbins, minx, maxx);
		for(auto& tt : times) hist->Fill(tt);
		mean4.push_back(hist->GetMean());
		stdev4.push_back(hist->GetStdDev());
		skew4.push_back(hist->GetSkewness());
		kurt4.push_back(hist->GetKurtosis());
		delete hist;
	}
	std::cout << std::endl;
	
	// Heun_naive
	std::vector<double> mean5, stdev5, skew5, kurt5;
	int percent5 = 0, ii5 = -1;
	for(auto& h : hvec) {
		int_lif::utils::printPercent(ii5++, percent5, hvec.size(), "Heun_naive: ");
		std::string filename = std::string("LIF4_all_timings_")+"Heun_naive"+"_h_"+std::to_string(h)+".txt";
		std::vector<double> times;
		textIo::textIn(filename, '\t', '#', &times);
		times.pop_back();
		double mean = int_lif::utils::mean(times);
		double stdev = int_lif::utils::stdev(times);
		double minx = mean - 2*stdev;
		double maxx = mean + 8*stdev;
		//double minx = 1e-8*floor(*std::min_element(times.begin(), times.end())*1e8);
		//double maxx = 1e-8*ceil(*std::max_element(times.begin(), times.end())*1e8);
		TH1D* hist = new TH1D("hist", (filename+";time [s];counts []").c_str(), nbins, minx, maxx);
		for(auto& tt : times) hist->Fill(tt);
		mean5.push_back(hist->GetMean());
		stdev5.push_back(hist->GetStdDev());
		skew5.push_back(hist->GetSkewness());
		kurt5.push_back(hist->GetKurtosis());
		delete hist;
	}
	std::cout << std::endl;
	
	// RK4_naive
	std::vector<double> mean6, stdev6, skew6, kurt6;
	int percent6 = 0, ii6 = -1;
	for(auto& h : hvec) {
		int_lif::utils::printPercent(ii6++, percent6, hvec.size(), "RK4_naive: ");
		std::string filename = std::string("LIF4_all_timings_")+"RK4_naive"+"_h_"+std::to_string(h)+".txt";
		std::vector<double> times;
		textIo::textIn(filename, '\t', '#', &times);
		times.pop_back();
		double mean = int_lif::utils::mean(times);
		double stdev = int_lif::utils::stdev(times);
		double minx = mean - 2*stdev;
		double maxx = mean + 8*stdev;
		//double minx = 1e-8*floor(*std::min_element(times.begin(), times.end())*1e8);
		//double maxx = 1e-8*ceil(*std::max_element(times.begin(), times.end())*1e8);
		TH1D* hist = new TH1D("hist", (filename+";time [s];counts []").c_str(), nbins, minx, maxx);
		for(auto& tt : times) hist->Fill(tt);
		mean6.push_back(hist->GetMean());
		stdev6.push_back(hist->GetStdDev());
		skew6.push_back(hist->GetSkewness());
		kurt6.push_back(hist->GetKurtosis());
		delete hist;
	}
	std::cout << std::endl;
	
	
	TCanvas* canvas1 = new TCanvas("canvas1", "Canvas1", 600, 400);
	TMultiGraph* multigraph1 = new TMultiGraph();
	TLegend* legend1 = new TLegend(0.75, 0.75, 0.85, 0.85);
	
	TCanvas* canvas2 = new TCanvas("canvas2", "Canvas2", 600, 400);
	TMultiGraph* multigraph2 = new TMultiGraph();
	TLegend* legend2 = new TLegend(0.75, 0.75, 0.85, 0.85);
	
	TCanvas* canvas3 = new TCanvas("canvas3", "Canvas3", 600, 400);
	TMultiGraph* multigraph3 = new TMultiGraph();
	TLegend* legend3 = new TLegend(0.75, 0.75, 0.85, 0.85);
	
	TCanvas* canvas4 = new TCanvas("canvas4", "Canvas4", 600, 400);
	TMultiGraph* multigraph4 = new TMultiGraph();
	TLegend* legend4 = new TLegend(0.75, 0.75, 0.85, 0.85);
	
	auto graph1 = addToMultigraph(multigraph1, legend1, hvec, mean1, hvec.size(), 1, "fwdEuler", "p", "");
	auto graph2 = addToMultigraph(multigraph1, legend1, hvec, mean2, hvec.size(), 2, "bwdEuler", "p", "");
	auto graph3 = addToMultigraph(multigraph1, legend1, hvec, mean3, hvec.size(), 3, "Heun", "p", "");
	auto graph4 = addToMultigraph(multigraph1, legend1, hvec, mean4, hvec.size(), 4, "RK4", "p", "");
	auto graph5 = addToMultigraph(multigraph1, legend1, hvec, mean5, hvec.size(), 5, "Heun_naive", "p", "");
	auto graph6 = addToMultigraph(multigraph1, legend1, hvec, mean6, hvec.size(), 6, "RK4_naive", "p", "");
	
	addToMultigraph(multigraph2, legend2, hvec, stdev1, hvec.size(), 1, "fwdEuler", "p", "PL");
	addToMultigraph(multigraph2, legend2, hvec, stdev2, hvec.size(), 2, "bwdEuler", "p", "");
	addToMultigraph(multigraph2, legend2, hvec, stdev3, hvec.size(), 3, "Heun", "p", "");
	addToMultigraph(multigraph2, legend2, hvec, stdev4, hvec.size(), 4, "RK4", "p", "");
	addToMultigraph(multigraph2, legend2, hvec, stdev5, hvec.size(), 5, "Heun_naive", "p", "");
	addToMultigraph(multigraph2, legend2, hvec, stdev6, hvec.size(), 6, "RK4_naive", "p", "PL");
	
	addToMultigraph(multigraph3, legend3, hvec, skew1, hvec.size(), 1, "fwdEuler", "p", "PL");
	addToMultigraph(multigraph3, legend3, hvec, skew2, hvec.size(), 2, "bwdEuler", "p", "");
	addToMultigraph(multigraph3, legend3, hvec, skew3, hvec.size(), 3, "Heun", "p", "");
	addToMultigraph(multigraph3, legend3, hvec, skew4, hvec.size(), 4, "RK4", "p", "");
	addToMultigraph(multigraph3, legend3, hvec, skew5, hvec.size(), 5, "Heun_naive", "p", "");
	addToMultigraph(multigraph3, legend3, hvec, skew6, hvec.size(), 6, "RK4_naive", "p", "PL");
	
	addToMultigraph(multigraph4, legend4, hvec, kurt1, hvec.size(), 1, "fwdEuler", "p", "PL");
	addToMultigraph(multigraph4, legend4, hvec, kurt2, hvec.size(), 2, "bwdEuler", "p", "");
	addToMultigraph(multigraph4, legend4, hvec, kurt3, hvec.size(), 3, "Heun", "p", "");
	addToMultigraph(multigraph4, legend4, hvec, kurt4, hvec.size(), 4, "RK4", "p", "");
	addToMultigraph(multigraph4, legend4, hvec, kurt5, hvec.size(), 5, "Heun_naive", "p", "");
	addToMultigraph(multigraph4, legend4, hvec, kurt6, hvec.size(), 6, "RK4_naive", "p", "PL");
	
	fs::current_path(fs::current_path() / "..");
	
	canvas1->cd();
	canvas1->SetGrid();
	multigraph1->Draw("AP");
	multigraph1->SetTitle("timings comparison;h [s];mean [s]");
	legend1->Draw();
	canvas1->SetLogx();
	canvas1->SetLogy();
	canvas1->SaveAs("5_analyze_timings_mean.pdf");
	
	canvas2->cd();
	canvas2->SetGrid();
	multigraph2->Draw("AP");
	multigraph2->SetTitle("timings comparison;h [s];standard deviation [s]");
	legend2->Draw();
	canvas2->SetLogx();
	canvas2->SetLogy();
	canvas2->SaveAs("5_analyze_timings_stdev.pdf");
	
	canvas3->cd();
	canvas3->SetGrid();
	multigraph3->SetMaximum(33);
	multigraph3->SetMinimum(0);
	multigraph3->Draw("AP");
	multigraph3->SetTitle("timings comparison;h [s];skewness [s]");
	legend3->Draw();
	canvas3->SetLogx();
	canvas3->SaveAs("5_analyze_timings_skew.pdf");
	
	canvas4->cd();
	canvas4->SetGrid();
	multigraph4->SetMaximum(pow(10, 3));
	multigraph4->SetMinimum(pow(10,-2));
	multigraph4->Draw("AP");
	multigraph4->SetTitle("timings comparison;h [s];kurtosis [s]");
	legend4->Draw();
	canvas4->SetLogx();
	canvas4->SetLogy();
	canvas4->SaveAs("5_analyze_timings_kurt.pdf");
	
	
	auto f1 = fitOneOverH("f1", 1e-4, 1e0, graph1, 1);
	auto f2 = fitOneOverH("f2", 1e-4, 1e0, graph2, 2);
	auto f3 = fitOneOverH("f3", 1e-4, 1e0, graph3, 3);
	auto f4 = fitOneOverH("f4", 1e-4, 1e0, graph4, 4);
	auto f5 = fitOneOverH("f5", 1e-4, 1e0, graph5, 5);
	auto f6 = fitOneOverH("f6", 1e-4, 1e0, graph6, 6);
	
	fs::current_path(fs::current_path() / "5_timings");
	
	std::vector<double> mse1(hvec.size(), 0.), mse2(hvec.size(), 0.), mse3(hvec.size(), 0.), mse4(hvec.size(), 0.), mse5(hvec.size(), 0.), mse6(hvec.size(), 0.);
	for(int i = 0; i < hvec.size(); i++) {
		std::string filename1 = std::string("LIF4_all_timings_")+"fwdEuler"+"_h_"+std::to_string(hvec[i])+".txt";
		std::vector<double> alltimevec1;
		textIo::textIn(filename1, '\t', '#', &alltimevec1);
		alltimevec1.pop_back();
		std::string filename2 = std::string("LIF4_all_timings_")+"bwdEuler"+"_h_"+std::to_string(hvec[i])+".txt";
		std::vector<double> alltimevec2;
		textIo::textIn(filename2, '\t', '#', &alltimevec2);
		alltimevec2.pop_back();
		std::string filename3 = std::string("LIF4_all_timings_")+"Heun"+"_h_"+std::to_string(hvec[i])+".txt";
		std::vector<double> alltimevec3;
		textIo::textIn(filename3, '\t', '#', &alltimevec3);
		alltimevec3.pop_back();
		std::string filename4 = std::string("LIF4_all_timings_")+"RK4"+"_h_"+std::to_string(hvec[i])+".txt";
		std::vector<double> alltimevec4;
		textIo::textIn(filename4, '\t', '#', &alltimevec4);
		alltimevec4.pop_back();
		std::string filename5 = std::string("LIF4_all_timings_")+"Heun_naive"+"_h_"+std::to_string(hvec[i])+".txt";
		std::vector<double> alltimevec5;
		textIo::textIn(filename5, '\t', '#', &alltimevec5);
		alltimevec5.pop_back();
		std::string filename6 = std::string("LIF4_all_timings_")+"RK4_naive"+"_h_"+std::to_string(hvec[i])+".txt";
		std::vector<double> alltimevec6;
		textIo::textIn(filename6, '\t', '#', &alltimevec6);
		alltimevec6.pop_back();
		for(int j = 0; j < alltimevec1.size(); j++) mse1[i] += pow(alltimevec1[j] - f1->Eval(hvec[i]), 2) / alltimevec1.size();
		for(int j = 0; j < alltimevec2.size(); j++) mse2[i] += pow(alltimevec2[j] - f2->Eval(hvec[i]), 2) / alltimevec2.size();
		for(int j = 0; j < alltimevec3.size(); j++) mse3[i] += pow(alltimevec3[j] - f3->Eval(hvec[i]), 2) / alltimevec3.size();
		for(int j = 0; j < alltimevec4.size(); j++) mse4[i] += pow(alltimevec4[j] - f4->Eval(hvec[i]), 2) / alltimevec4.size();
		for(int j = 0; j < alltimevec5.size(); j++) mse5[i] += pow(alltimevec5[j] - f5->Eval(hvec[i]), 2) / alltimevec5.size();
		for(int j = 0; j < alltimevec6.size(); j++) mse6[i] += pow(alltimevec6[j] - f6->Eval(hvec[i]), 2) / alltimevec6.size();
	}
	
	fs::current_path(fs::current_path() / "..");
	
	TCanvas* canvasm1 = new TCanvas("canvasm1", "Canvasm1", 600, 400);
	TMultiGraph* multigraphm1 = new TMultiGraph();
	TLegend* legendm1 = new TLegend(0.75, 0.75, 0.85, 0.85);
	
	addToMultigraph(multigraphm1, legendm1, hvec, mse1, hvec.size(), 1, "fwdEuler", "p", "");
	
	canvasm1->cd();
	canvasm1->SetGrid();
	multigraphm1->Draw("AP");
	multigraphm1->SetTitle("mse fit-timings comparison;h [s];MSE(time, fit) [s^2]");
	legendm1->Draw();
	multigraphm1->GetXaxis()->SetLimits(0.5e-4, 2e0);
	multigraphm1->SetMaximum(pow(10, ceil(log10(multigraphm1->GetHistogram()->GetMaximum()))));
	multigraphm1->SetMinimum(pow(10, -17));
	canvasm1->SetLogx();
	canvasm1->SetLogy();
	canvasm1->SaveAs("5_analyze_timings_mse_fit_fwdEuler.pdf");
	
	TCanvas* canvasm2 = new TCanvas("canvasm2", "Canvasm2", 600, 400);
	TMultiGraph* multigraphm2 = new TMultiGraph();
	TLegend* legendm2 = new TLegend(0.75, 0.75, 0.85, 0.85);
	
	addToMultigraph(multigraphm2, legendm2, hvec, mse2, hvec.size(), 2, "bwdEuler", "p", "");
	
	canvasm2->cd();
	canvasm2->SetGrid();
	multigraphm2->Draw("AP");
	multigraphm2->SetTitle("mse fit-timings comparison;h [s];MSE(time, fit) [s^2]");
	legendm2->Draw();
	multigraphm2->GetXaxis()->SetLimits(0.5e-4, 2e0);
	multigraphm2->SetMaximum(pow(10, ceil(log10(multigraphm2->GetHistogram()->GetMaximum()))));
	multigraphm2->SetMinimum(pow(10, -17));
	canvasm2->SetLogx();
	canvasm2->SetLogy();
	canvasm2->SaveAs("5_analyze_timings_mse_fit_bwdEuler.pdf");
	
	TCanvas* canvasm3 = new TCanvas("canvasm3", "Canvasm3", 600, 400);
	TMultiGraph* multigraphm3 = new TMultiGraph();
	TLegend* legendm3 = new TLegend(0.75, 0.75, 0.85, 0.85);
	
	addToMultigraph(multigraphm3, legendm3, hvec, mse3, hvec.size(), 3, "Heun", "p", "");
	
	canvasm3->cd();
	canvasm3->SetGrid();
	multigraphm3->Draw("AP");
	multigraphm3->SetTitle("mse fit-timings comparison;h [s];MSE(time, fit) [s^2]");
	legendm3->Draw();
	multigraphm3->GetXaxis()->SetLimits(0.5e-4, 2e0);
	multigraphm3->SetMaximum(pow(10, ceil(log10(multigraphm3->GetHistogram()->GetMaximum()))));
	multigraphm3->SetMinimum(pow(10, -17));
	canvasm3->SetLogx();
	canvasm3->SetLogy();
	canvasm3->SaveAs("5_analyze_timings_mse_fit_Heun.pdf");
	
	TCanvas* canvasm4 = new TCanvas("canvasm4", "Canvasm4", 600, 400);
	TMultiGraph* multigraphm4 = new TMultiGraph();
	TLegend* legendm4 = new TLegend(0.75, 0.75, 0.85, 0.85);
	
	addToMultigraph(multigraphm4, legendm4, hvec, mse4, hvec.size(), 4, "RK4", "p", "");
	
	canvasm4->cd();
	canvasm4->SetGrid();
	multigraphm4->Draw("AP");
	multigraphm4->SetTitle("mse fit-timings comparison;h [s];MSE(time, fit) [s^2]");
	legendm4->Draw();
	multigraphm4->GetXaxis()->SetLimits(0.5e-4, 2e0);
	multigraphm4->SetMaximum(pow(10, ceil(log10(multigraphm4->GetHistogram()->GetMaximum()))));
	multigraphm4->SetMinimum(pow(10, -17));
	canvasm4->SetLogx();
	canvasm4->SetLogy();
	canvasm4->SaveAs("5_analyze_timings_mse_fit_RK4.pdf");
	
	TCanvas* canvasm5 = new TCanvas("canvasm5", "Canvasm5", 600, 400);
	TMultiGraph* multigraphm5 = new TMultiGraph();
	TLegend* legendm5 = new TLegend(0.75, 0.75, 0.85, 0.85);
	
	addToMultigraph(multigraphm5, legendm5, hvec, mse5, hvec.size(), 5, "Heun_naive", "p", "");
	
	canvasm5->cd();
	canvasm5->SetGrid();
	multigraphm5->Draw("AP");
	multigraphm5->SetTitle("mse fit-timings comparison;h [s];MSE(time, fit) [s^2]");
	legendm5->Draw();
	multigraphm5->GetXaxis()->SetLimits(0.5e-4, 2e0);
	multigraphm5->SetMaximum(pow(10, ceil(log10(multigraphm5->GetHistogram()->GetMaximum()))));
	multigraphm5->SetMinimum(pow(10, -17));
	canvasm5->SetLogx();
	canvasm5->SetLogy();
	canvasm5->SaveAs("5_analyze_timings_mse_fit_Heun_naive.pdf");
	
	TCanvas* canvasm6 = new TCanvas("canvasm6", "Canvasm6", 600, 400);
	TMultiGraph* multigraphm6 = new TMultiGraph();
	TLegend* legendm6 = new TLegend(0.75, 0.75, 0.85, 0.85);
	
	addToMultigraph(multigraphm6, legendm6, hvec, mse6, hvec.size(), 6, "RK4_naive", "p", "");
	
	canvasm6->cd();
	canvasm6->SetGrid();
	multigraphm6->Draw("AP");
	multigraphm6->SetTitle("mse fit-timings comparison;h [s];MSE(time, fit) [s^2]");
	legendm6->Draw();
	multigraphm6->GetXaxis()->SetLimits(0.5e-4, 2e0);
	multigraphm6->SetMaximum(pow(10, ceil(log10(multigraphm6->GetHistogram()->GetMaximum()))));
	multigraphm6->SetMinimum(pow(10, -17));
	canvasm6->SetLogx();
	canvasm6->SetLogy();
	canvasm6->SaveAs("5_analyze_timings_mse_fit_RK4_naive.pdf");
	
}