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
#include <TH1F.h>

#include "mystyle.h"
#include "integrate_LIF.h"
#include "text_io.h"
#include "vec_hist.h"

void addToMultigraph(TMultiGraph* multigraph, TLegend* legend, std::vector<double>& x, std::vector<double>& y, int n, int color, const char* name, const char* labeltype = "p", const char* drawoption = "") {
	TGraph* graph = new TGraph(n, x.data(), y.data());
	graph->SetMarkerColor(color);
	graph->SetLineColor(color);
	multigraph->Add(graph, drawoption);
	legend->AddEntry(graph, name, labeltype);
}

int main() {
	fs::current_path(fs::current_path() / "measures");
	std::string filename11 = "stimulate_response_stim5_neur1.txt", filename21 = "stimulate_response_stim6_neur1.txt", filename12 = "stimulate_response_stim5_neur2.txt", filename22 = "stimulate_response_stim6_neur2.txt";
	std::vector<double> s1r1, s2r1, s1r2, s2r2;
	textIo::textIn(filename11, '\t', '#', &s1r1);
	textIo::textIn(filename21, '\t', '#', &s2r1);
	textIo::textIn(filename12, '\t', '#', &s1r2);
	textIo::textIn(filename22, '\t', '#', &s2r2);
	s1r1.pop_back();
	s2r1.pop_back();
	s1r2.pop_back();
	s2r2.pop_back();
	int N = s1r1.size();
	int nbins = 3;
	double P_Y[2] = {0.5, 0.5};
	double P_X[nbins*nbins], P_XY[nbins*nbins][2];
	double minr1, maxr1, minr2, maxr2;
	minr1 = std::min(*std::min_element(s1r1.begin(), s1r1.end()), *std::min_element(s2r1.begin(), s2r1.end()));
	maxr1 = std::max(*std::max_element(s1r1.begin(), s1r1.end()), *std::max_element(s2r1.begin(), s2r1.end()));
	minr2 = std::min(*std::min_element(s1r2.begin(), s1r2.end()), *std::min_element(s2r2.begin(), s2r2.end()));
	maxr2 = std::max(*std::max_element(s1r2.begin(), s1r2.end()), *std::max_element(s2r2.begin(), s2r2.end()));
	auto hists1r1 = vecHist::makeHist(s1r1, nbins, minr1, maxr1);
	auto hists1r2 = vecHist::makeHist(s1r2, nbins, minr2, maxr2);
	auto hists2r1 = vecHist::makeHist(s2r1, nbins, minr1, maxr1);
	auto hists2r2 = vecHist::makeHist(s2r2, nbins, minr2, maxr2);
	for(int ii = 0; ii < nbins; ii++) for(int jj = 0; jj < nbins; jj++) P_XY[nbins*ii+jj][0] = 1. * (hists1r1[ii] * hists1r2[jj]) / N / N;
	for(int ii = 0; ii < nbins; ii++) for(int jj = 0; jj < nbins; jj++) P_XY[nbins*ii+jj][1] = 1. * (hists2r1[ii] * hists2r2[jj]) / N / N;
	for(int ii = 0; ii < nbins*nbins; ii++) P_X[ii] = 0.5 * (P_XY[ii][0] + P_XY[ii][1]);
	std::cout << minr1 << " " << maxr1 << std::endl;
	for(int ii = 0; ii < nbins; ii++) std::cout << hists1r1[ii] << std::endl;
	for(int ii = 0; ii < nbins; ii++) std::cout << hists2r1[ii] << std::endl;
	std::cout << minr2 << " " << maxr2 << std::endl;
	for(int ii = 0; ii < nbins; ii++) std::cout << hists1r2[ii] << std::endl;
	for(int ii = 0; ii < nbins; ii++) std::cout << hists2r2[ii] << std::endl;
	for(int ii = 0; ii < nbins*nbins; ii++) std::cout << "P_X[" << ii << "] = " << P_X[ii] << std::endl;
	
	double information = 0.;
	for(int ii = 0; ii < nbins*nbins; ii++) for(int jj = 0; jj < 2; jj++) if(P_XY[ii][jj] != 0) information += P_Y[jj] * P_XY[ii][jj] * log(P_XY[ii][jj] / P_X[ii]) / log(2.);
	std::cout << "info = " << information << " bits" << std::endl;
	
	double info1 = 0.;
	for(int ii = 0; ii < nbins; ii++) for(int jj = 0; jj < 2; jj++) {
		if(jj == 0) if(hists1r1[ii] != 0) info1 += P_Y[jj] * hists1r1[ii]/N * log(1. * hists1r1[ii] / (hists1r1[ii]+hists2r1[ii]) * 2) / log(2.);
		if(jj == 1) if(hists2r1[ii] != 0) info1 += P_Y[jj] * hists2r1[ii]/N * log(1. * hists2r1[ii] / (hists1r1[ii]+hists2r1[ii]) * 2) / log(2.);
	}
	std::cout << "info1 = " << info1 << " bits" << std::endl;
	
	double info2 = 0.;
	for(int ii = 0; ii < nbins; ii++) for(int jj = 0; jj < 2; jj++) {
		if(jj == 0) if(hists1r2[ii] != 0) info2 += P_Y[jj] * hists1r2[ii]/N * log(1. * hists1r2[ii] / (hists1r2[ii]+hists2r2[ii]) * 2) / log(2.);
		if(jj == 1) if(hists2r2[ii] != 0) info2 += P_Y[jj] * hists2r2[ii]/N * log(1. * hists2r2[ii] / (hists1r2[ii]+hists2r2[ii]) * 2) / log(2.);
	}
	std::cout << "info2 = " << info2 << " bits" << std::endl;
	
	myStyle();
	
	TCanvas* canvas1 = new TCanvas("canvas1", "Canvas1", 600, 400);
	TMultiGraph* multigraph1 = new TMultiGraph();
	TLegend* legend1 = new TLegend(0.15, 0.75, 0.25, 0.85);
	
	addToMultigraph(multigraph1, legend1, s1r1, s1r2, N, 2, "stim. 1", "p", "");
	addToMultigraph(multigraph1, legend1, s2r1, s2r2, N, 3, "stim. 2", "p", "");
	
	canvas1->cd();
	canvas1->SetGrid();
	multigraph1->Draw("AP");
	multigraph1->SetTitle(";r1 [Cm/g];r2 [Cm/g]");
	legend1->Draw();
	canvas1->SaveAs("information_analysis.pdf");
}