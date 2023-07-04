#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <filesystem>
namespace fs = std::filesystem;
#include <fstream>
#include <algorithm>
#include <stdio.h>

#include <TGraphErrors.h>
#include <TF1.h>
#include <TMatrixDSym.h>
#include <TFitResult.h>

#include "bound_cond_vecs.h"
#include "wave_plots.h"
#include "mystyle.h"
#include "text_io.h"

int main() {
	fs::current_path(fs::current_path() / "measures");
	myStyle();
	
	std::string filename = "HH_max_gaussian_fitres.txt";
	
	std::vector<double> ampl, stdev, m, dm, q, dq, c, dc;
	textIo::textIn(filename, '\t', '#', &ampl, &stdev, &m, &dm, &q, &dq, &c, &dc);
	ampl.pop_back(); stdev.pop_back(); m.pop_back(); dm.pop_back(); q.pop_back(); dq.pop_back(); c.pop_back(); dc.pop_back();
	int N = ampl.size();
	std::vector<double> ampl_unique = ampl, stdev_unique = stdev;
	std::sort(ampl_unique.begin(), ampl_unique.end());
	std::sort(stdev_unique.begin(), stdev_unique.end());
	auto it1 = std::unique(ampl_unique.begin(), ampl_unique.end());
	ampl_unique.resize(std::distance(ampl_unique.begin(), it1));
	auto it2 = std::unique(stdev_unique.begin(), stdev_unique.end());
	stdev_unique.resize(std::distance(stdev_unique.begin(), it2));
	for(auto& aa : ampl_unique) std::cout << aa << std::endl;
	for(auto& aa : stdev_unique) std::cout << aa << std::endl;
	
	{
	TCanvas* canvas = new TCanvas("canvas", "Canvas", 600, 400);
	TMultiGraph* multigraph = new TMultiGraph();
	TLegend* legend = new TLegend(0.75, 0.15, 0.90, 0.35);
	for(int ii = 0; ii < stdev_unique.size(); ii++) {
		std::vector<double> x, y, dy;
		for(int jj = 0; jj < N; jj++) {
			if(stdev[jj] == stdev_unique[ii]) {
				x.push_back(ampl[jj]);
				y.push_back(m[jj]);
				dy.push_back(dm[jj]);
			}
		}
		TGraphErrors* graph = new TGraphErrors(x.size(), x.data(), y.data(), nullptr, nullptr);
		graph->Sort();
		graph->SetMarkerColor(ii+1);
		graph->SetLineColor(ii+1);
		multigraph->Add(graph);
		char leg_entry[30];
		sprintf(leg_entry, "stdev = %.0f cm", stdev_unique[ii]);
		legend->AddEntry(graph, leg_entry, "p");
	}
	multigraph->SetTitle(";ampl [mV];m [s/dam]");
	canvas->cd();
	canvas->SetGrid();
	// canvas->SetLogx();
	// canvas->SetLogy();
	multigraph->Draw("APL");
	legend->Draw();
	
	canvas->SaveAs("HH_test_max_fitres_ampl_m.pdf"); 
	}
	
	{
	TCanvas* canvas = new TCanvas("canvas", "Canvas", 600, 400);
	TMultiGraph* multigraph = new TMultiGraph();
	TLegend* legend = new TLegend(0.75, 0.15, 0.90, 0.35);
	for(int ii = 0; ii < stdev_unique.size(); ii++) {
		std::vector<double> x, y, dy;
		for(int jj = 0; jj < N; jj++) {
			if(stdev[jj] == stdev_unique[ii]) {
				x.push_back(ampl[jj]);
				y.push_back(q[jj]);
				dy.push_back(dq[jj]);
			}
		}
		TGraphErrors* graph = new TGraphErrors(x.size(), x.data(), y.data(), nullptr, dy.data());
		graph->Sort();
		graph->SetMarkerColor(ii+1);
		graph->SetLineColor(ii+1);
		multigraph->Add(graph);
		char leg_entry[30];
		sprintf(leg_entry, "stdev = %.0f cm", stdev_unique[ii]);
		legend->AddEntry(graph, leg_entry, "p");
	}
	multigraph->SetTitle(";ampl [mV];q [ms]");
	canvas->cd();
	canvas->SetGrid();
	// canvas->SetLogx();
	// canvas->SetLogy();
	multigraph->Draw("APL");
	legend->Draw();
	
	canvas->SaveAs("HH_test_max_fitres_ampl_q.pdf"); 
	}
	
	{
	TCanvas* canvas = new TCanvas("canvas", "Canvas", 600, 400);
	TMultiGraph* multigraph = new TMultiGraph();
	TLegend* legend = new TLegend(0.15, 0.65, 0.30, 0.85);
	for(int ii = 0; ii < stdev_unique.size(); ii++) {
		std::vector<double> x, y, dy;
		for(int jj = 0; jj < N; jj++) {
			if(stdev[jj] == stdev_unique[ii]) {
				x.push_back(ampl[jj]);
				y.push_back(c[jj]);
				dy.push_back(dc[jj]);
			}
		}
		TGraphErrors* graph = new TGraphErrors(x.size(), x.data(), y.data(), nullptr, dy.data());
		graph->Sort();
		graph->SetMarkerColor(ii+1);
		graph->SetLineColor(ii+1);
		multigraph->Add(graph);
		char leg_entry[30];
		sprintf(leg_entry, "stdev = %.0f cm", stdev_unique[ii]);
		legend->AddEntry(graph, leg_entry, "p");
	}
	multigraph->SetTitle(";ampl [mV];c [ms]");
	canvas->cd();
	canvas->SetGrid();
	// canvas->SetLogx();
	// canvas->SetLogy();
	multigraph->Draw("APL");
	legend->Draw();
	
	canvas->SaveAs("HH_test_max_fitres_ampl_c.pdf"); 
	}
	
	{
	TCanvas* canvas = new TCanvas("canvas", "Canvas", 600, 400);
	TMultiGraph* multigraph = new TMultiGraph();
	TLegend* legend = new TLegend(0.15, 0.15, 0.30, 0.35);
	for(int ii = 0; ii < stdev_unique.size(); ii++) {
		std::vector<double> x, y, dy;
		for(int jj = 0; jj < N; jj++) {
			if(ampl[jj] == ampl_unique[ii]) {
				x.push_back(stdev[jj]);
				y.push_back(m[jj]);
				dy.push_back(dm[jj]);
			}
		}
		TGraphErrors* graph = new TGraphErrors(x.size(), x.data(), y.data(), nullptr, nullptr);
		graph->Sort();
		graph->SetMarkerColor(ii+1);
		graph->SetLineColor(ii+1);
		multigraph->Add(graph);
		char leg_entry[30];
		sprintf(leg_entry, "ampl = %.0f cm", ampl_unique[ii]);
		legend->AddEntry(graph, leg_entry, "p");
	}
	multigraph->SetTitle(";stdev [mV];m [s/dam]");
	canvas->cd();
	canvas->SetGrid();
	// canvas->SetLogx();
	// canvas->SetLogy();
	multigraph->Draw("APL");
	legend->Draw();
	
	canvas->SaveAs("HH_test_max_fitres_stdev_m.pdf"); 
	}
	
	{
	TCanvas* canvas = new TCanvas("canvas", "Canvas", 600, 400);
	TMultiGraph* multigraph = new TMultiGraph();
	TLegend* legend = new TLegend(0.15, 0.65, 0.30, 0.85);
	for(int ii = 0; ii < stdev_unique.size(); ii++) {
		std::vector<double> x, y, dy;
		for(int jj = 0; jj < N; jj++) {
			if(ampl[jj] == ampl_unique[ii]) {
				x.push_back(stdev[jj]);
				y.push_back(q[jj]);
				dy.push_back(dq[jj]);
			}
		}
		TGraphErrors* graph = new TGraphErrors(x.size(), x.data(), y.data(), nullptr, dy.data());
		graph->Sort();
		graph->SetMarkerColor(ii+1);
		graph->SetLineColor(ii+1);
		multigraph->Add(graph);
		char leg_entry[30];
		sprintf(leg_entry, "ampl = %.0f cm", ampl_unique[ii]);
		legend->AddEntry(graph, leg_entry, "p");
	}
	multigraph->SetTitle(";stdev [mV];q [ms]");
	canvas->cd();
	canvas->SetGrid();
	// canvas->SetLogx();
	// canvas->SetLogy();
	multigraph->Draw("APL");
	legend->Draw();
	
	canvas->SaveAs("HH_test_max_fitres_stdev_q.pdf"); 
	}
	
	{
	TCanvas* canvas = new TCanvas("canvas", "Canvas", 600, 400);
	TMultiGraph* multigraph = new TMultiGraph();
	TLegend* legend = new TLegend(0.15, 0.65, 0.30, 0.85);
	for(int ii = 0; ii < stdev_unique.size(); ii++) {
		std::vector<double> x, y, dy;
		for(int jj = 0; jj < N; jj++) {
			if(ampl[jj] == ampl_unique[ii]) {
				x.push_back(stdev[jj]);
				y.push_back(c[jj]);
				dy.push_back(dc[jj]);
			}
		}
		TGraphErrors* graph = new TGraphErrors(x.size(), x.data(), y.data(), nullptr, dy.data());
		graph->Sort();
		graph->SetMarkerColor(ii+1);
		graph->SetLineColor(ii+1);
		multigraph->Add(graph);
		char leg_entry[30];
		sprintf(leg_entry, "ampl = %.0f cm", ampl_unique[ii]);
		legend->AddEntry(graph, leg_entry, "p");
	}
	multigraph->SetTitle(";stdev [mV];c [ms]");
	canvas->cd();
	canvas->SetGrid();
	// canvas->SetLogx();
	// canvas->SetLogy();
	multigraph->Draw("APL");
	legend->Draw();
	
	canvas->SaveAs("HH_test_max_fitres_stdev_c.pdf"); 
	}
}