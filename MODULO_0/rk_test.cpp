#include <iostream>
#include <vector>
#include <chrono>
#include <math.h>
#include <string>

#include "runge_kutta.h"

#include <TGraph.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TFrame.h>
#include <TLegend.h>

std::vector<double> f(std::vector<double> y, double x, int a, int b) {
	std::vector<double> res;
	double tmp = b / (1. + x * x);
	for(auto yy: y) res.push_back(tmp - a * yy * yy);
	return res;
}

int main() {
	auto rk = rungekutta::RK4<int, int>();
	rk.printTableau();
	
	int p1 = 2, p2 = 1;
	rk.setEquation(f, p1, p2);
	
	double init_x = 0.;
	std::vector<double> init_y = {0.};
	rk.setInit(init_x, init_y);
	
	for(int steps = 10; steps <= 10000; steps*=10) {
		std::cout << steps << " steps : ";
		double h = 0.01;
		rk.clear();
		rk.setInit(init_x, init_y);
		const auto start = std::chrono::steady_clock::now();
		std::vector<std::vector<double>> y = rk.run(steps, h);
		const auto end = std::chrono::steady_clock::now();
		std::cout << (end - start).count() / 1000000. << " ms" << std::endl;
	}
	rk.printToFile();
	rk.clear();
	
	double h = 0.01;
	int nsteps = 1000;
	std::vector<double> xh;
	std::vector<std::vector<double>> yh, dyh;
	rungekutta::runWithRichardsonError(rk, init_x, init_y, nsteps, h, xh, yh, dyh);
	textIo::textOut("rk_richardson.txt", '\t', '#', "x \t y \t dy", nsteps, false, xh, yh, dyh);
	
	double max_x = 5., var_h = 1.;
	int n_h = 4;
	TCanvas* canvas1 = new TCanvas("canvas1", "Canvas1", 600, 400);
	TGraph* graph[n_h];
	TMultiGraph* multigraph = new TMultiGraph();
	TLegend* legend = new TLegend(0.75, 0.75, 0.85, 0.85);
	
	for(int i = 0; i < n_h; i++) {
		int steps = max_x / var_h;
		rk.clear();
		rk.setInit(init_x, init_y);
		auto y_var_h_tmp = rk.run(steps, var_h);
		auto x_var_h = rk.x;
		std::vector<double> y_var_h;
		for(auto& ii: y_var_h_tmp) for(auto& jj: ii) y_var_h.push_back(jj);
		
		std::ostringstream string_stream;
		string_stream << "rkroot" << i << ".txt";
		std::string filename = string_stream.str();
		rk.printToFile(filename);
		
		graph[i] = new TGraph(x_var_h.size(), x_var_h.data(), y_var_h.data());
		graph[i]->SetMarkerColor(i+1);
		multigraph->Add(graph[i]);
		std::ostringstream string_stream2;
		string_stream2 << "h = " << var_h;
		std::string legend_entry = string_stream2.str();
		legend->AddEntry(graph[i], legend_entry.c_str(), "p");
		var_h /= 2.;
	}
	
	canvas1->SetGrid();
	multigraph->Draw("A*");
	legend->Draw();
	canvas1->SaveAs("rkroot.pdf");
}