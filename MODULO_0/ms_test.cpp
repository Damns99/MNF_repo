#include <iostream>
#include <vector>
#include <chrono>
#include <math.h>
#include <string>

#include "multistep.h"

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
	auto ms = multistep::AM2<int, int>();
	ms.print();
	
	int p1 = 2, p2 = 1;
	ms.setEquation(f, p1, p2);
	
	double init_x = 0.;
	std::vector<double> init_y = {0.};
	ms.setInit(init_x, init_y);
	
	for(int steps = 10; steps <= 10000; steps*=10) {
		double h = 0.01;
		ms.clear();
		ms.setInit(init_x, init_y);
		const auto start = std::chrono::steady_clock::now();
		std::vector<std::vector<double>> y = ms.run(steps, h);
		const auto end = std::chrono::steady_clock::now();
		std::cout << steps << " steps : " << (end - start).count() / 1000000. << " ms" << std::endl;
	}
	std::cout << std::endl;
	ms.printToFile();
	ms.clear();
	
	double h = 0.01;
	int nsteps = 1000;
	std::vector<double> xh;
	std::vector<std::vector<double>> yh, y1h, dyh;
	multistep::runWithRichardsonError(ms, init_x, init_y, nsteps, h, xh, yh, y1h, dyh);
	textIo::textOut("ms_richardson.txt", '\t', '#', "x \t y \t y1 \t dy", nsteps, false, xh, yh, y1h, dyh);
	std::cout << std::endl;
	
	double max_x = 5., var_h = 1.;
	int n_h = 4;
	TCanvas* canvas1 = new TCanvas("canvas1", "Canvas1", 600, 400);
	TGraph* graph[n_h];
	TMultiGraph* multigraph = new TMultiGraph();
	TLegend* legend = new TLegend(0.75, 0.75, 0.85, 0.85);
	
	for(int i = 0; i < n_h; i++) {
		int steps = max_x / var_h;
		ms.clear();
		ms.setInit(init_x, init_y);
		auto y_var_h_tmp = ms.run(steps, var_h);
		auto x_var_h = ms.x;
		std::vector<double> y_var_h;
		for(auto& ii: y_var_h_tmp) for(auto& jj: ii) y_var_h.push_back(jj);
		
		std::ostringstream string_stream;
		string_stream << "msroot" << i << ".txt";
		std::string filename = string_stream.str();
		ms.printToFile(filename);
		
		graph[i] = new TGraph(x_var_h.size(), x_var_h.data(), y_var_h.data());
		graph[i]->SetMarkerColor(i+1);
		multigraph->Add(graph[i]);
		std::ostringstream string_stream2;
		string_stream2 << "h = " << var_h;
		std::string legend_entry = string_stream2.str();
		legend->AddEntry(graph[i], legend_entry.c_str(), "p");
		var_h /= 2.;
	}
	std::cout << std::endl;
	
	canvas1->SetGrid();
	multigraph->Draw("A*");
	legend->Draw();
	canvas1->SaveAs("msroot.pdf");
	std::cout << std::endl;
}