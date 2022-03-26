#include <iostream>
#include <vector>
#include <chrono>
#include <math.h>
#include <string>

#include "runge_kutta.h"
#include "multistep.h"

#include <TGraph.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TFrame.h>
#include <TLegend.h>

std::vector<double> f(std::vector<double> y, double x) {
	std::vector<double> res;
	double tmp = -2. * exp(-x);
	for(auto yy: y) res.push_back(yy - tmp);
	return res;
}

int main() {
	auto rk4 = rungekutta::RK4<>();
	auto ab4 = multistep::AB4<>();
	auto am3 = multistep::AM3<>();
	
	rk4.setEquation(f);
	ab4.setEquation(f);
	am3.setEquation(f);
	
	double init_x = 0.;
	std::vector<double> init_y = {1.};
	rk4.setInit(init_x, init_y);
	ab4.setInit(init_x, init_y);
	am3.setInit(init_x, init_y);
	
	TCanvas* canvas1 = new TCanvas("canvas1", "Canvas1", 600, 400);
	TGraph* graph[3];
	TMultiGraph* multigraph = new TMultiGraph();
	TLegend* legend = new TLegend(0.75, 0.75, 0.85, 0.85);
	
	double max_x = 2.;
	{
		double h = 0.1;
		int steps = max_x / h;
		const auto start = std::chrono::steady_clock::now();
		auto y_tmp = rk4.run(steps, h);
		const auto end = std::chrono::steady_clock::now();
		auto x = rk4.x;
		std::vector<double> y;
		for(auto& ii: y_tmp) for(auto& jj: ii) y.push_back(jj);
		std::cout << "RK4: " << steps << " steps, h = " << h << " -> " << (end - start).count() / 1000000. << " ms" << std::endl;
		graph[0] = new TGraph(x.size(), x.data(), y.data());
		graph[0]->SetMarkerColor(1);
		multigraph->Add(graph[0]);
		std::ostringstream string_stream2;
		string_stream2 << "RK4: h = " << h;
		std::string legend_entry = string_stream2.str();
		legend->AddEntry(graph[0], legend_entry.c_str(), "p");
	}
	{
		double h = 0.025;
		int steps = max_x / h;
		const auto start = std::chrono::steady_clock::now();
		auto y_tmp = ab4.run(steps, h);
		const auto end = std::chrono::steady_clock::now();
		auto x = ab4.x;
		std::vector<double> y;
		for(auto& ii: y_tmp) for(auto& jj: ii) y.push_back(jj);
		std::cout << "AB4: " << steps << " steps, h = " << h << " -> " << (end - start).count() / 1000000. << " ms" << std::endl;
		graph[1] = new TGraph(x.size(), x.data(), y.data());
		graph[1]->SetMarkerColor(2);
		multigraph->Add(graph[1]);
		std::ostringstream string_stream2;
		string_stream2 << "AB4: h = " << h;
		std::string legend_entry = string_stream2.str();
		legend->AddEntry(graph[1], legend_entry.c_str(), "p");
	}
	{
		double h = 0.05;
		int steps = max_x / h;
		const auto start = std::chrono::steady_clock::now();
		auto y_tmp = am3.run(steps, h);
		const auto end = std::chrono::steady_clock::now();
		auto x = am3.x;
		std::vector<double> y;
		for(auto& ii: y_tmp) for(auto& jj: ii) y.push_back(jj);
		std::cout << "AM3: " << steps << " steps, h = " << h << " -> " << (end - start).count() / 1000000. << " ms" << std::endl;
		graph[2] = new TGraph(x.size(), x.data(), y.data());
		graph[2]->SetMarkerColor(3);
		multigraph->Add(graph[2]);
		std::ostringstream string_stream2;
		string_stream2 << "AM3: h = " << h;
		std::string legend_entry = string_stream2.str();
		legend->AddEntry(graph[2], legend_entry.c_str(), "p");
	}
	
	canvas1->SetGrid();
	multigraph->Draw("A*");
	legend->Draw();
	canvas1->SaveAs("comparison.pdf");
	std::cout << std::endl;
}