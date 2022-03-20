#include <iostream>
#include <vector>
#include <chrono>
#include <math.h>

#include "runge_kutta.h"

/* Build ROOT with c++17 o c++2a
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
*/

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
		const auto start = std::chrono::steady_clock::now();
		std::vector<std::vector<double>> y = rk.run(steps, h);
		const auto end = std::chrono::steady_clock::now();
		std::cout << (end - start).count() / 1000000. << " ms" << std::endl;
	}
	rk.printToFile();
	rk.clear();
	
	double h = 0.1;
	int nsteps = 1000;
	int order = 4; // qual'è? bho?
	std::vector<double> xh;
	std::vector<std::vector<double>> yh, dyh;
	runWithRichardsonError(rk, init_x, init_y, nsteps, h, order, xh, yh, dyh);
	textIo::textOut("rk_richardson.txt", '\t', '#', "x \t y \t dy", nsteps, false, xh, yh, dyh);
	/*
	TCanvas* canvas1 = new TCanvas("canvas1", "Canvas1");
	TMultiGraph* multigraph = new TMultiGraph();
	
	double max_x = 10.;
	for(double var_h = 1.; var_h >= 0.001; var_h++) {
		int steps = max_x / var_h;
		rk.clear();
		rk.setInit(init_x, init_y);
		auto y_var_h = rk.run(steps, var_h);
		auto x_var_h = rk.x;
		TGraph* graph1 = new TGraph(x_var_h.size(), x_var_h.data(), y_var_h[0].data());
		multigraph->Add(graph1);
	}
	
    canvas1->cd();
	multigraph->Draw();
	*/
}