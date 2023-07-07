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
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TFrame.h>
#include <TLegend.h>

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

int main() {
	// Voltages in [Vrest = -70e-3 V], Currents in [g*Vrest = -7e-10 A], Times in [Cm/g = 1e-2 s]
	double params[2] = {
		0.7, // yth
		1.,  // yreset
	};
	double simtime = 1.e1;
	int N = 1000;
	double h = simtime/(N-1);
	double x0 = 0.;
	double y0 = 1.;
	int nrep = 100;
	
	myStyle();
	
	TCanvas* canvas = new TCanvas("canvas0", "Canvas0", 600, 400);
	TMultiGraph* multigraph = new TMultiGraph();
	TLegend* legend = new TLegend(0.75, 0.75, 0.85, 0.85);
	
	auto x = int_lif::utils::linspace(x0, x0 + simtime, N);
	auto z1 = int_lif::currents::square_wave(N, N/7, 3., 0, 0.);
	auto z2 = int_lif::currents::sine_wave(N, N/5, 2.5, N/10, -0.5);
	auto z3 = int_lif::currents::constant_pulse(N, N/4, N/3, -2.);
	auto z4 = int_lif::currents::white_noise(N, 0., 0.2, -564564);
	auto z5 = int_lif::currents::ou_noise(N, 0., 0.2, 1., h, -645645);
	
	addToMultigraph(multigraph, legend, x, z1, N, 1, "square_wave", "p", "L");
	addToMultigraph(multigraph, legend, x, z2, N, 2, "sine_wave", "p", "L");
	addToMultigraph(multigraph, legend, x, z3, N, 3, "constant_pulse", "p", "L");
	addToMultigraph(multigraph, legend, x, z4, N, 4, "white_noise", "p", "L");
	addToMultigraph(multigraph, legend, x, z5, N, 5, "ou_noise", "p", "L");
	
	canvas->cd();
	canvas->SetGrid();
	multigraph->Draw("AP");
	multigraph->SetTitle("currents tes;t [Cm/g];I [g*Vrest]");
	legend->Draw();
	canvas->SaveAs("test_currents.pdf");
	
}