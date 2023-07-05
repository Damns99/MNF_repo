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
#include "cmdline_parser.h"

int main(int argc, char* argv[]) {
	// Voltages in [Vrest = -70e-3 V], Currents in [g*Vrest = -7e-10 A], Times in [Cm/g = 1e-2 s]
	double params[2] = {
		0.7, // yth
		1.,  // yreset
	};
	double simtime = 1e1;
	double az = -1.;
	// double startz = 5, durz = 50, intz = 20;
	// int numz = 25;
	// std::vector<double> z = int_lif::currents::pulse_train(N, int(startz/h), int(durz/h), az, int(intz/h), numz);
	// double periodz = 40, phasez = 3.923, offsetz = -0.20;
	double periodz = 1, phasez = 0., offsetz = -0.5;
	double x0 = 0.;
	double y0 = 1.;
	double ou_sigma = 0.5, ou_tau = 0.1;
	std::string filename;
	
	cmdlineParser::CmdlineParser parser;
    parser.addOptParameter<double>("periodz", &periodz, 1., "input period [double]");
    parser.addOptParameter<double>("offsetz", &offsetz, 0.5, "input offset [double]");
    parser.addOptParameter<double>("ou_sigma", &ou_sigma, 0.5, "OU noise sigma [double]");
    parser.addOptParameter<double>("ou_tau", &ou_tau, 0.1, "OU noise tau [double]");
    parser.addPosParameter<std::string>("filename", &filename, "stimulate_response_stim1_neur2.txt", "output file name [string]");
	if(parser.parseAll(argc, argv) == HELP_RETURN) return 0;
    parser.kickOff(argv[0]);
	offsetz = -offsetz;
	
	int N = 1e4;
	double h = simtime/(N-1);
	std::vector<double> z0 = int_lif::currents::sine_wave(N, h, periodz, az, phasez, offsetz);
	
	std::vector<double> z1 = int_lif::currents::ou_noise(N, 0., ou_sigma, ou_tau, h, -long(4243332235*ou_sigma*ou_tau+8763));
	std::vector<double> z(N);
	for(int ii = 0; ii < N; ii++) z[ii] = z0[ii] + z1[ii];
	
	std::vector<double> x, y, spkt;
	x = int_lif::utils::linspace(x0, x0 + simtime, N);
	y = int_lif::Heun(y0, h, N, z, params, &spkt);
	
	myStyle();
	
	TCanvas* canvas1 = new TCanvas("canvas1", "Canvas1", 600, 400);
	TGraph* graph1 = new TGraph(N, x.data(), y.data());
	canvas1->cd();
	canvas1->SetGrid();
	graph1->Draw("AP");
	graph1->SetTitle("example;x [Cm/g];y [Vrest]");
	canvas1->SaveAs("stimulate_response_example.pdf");
	
	double mean_isi = 0.;
	for(int ii = 1; ii < spkt.size(); ii++) mean_isi += spkt[ii] - spkt[ii-1];
	mean_isi /= (spkt.size()-1);
	std::cout << "Mean ISI = " << mean_isi << " Cm/g" << std::endl;
	
	
	int nrep = 100;
	std::vector<double> mean_isi_vec;
	for(int jj = 0; jj < nrep; jj++) {
		z1 = int_lif::currents::ou_noise(N, 0., ou_sigma, ou_tau, h, -4233235*jj);
		for(int ii = 0; ii < N; ii++) z[ii] = z0[ii] + z1[ii];
		spkt.clear();
		x = int_lif::utils::linspace(x0, x0 + simtime, N);
		y = int_lif::Heun(y0, h, N, z, params, &spkt);
		mean_isi = 0.;
		for(int ii = 1; ii < spkt.size(); ii++) mean_isi += spkt[ii] - spkt[ii-1];
		mean_isi /= (spkt.size()-1);
		mean_isi_vec.push_back(mean_isi);
	}
	
	fs::current_path(fs::current_path() / "measures");
	textIo::textOut(filename, '\t', '#', "mean_isi", mean_isi_vec.size(), false, mean_isi_vec);
}