#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <fstream>
#include <filesystem>
namespace fs = std::filesystem;

#include "text_io.h"
#include "cmdline_parser.h"

#include <TMultiGraph.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TFrame.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TNtuple.h>

#define N_MAX_LENGTHS 10

int main(int argc, char* argv[]) {
	int nlengths, lengths[N_MAX_LENGTHS];
	std::string infilename, critfilename, outname;
	
	cmdlineParser::CmdlineParser parser;
	parser.addPosParameter<int>("nl", &nlengths, N_MAX_LENGTHS, "[int] Number of lengths to consider.");
    parser.addPosArrayParameter<int>("lengths", lengths, 0, "[int []] Lengths list.", N_MAX_LENGTHS);
	parser.addOptParameter<std::string>("infile", &infilename, "beta_results.txt", "[std::string] Input file name to search.");
	parser.addOptParameter<std::string>("critfile", &critfilename, "critical_fit_results.txt", "[std::string] Input file name to search.");
	parser.addOptParameter<std::string>("outname", &outname, "scalings_K", "[std::string] Output file name prefix for the results.");
	if (parser.parseAll(argc, argv) == HELP_RETURN) return 0;
    parser.kickOff(argv[0]);
	
	assert(nlengths < N_MAX_LENGTHS);
	
	double alpha_nu, dalpha_nu, gamma_nu, dgamma_nu, beta_c, dbeta_c, _nu, d_nu;
	std::ifstream tmpif;
	std::string tmp;
	tmpif.open(critfilename, std::fstream::in);
	tmpif >> tmp >> tmp >> alpha_nu >> tmp >> dalpha_nu;
	tmpif >> tmp >> tmp >> gamma_nu >> tmp >> dgamma_nu;
	tmpif >> tmp >> tmp >> beta_c >> tmp >> dbeta_c;
	tmpif >> tmp >> tmp >> _nu >> tmp >> d_nu;
	tmpif.close();
	
	//_nu = 1.; 
	//gamma_nu = 7./4.;
	//beta_c = 0.66;
	
	TCanvas* K_canvas = new TCanvas("K_canvas", "K_canvas", 600, 400);
	K_canvas->cd();
	K_canvas->SetGrid();
	TGraph* K_graphs[nlengths];
	TMultiGraph* K_multigraph = new TMultiGraph();
	
	fs::current_path(fs::current_path() / "measures");
    auto measures_folder = fs::current_path();
	for(int i = 0; i < nlengths; i++) {
		fs::current_path(fs::current_path() / ("L" + std::to_string(lengths[i])));
		
		std::vector<double> beta, E, dE, M, dM, K, dK, C, dC;
		int nmeas = textIo::textIn(infilename, '\t', '#', &beta, &E, &dE, &M, &dM, &C, &dC, &K, &dK);
		
		fs::current_path(measures_folder);
		
		double x[nmeas], y[nmeas];
		for(int j = 0; j < nmeas; j++) {
			x[j] = (beta[j] - beta_c) * pow(lengths[i], _nu);
			y[j] = K[j] / pow(lengths[i], gamma_nu);
		}
		
		K_graphs[i] = new TGraph(nmeas, x, y);
		K_graphs[i]->SetMarkerStyle(kFullCircle);
		K_graphs[i]->SetMarkerSize(0.25);
		K_graphs[i]->SetMarkerColor(kBlue+i);
		K_multigraph->Add(K_graphs[i]);
	}
	
	K_multigraph->SetTitle(";;");
	K_multigraph->Draw("AP");
	K_canvas->SaveAs((outname + ".pdf").c_str());
}