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

int main(int argc, char* argv[]) {
	myStyle();
	
	std::string filename;
	int nbins;
	
	cmdlineParser::CmdlineParser parser;
    parser.addPosParameter<std::string>("filename", &filename, "LIF4_all_timings_bwdEuler_h_0.000100", "input filename [std::string]");
    parser.addPosParameter<int>("nbins", &nbins, 10, "number of bins [int]");
	
	if(parser.parseAll(argc, argv) == HELP_RETURN) return 0;
    parser.kickOff(argv[0]);
	
	fs::current_path(fs::current_path() / "measures" / "5_timings");
	
	std::vector<double> times;
	textIo::textIn(filename, '\t', '#', &times);
	times.pop_back();
	
	double mean = int_lif::utils::mean(times);
	double stdev = int_lif::utils::stdev(times);
	
	double minx = mean - 2*stdev;
	double maxx = mean + 8*stdev;
	
	//double minx = 1e-8*floor(*std::min_element(times.begin(), times.end())*1e8);
	//double maxx = 1e-8*ceil(*std::max_element(times.begin(), times.end())*1e8);
	
	TCanvas* canvas = new TCanvas("canvas0", "Canvas0", 600, 400);
	TH1D* hist = new TH1D("hist", (filename+";time [s];counts []").c_str(), nbins*3, minx, maxx);
	
	for(auto& tt : times) hist->Fill(tt);
	
	gStyle->SetOptStat(110001111);
	canvas->cd();
	canvas->SetGrid();
	hist->Draw();
	canvas->SetLogy();
	
	fs::current_path(fs::current_path() / "pngs");
	canvas->SaveAs((filename+"_hist.png").c_str());
	
}