#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <filesystem>
namespace fs = std::filesystem;

#include "text_io.h"
#include "cmdline_parser.h"

#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TFrame.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TNtuple.h>
#include <TMultiGraph.h>

template<typename T>
T max(std::vector<T> v) {
	T res = v[0];
	for(auto& vv: v) if(vv > res) res = vv;
	return res;
}

int main(int argc, char* argv[]) {
    std::string folder, measfilename, outname;
    std::vector<double> particle, tau, corr, dcorr;
	int npart;
	
	cmdlineParser::CmdlineParser parser;
    parser.addPosParameter<std::string>("folder", &folder, "0220520212823", "[std::string] Working folder name.");
    parser.addPosParameter<std::string>("measfile", &measfilename, "two_point_results.txt", "[std::string] Input file for the measures.");
    parser.addOptParameter<std::string>("outname", &outname, "", "[std::string] Output files name prefix.");
	if (parser.parseAll(argc, argv) == HELP_RETURN) return 0;
    parser.kickOff(argv[0]);
	
    fs::current_path(fs::current_path() / "measures" / folder);
    int ntau = textIo::textIn(measfilename, '\t', '#', &particle, &tau, &corr, &dcorr);
	int nparticles = int(max(particle)) + 1;
	
	std::vector<double> pcorr[nparticles], dpcorr[nparticles];
	for(int ii = 0; ii < ntau; ii++) {
		pcorr[int(particle[ii])].push_back(corr[ii]);
		dpcorr[int(particle[ii])].push_back(dcorr[ii]);
	}
	int pntau = ntau / nparticles;
	double ptau[pntau];
	for(int ii = 0; ii < pntau; ii++) ptau[ii] = ii;
	
	TCanvas* c = new TCanvas("two_point_canvas", "two_point_canvas", 2000, 1000);
	TGraphErrors* graph[nparticles];
	TMultiGraph* multigraph = new TMultiGraph();
	TLegend* legend = new TLegend(0.75, 0.75, 0.85, 0.85);
	for(int i = 0; i < nparticles; i++) {
		graph[i] = new TGraphErrors(pntau, ptau, pcorr[i].data(), nullptr, dpcorr[i].data());
		graph[i]->SetMarkerColor(i + 1);
		graph[i]->SetLineColor(i + 1);
		multigraph->Add(graph[i]);
		std::ostringstream string_stream2;
		string_stream2 << "particle " << i;
		std::string legend_entry = string_stream2.str();
		legend->AddEntry(graph[i], legend_entry.c_str(), "p");
	}
	c->SetGrid();
	multigraph->SetTitle((folder + ";#tau;conn_corr").c_str());
	multigraph->Draw("A*");
	c->SaveAs((outname + "two_point_plot.pdf").c_str());
	delete c;
	delete multigraph;
	delete legend;
}