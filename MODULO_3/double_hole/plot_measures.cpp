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

#include <TGraph.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TFrame.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TNtuple.h>
#include <TMultiGraph.h>

int main(int argc, char* argv[]) {
    std::string folder, measfilename, outname;
    std::vector<double> particle, acc;
    std::vector<double> obs1, obs2;
	int npart;
	
	cmdlineParser::CmdlineParser parser;
    parser.addPosParameter<std::string>("folder", &folder, "0220520212823", "[std::string] Working folder name.");
    parser.addPosParameter<std::string>("measfile", &measfilename, "double_hole_meas.txt", "[std::string] Input file for the measures.");
    parser.addOptParameter<std::string>("outname", &outname, "", "[std::string] Output files name prefix.");
	if (parser.parseAll(argc, argv) == HELP_RETURN) return 0;
    parser.kickOff(argv[0]);
	
    fs::current_path(fs::current_path() / "measures" / folder);
    int nmeas = textIo::textIn(measfilename, '\t', '#', &particle, &obs1, &obs2, &acc);
	
	std::string tmp;
	std::ifstream tmpif;
	tmpif.open(measfilename, std::fstream::in);
	tmpif >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> npart;
	tmpif.close();
	int pnmeas = nmeas / npart;
	
	double time[pnmeas];
	for(int ii = 0; ii < pnmeas; ii++) time[ii] = ii;
	
	std::vector<double> pobs1[npart], pobs2[npart];
	for(int ii = 0; ii < npart; ii++) {
		pobs1[ii].reserve(pnmeas);
		pobs2[ii].reserve(pnmeas);
	}
	for(int ii = 0; ii < nmeas; ii++) {
		pobs1[int(particle[ii])].push_back(obs1[ii]);
		pobs2[int(particle[ii])].push_back(obs2[ii]);
	}
	
	TCanvas* c1 = new TCanvas("obs1_canvas", "obs1_canvas", 2000, 1000);
	TGraph* graph1[npart];
	TMultiGraph* multigraph1 = new TMultiGraph();
	TLegend* legend1 = new TLegend(0.75, 0.75, 0.85, 0.85);
	for(int i = 0; i < npart; i++) {
		graph1[i] = new TGraph(pnmeas, time, pobs1[i].data());
		graph1[i]->SetMarkerColor(i + 1);
		graph1[i]->SetLineColor(i + 1);
		multigraph1->Add(graph1[i]);
		std::ostringstream string_stream2;
		string_stream2 << "part " << i;
		std::string legend_entry = string_stream2.str();
		legend1->AddEntry(graph1[i], legend_entry.c_str(), "p");
	}
	c1->SetGrid();
	multigraph1->SetTitle((folder + ";time;obs1").c_str());
	multigraph1->Draw("AL*");
	legend1->SetTextSize(0.02);
	legend1->Draw();
	c1->SaveAs((outname + "obs1_plot.pdf").c_str());
	delete c1;
	delete multigraph1;
	delete legend1;
	
	TCanvas* c2 = new TCanvas("obs2_canvas", "obs2_canvas", 2000, 1000);
	TGraph* graph2[npart];
	TMultiGraph* multigraph2 = new TMultiGraph();
	TLegend* legend2 = new TLegend(0.75, 0.75, 0.85, 0.85);
	for(int i = 0; i < npart; i++) {
		graph2[i] = new TGraph(pnmeas, time, pobs2[i].data());
		graph2[i]->SetMarkerColor(i + 1);
		graph2[i]->SetLineColor(i + 1);
		multigraph2->Add(graph2[i]);
		std::ostringstream string_stream2;
		string_stream2 << "part " << i;
		std::string legend_entry = string_stream2.str();
		legend2->AddEntry(graph2[i], legend_entry.c_str(), "p");
	}
	c2->SetGrid();
	multigraph2->SetTitle((folder + ";time;obs2").c_str());
	multigraph2->Draw("AL*");
	legend2->SetTextSize(0.02);
	legend2->Draw();
	c2->SaveAs((outname + "obs2_plot.pdf").c_str());
	delete c2;
	delete multigraph2;
	delete legend2;
}