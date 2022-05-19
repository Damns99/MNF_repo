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

int main(int argc, char* argv[]) {
    std::string folder, measfilename, outname;
    std::vector<double> energy, magnetization, acceptance;
	
	cmdlineParser::CmdlineParser parser;
    parser.addPosParameter<std::string>("folder", &folder, "0220419144203", "[std::string] Working folder name.");
    parser.addPosParameter<std::string>("measfile", &measfilename, "metro_ising_meas.txt", "[std::string] Input file for the measures.");
    parser.addOptParameter<std::string>("outname", &outname, "", "[std::string] Output files name prefix.");
	if (parser.parseAll(argc, argv) == HELP_RETURN) return 0;
    parser.kickOff(argv[0]);
	
    fs::current_path(fs::current_path() / "measures" / folder);
    int nmeas = textIo::textIn(measfilename, '\t', '#', &energy, &magnetization, &acceptance);
	
	double time[nmeas];
	for(int ii = 0; ii < nmeas; ii++) time[ii] = ii;
	
	TCanvas* e_canvas = new TCanvas("e_canvas", "e_canvas", 600, 400);
	TGraph* e_graph = new TGraph(nmeas, time, energy.data());
	e_canvas->cd();
	e_canvas->SetGrid();
	e_graph->SetTitle((folder + ";time;#varepsilon").c_str());
	e_graph->Draw("A*");
	e_canvas->SaveAs((outname + "e_plot.pdf").c_str());
	
	TCanvas* m_canvas = new TCanvas("m_canvas", "m_canvas", 600, 400);
	TGraph* m_graph = new TGraph(nmeas, time, magnetization.data());
	m_canvas->cd();
	m_canvas->SetGrid();
	m_graph->SetTitle((folder + ";time;M").c_str());
	m_graph->Draw("A*");
	m_canvas->SaveAs((outname + "m_plot.pdf").c_str());
	
	TCanvas* eh_canvas = new TCanvas("eh_canvas", "eh_canvas", 600, 400);
	TH1F* e_hist = new TH1F("e_hist", "e_hist", 100, -2., 2.);
	e_hist->FillN(energy.size(), energy.data(), NULL);
	eh_canvas->cd();
	eh_canvas->SetGrid();
	e_hist->SetTitle((folder + ";#varepsilon;").c_str());
	e_hist->Draw("HIST");
	eh_canvas->SaveAs((outname + "e_hist.pdf").c_str());
	
	TCanvas* mh_canvas = new TCanvas("mh_canvas", "mh_canvas", 600, 400);
	TH1F* m_hist = new TH1F("m_hist", "m_hist", 100, -1., 1.);
	m_hist->FillN(magnetization.size(), magnetization.data(), NULL);
	mh_canvas->cd();
	mh_canvas->SetGrid();
	m_hist->SetTitle((folder + ";M;").c_str());
	m_hist->Draw("HIST");
	mh_canvas->SaveAs((outname + "m_hist.pdf").c_str());
}