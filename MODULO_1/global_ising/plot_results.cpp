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
#include <TNtuple.h>

void plotGraph(const int n, double* x, double* y, double* dx, double* dy, const std::string xname, const std::string yname, const std::string title, const std::string outname) {
	TCanvas* canvas = new TCanvas((outname + "canvas").c_str(), (outname + "canvas").c_str(), 600, 400);
	TGraphErrors* graph = new TGraphErrors(n, x, y, dx, dy);
	canvas->cd();
	canvas->SetGrid();
	graph->SetTitle((title + ";" + xname + ";" + yname).c_str());
	graph->Draw("A*");
	canvas->SaveAs((outname + "_plot.pdf").c_str());
}

int main(int argc, char* argv[]) {
    std::string folder, infilename, outname, xname;
    std::vector<double> x, E, dE, M, dM, C, dC, K, dK, BE, dBE, BM, dBM;
	
	cmdlineParser::CmdlineParser parser;
    parser.addPosParameter<std::string>("folder", &folder, "0220419144203", "[std::string] Working folder name.");
    parser.addPosParameter<std::string>("infile", &infilename, "results.txt", "[std::string] Input file for the results to plot.");
    parser.addOptParameter<std::string>("outname", &outname, "", "[std::string] Output files name prefix.");
    parser.addOptParameter<std::string>("xname", &xname, "", "[std::string] Name for the x axis.");
	if (parser.parseAll(argc, argv) == HELP_RETURN) return 0;
    parser.kickOff(argv[0]);
	
    fs::current_path(fs::current_path() / "measures" / folder);
    int nmeas = textIo::textIn(infilename, '\t', '#', &x, &E, &dE, &M, &dM, &C, &dC, &K, &dK, &BE, &dBE, &BM, &dBM);
	
	plotGraph(nmeas, x.data(), E.data(), nullptr, dE.data(), xname, "#varepsilon", folder, outname + "E");
	plotGraph(nmeas, x.data(), M.data(), nullptr, dM.data(), xname, "|M|", folder, outname + "M");
	plotGraph(nmeas, x.data(), C.data(), nullptr, dC.data(), xname, "d#varepsilon/dT", folder, outname + "C");
	plotGraph(nmeas, x.data(), K.data(), nullptr, dK.data(), xname, "d|M|/dh", folder, outname + "K");
	plotGraph(nmeas, x.data(), BE.data(), nullptr, dBE.data(), xname, "B(#varepsilon)", folder, outname + "BE");
	plotGraph(nmeas, x.data(), BM.data(), nullptr, dBM.data(), xname, "B(|M|)", folder, outname + "BM");
}