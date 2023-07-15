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
#include "lattice.h"

#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TFrame.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TNtuple.h>

#define min(a,b) (a) < (b) ? (a) : (b)

int distsq(int i, int j) {
	int dx = abs((i % length) - (j % length));
	int dy = abs((i / length) - (j / length));
	dx = min(dx, length - dx);
	dy = min(dy, length - dy);
	return dx * dx + dy * dy;
}

int main(int argc, char* argv[]) {
    std::string folder, measfilename, outname;
    std::vector<double> energy, magnetization, acceptance;
	int iscuda;
	
	cmdlineParser::CmdlineParser parser;
    parser.addPosParameter<std::string>("folder", &folder, "0220419144203", "[std::string] Working folder name.");
    parser.addPosParameter<std::string>("measfile", &measfilename, "metro_ising_conf", "[std::string] Input file for the configuration.");
    parser.addOptParameter<std::string>("outfile", &outname, "", "[std::string] Output files name for the results.");
    parser.addOptParameter<int>("iscuda", &iscuda, 0, "[int] If != 0 searches in cudamesures/, if 0 in measures/.");
	if (parser.parseAll(argc, argv) == HELP_RETURN) return 0;
    parser.kickOff(argv[0]);
	
	if(iscuda != 0) fs::current_path(fs::current_path() / "cudameasures");
	else fs::current_path(fs::current_path() / "measures");
    fs::current_path(fs::current_path() / folder);
    load(measfilename + ".txt");
	
	int dsq_max = 2 * (length / 2) * (length / 2);
	std::vector<double> x(dsq_max), y(dsq_max, 0.), dy(dsq_max, 0.), logy(dsq_max, 0.);
	std::vector<int> cy(dsq_max, 0);
	
	for(int i = 0; i < length * length; i++) for(int j = 0; j <= i; j++) {
		int index = distsq(i, j);
		cy[index]++;
		double newvalue = 1. * spin[i] * spin[j];
		double delta1 = newvalue - y[index];
		y[index] += delta1 / cy[index];
		double delta2 = newvalue - y[index];
		dy[index] += delta1 * delta2;
	}
	for(int i = dsq_max - 1; i >= 0; i--) {
		x[i] = sqrt(i);
		if(cy[i] > 1 && dy[i] > 0.) {
			logy[i] = -log(y[i]);
			dy[i] = sqrt(dy[i] / cy[i] / (cy[i] - 1));
		}
		else {
			x.erase(x.begin() + i);
			y.erase(y.begin() + i);
			dy.erase(dy.begin() + i);
			cy.erase(cy.begin() + i);
			logy.erase(logy.begin() + i);
		}
	}
	
	TCanvas* canvas = new TCanvas((outname + "spatial_corr").c_str(), (outname + "spatial_corr").c_str(), 600, 400);
	TGraphErrors* graph = new TGraphErrors(x.size(), x.data(), y.data(), nullptr, dy.data());
	canvas->cd();
	canvas->SetGrid();
	graph->SetTitle((outname + ";" + "|i - j|" + ";" + "<si sj>").c_str());
	graph->Draw("AP");
	graph->SetLineWidth(1);
	graph->SetMarkerStyle(kFullCircle);
	graph->SetMarkerSize(0.2);
	graph->SetMarkerColor(kMagenta+3);
	graph->SetLineColor(kMagenta+3);
	canvas->SaveAs((outname + "spatial_corr.pdf").c_str());
	
	TCanvas* log_canvas = new TCanvas((outname + "spatial_corr_log").c_str(), (outname + "spatial_corr_log").c_str(), 600, 400);
	TGraph* log_graph = new TGraph(x.size(), x.data(), logy.data());
	log_canvas->cd();
	log_canvas->SetGrid();
	log_graph->SetTitle((outname + ";" + "|i - j|" + ";" + "-ln(<si sj>)").c_str());
	log_graph->Draw("AP");
	log_graph->SetMarkerStyle(kFullCircle);
	log_graph->SetMarkerSize(0.2);
	log_graph->SetMarkerColor(kMagenta+3);
	log_canvas->SaveAs((outname + "spatial_corr_log.pdf").c_str());
	
	textIo::textOut(outname + "spatial_corr.txt", '\t', '#', " |i - j| \t<si sj> \terror \t-ln(<si sj>)", x.size(), false, x, y, dy, logy);
}