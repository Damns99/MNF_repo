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

double mean(const std::vector<double>& vec) {
	double s = 0.;
	int length = vec.size();
	for(auto& ii: vec) s += ii / length;
	return s;
}

double stdev(const std::vector<double>& vec) {
	double s1 = 0., s2 = 0.;
	int length = vec.size();
	for(auto& ii: vec) {
		s1 += ii / length;
		s2 += (ii * ii) / length;
	}
	return sqrt(s2 - s1 * s1);
}

double correlation(int k, const std::vector<double>& vec, double vec_mean, double vec_std) {
	double c = 0.;
	int N = vec.size();
	for(int i = 0; i < N - k; i++) c += (vec[i] - vec_mean) * (vec[i+k] - vec_mean);
	c = c / (N - k) / vec_std;
	return c;
}

int main(int argc, char* argv[]) {
    std::string folder, measfilename, outname;
    std::vector<double> energy, magnetization, acceptance;
	int iscuda;
	
	cmdlineParser::CmdlineParser parser;
    parser.addPosParameter<std::string>("folder", &folder, "0220419144203", "[std::string] Working folder name.");
    parser.addPosParameter<std::string>("measfile", &measfilename, "metro_ising_meas", "[std::string] Input file for the measures.");
    parser.addOptParameter<std::string>("outfile", &outname, "", "[std::string] Output files name for the results.");
    parser.addOptParameter<int>("iscuda", &iscuda, 0, "[int] If != 0 searches in cudamesures/, if 0 in measures/.");
	if (parser.parseAll(argc, argv) == HELP_RETURN) return 0;
    parser.kickOff(argv[0]);
	
	if(iscuda != 0) fs::current_path(fs::current_path() / "cudameasures");
	else fs::current_path(fs::current_path() / "measures");
    fs::current_path(fs::current_path() / folder);
    int nmeas = textIo::textIn(measfilename + ".txt", '\t', '#', &energy, &magnetization, &acceptance);
	// for(auto& ii: magnetization) ii = abs(ii);
    
    double energy_mean = mean(energy), energy_std = stdev(energy);
	double magnetization_mean = mean(magnetization), magnetization_std = stdev(magnetization);
	
	std::ofstream outfile;
	outfile.open(outname + "corr.txt", std::fstream::out);
	outfile << "# k \t C_E(k) \t C_M(k)" << std::endl;
	
	std::vector<double> e_corr(nmeas), m_corr(nmeas);
	double time[nmeas];
	for(int k = 0; k < nmeas; k++) {
		time[k] = k;
		e_corr[k] = correlation(k, energy, energy_mean, energy_std);
		m_corr[k] = correlation(k, magnetization, magnetization_mean, magnetization_std);
		outfile << k << '\t' << e_corr[k] << '\t' << m_corr[k] << std::endl;
	}
	
	outfile.close();
	
	TCanvas* e_canvas = new TCanvas("e_canvas", "e_canvas", 600, 400);
	TGraph* e_graph = new TGraph(nmeas, time, e_corr.data());
	e_canvas->cd();
	e_canvas->SetGrid();
	e_graph->SetTitle((folder + ";k;#C_E(k)").c_str());
	e_graph->Draw("AP");
	e_graph->SetMarkerStyle(kFullCircle);
	e_graph->SetMarkerSize(0.2);
	e_graph->SetMarkerColor(kBlue+1);
	e_canvas->SaveAs((outname + "corr_e_plot.pdf").c_str());
	
	TCanvas* m_canvas = new TCanvas("m_canvas", "m_canvas", 600, 400);
	TGraph* m_graph = new TGraph(nmeas, time, m_corr.data());
	m_canvas->cd();
	m_canvas->SetGrid();
	m_graph->SetTitle((folder + ";k;#C_M(k)").c_str());
	m_graph->Draw("AP");
	m_graph->SetMarkerStyle(kFullCircle);
	m_graph->SetMarkerSize(0.2);
	m_graph->SetMarkerColor(kRed+1);
	m_canvas->SaveAs((outname + "corr_m_plot.pdf").c_str());
}