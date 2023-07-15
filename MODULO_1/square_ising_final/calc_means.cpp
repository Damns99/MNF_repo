#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <fstream>
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

int main(int argc, char* argv[]) {
    std::string folder, measfilename, outname;
    std::vector<double> energy, magnetization, acceptance;
	int iscuda;
	
	cmdlineParser::CmdlineParser parser;
    parser.addPosParameter<std::string>("folder", &folder, "0220419144203", "[std::string] Working folder name.");
    parser.addPosParameter<std::string>("measfile", &measfilename, "metro_ising_meas.txt", "[std::string] Input file for the measures.");
    parser.addOptParameter<std::string>("outname", &outname, "", "[std::string] Output files name prefix.");
    parser.addOptParameter<int>("iscuda", &iscuda, 0, "[int] If != 0 searches in cudamesures/, if 0 in measures/.");
	if (parser.parseAll(argc, argv) == HELP_RETURN) return 0;
    parser.kickOff(argv[0]);
	
	if(iscuda != 0) fs::current_path(fs::current_path() / "cudameasures");
	else fs::current_path(fs::current_path() / "measures");
    fs::current_path(fs::current_path() / folder);
    int nmeas = textIo::textIn(measfilename, '\t', '#', &energy, &magnetization, &acceptance);
	
	int niter = log2(nmeas) - 1;
	double iterations[niter];
	for(int ii = 0; ii < niter; ii++) iterations[ii] = ii;
	
	double e_mean = mean(energy), m_mean = mean(magnetization);
	std::vector<double> e_errors, m_errors;
	e_errors.push_back(stdev(energy) / sqrt(nmeas - 1));
	m_errors.push_back(stdev(magnetization) / sqrt(nmeas - 1));
	
	for(int i = 1; i < niter; i++) {
		nmeas /= 2;
		double e_err_i = 0., m_err_i = 0.;
		double e_mean_i = 0., m_mean_i = 0.;
		for(int j = 0; j < nmeas; j++) {
			energy[j] = (energy[2*j] + energy[2*j+1]) / 2;
			double e_delta1 = energy[j] - e_mean_i;
			e_mean_i += e_delta1 / (j+1);
			double e_delta2 = energy[j] - e_mean_i;
			e_err_i += e_delta1 * e_delta2;
			magnetization[j] = (magnetization[2*j] + magnetization[2*j+1]) / 2;
			double m_delta1 = magnetization[j] - m_mean_i;
			m_mean_i += m_delta1 / (j+1);
			double m_delta2 = magnetization[j] - m_mean_i;
			m_err_i += m_delta1 * m_delta2;
		}
		e_errors.push_back(sqrt(e_err_i / (nmeas * (nmeas - 1.))));
		m_errors.push_back(sqrt(m_err_i / (nmeas * (nmeas - 1.))));
	}
	
	std::string header = " iter \terror(<E>) \terror(<M>) \t <E> = "+ std::to_string(e_mean) + " \t <M> = " + std::to_string(m_mean);
	textIo::textOut(outname + "blocking_errors.txt", '\t', '#', header, niter, false, iterations, e_errors, m_errors);
	
	TCanvas* e_canvas = new TCanvas("e_canvas", "e_canvas", 600, 400);
	TGraph* e_graph = new TGraph(niter, iterations, e_errors.data());
	e_canvas->cd();
	e_canvas->SetGrid();
	e_graph->SetTitle((folder + " <E> = " + std::to_string(e_mean) + ";iterations;error(<E>)").c_str());
	e_graph->Draw("AP");
	e_graph->SetMarkerStyle(kFullCircle);
	e_graph->SetMarkerSize(0.5);
	e_graph->SetMarkerColor(kBlue+1);
	e_canvas->SaveAs((outname + "e_blocking.pdf").c_str());
	
	TCanvas* m_canvas = new TCanvas("m_canvas", "m_canvas", 600, 400);
	TGraph* m_graph = new TGraph(niter, iterations, m_errors.data());
	m_canvas->cd();
	m_canvas->SetGrid();
	m_graph->SetTitle((folder + " <M> = " + std::to_string(m_mean) + ";iterations;error(<M>)").c_str());
	m_graph->Draw("AP");
	m_graph->SetMarkerStyle(kFullCircle);
	m_graph->SetMarkerSize(0.5);
	m_graph->SetMarkerColor(kRed+1);
	m_canvas->SaveAs((outname + "m_blocking.pdf").c_str());
	
	std::ofstream resfile;
	resfile.open("means.txt", std::fstream::out);
	resfile << "<E> = " << e_mean << " +- " << e_errors[niter-2] << std::endl;
	resfile << "<M> = " << m_mean << " +- " << m_errors[niter-2] << std::endl;
	resfile.close();
}