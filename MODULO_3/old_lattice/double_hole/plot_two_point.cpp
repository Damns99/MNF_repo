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
#include <TF1.h>
#include <TMatrixDSym.h>
#include <TFitResult.h>

template<typename T>
T max(std::vector<T> v) {
	T res = v[0];
	for(auto& vv: v) if(vv > res) res = vv;
	return res;
}

int main(int argc, char* argv[]) {
    std::string folder, measfilename, outname;
    std::vector<double> particle, tau, corr, dcorr, corr2, dcorr2;
	int npart;
	double beta, lambda;
	
	cmdlineParser::CmdlineParser parser;
    parser.addPosParameter<std::string>("folder", &folder, "0220520212823", "[std::string] Working folder name.");
    parser.addPosParameter<std::string>("measfile", &measfilename, "two_point_results.txt", "[std::string] Input file for the measures.");
    parser.addOptParameter<std::string>("outname", &outname, "", "[std::string] Output files name prefix.");
    parser.addOptParameter<double>("beta", &beta, 3., "[double] One over temperature of the system.");
	parser.addOptParameter<double>("lambda", &lambda, 1., "[double] Height of the wall.");
	if (parser.parseAll(argc, argv) == HELP_RETURN) return 0;
    parser.kickOff(argv[0]);
	
    fs::current_path(fs::current_path() / "measures" / folder);
    int ntau = textIo::textIn(measfilename, '\t', '#', &particle, &tau, &corr, &dcorr, &corr2, &dcorr2);
	int nparticles = int(max(particle)) + 1;
	
	std::vector<double> pcorr[nparticles], dpcorr[nparticles], pcorr2[nparticles], dpcorr2[nparticles];
	for(int ii = 0; ii < ntau; ii++) {
		pcorr[int(particle[ii])].push_back(corr[ii]);
		dpcorr[int(particle[ii])].push_back(dcorr[ii]);
		pcorr2[int(particle[ii])].push_back(corr2[ii]);
		dpcorr2[int(particle[ii])].push_back(dcorr2[ii]);
	}
	int pntau = ntau / nparticles;
	double ptau[pntau];
	for(int ii = 0; ii < pntau; ii++) ptau[ii] = ii;
	
	{
		TCanvas* c = new TCanvas("two_point_canvas", "two_point_canvas", 2000, 1000);
		TGraphErrors* graph[nparticles];
		TMultiGraph* multigraph = new TMultiGraph();
		TLegend* legend = new TLegend(0.75, 0.75, 0.9, 0.85);
		
		for(int i = 0; i < nparticles; i++) {
			graph[i] = new TGraphErrors(pntau, ptau, pcorr[i].data(), nullptr, dpcorr[i].data());
			graph[i]->SetMarkerColorAlpha(i + 1, 0.65);
			graph[i]->SetLineColorAlpha(i + 1, 0.5);
			
			TF1* func = new TF1("func", "[0] + [1] * exp(- [2] * x)", 0., pntau - 1.);
			func->SetParameters(-0.05, 0.5, beta / pntau);
			func->SetParNames("cinf","k","1/csi");
			func->SetNpx(pntau * 10);
			func->SetLineColorAlpha(i + 1, 1.);
			func->SetLineStyle(kDashed);
			
			TFitResultPtr result = graph[i]->Fit("func", "S+", "", 0., (pntau - 1.) / 2.);
			TMatrixDSym covm = result->GetCovarianceMatrix();
			double chi2 = result->Chi2(), ndof = result->Ndf();
			double cinf = result->Parameter(0), k = result->Parameter(1), csi_1 = result->Parameter(2);
			double dcinf = result->ParError(0), dk = result->ParError(1), dcsi_1 = result->ParError(2);
			double csi = 1. / csi_1 * beta / pntau, dcsi = dcsi_1 / (csi_1 * csi_1) * beta / pntau;
			graph[i]->GetFunction("func")->SetRange(0., pntau - 1.);
			
			std::fstream fstr;
			fstr.open((outname + "two_point_fit.txt").c_str(), std::fstream::out | std::fstream::app);
			std::streambuf* sb_cout = std::cout.rdbuf();
			std::streambuf* sb_file = fstr.rdbuf();
			std::cout.rdbuf(sb_file);
			result->Print("V");
			std::cout.rdbuf(sb_cout);
			
			multigraph->Add(graph[i]);
			
			std::ostringstream string_stream2;
			string_stream2 << "part " << i << std::fixed << std::setprecision(2) << "  csi = " << csi << "+-" << dcsi;
			std::string legend_entry = string_stream2.str();
			legend->AddEntry(graph[i], legend_entry.c_str(), "p");
		}
		
		c->SetGrid();
		// gPad->SetLogy();
		multigraph->SetTitle((folder + ";#tau;conn_corr").c_str());
		multigraph->Draw("A*");
		legend->SetTextSize(0.02);
		legend->Draw();
		c->SaveAs((outname + "two_point_plot.pdf").c_str());
		delete c;
		delete multigraph;
		delete legend;
	}
	{
		TCanvas* c = new TCanvas("two_point_canvas", "two_point_canvas", 2000, 1000);
		TGraphErrors* graph[nparticles];
		TMultiGraph* multigraph = new TMultiGraph();
		TLegend* legend = new TLegend(0.75, 0.75, 0.9, 0.85);
		
		for(int i = 0; i < nparticles; i++) {
			graph[i] = new TGraphErrors(pntau, ptau, pcorr2[i].data(), nullptr, dpcorr2[i].data());
			graph[i]->SetMarkerColorAlpha(i + 1, 0.65);
			graph[i]->SetLineColorAlpha(i + 1, 0.5);
			
			TF1* func = new TF1("func", "[0] + [1] * exp(- [2] * x)", 0., pntau - 1.);
			func->SetParameters(-0.05, 0.5, beta / pntau);
			func->SetParNames("c2inf","k","1/csi");
			func->SetNpx(pntau * 10);
			func->SetLineColorAlpha(i + 1, 1.);
			func->SetLineStyle(kDashed);
			
			TFitResultPtr result = graph[i]->Fit("func", "S+", "", 0., (pntau - 1.) / 3.);
			TMatrixDSym covm = result->GetCovarianceMatrix();
			double chi2 = result->Chi2(), ndof = result->Ndf();
			double c2inf = result->Parameter(0), k = result->Parameter(1), csi2_1 = result->Parameter(2);
			double dc2inf = result->ParError(0), dk = result->ParError(1), dcsi2_1 = result->ParError(2);
			double csi2 = 1. / csi2_1 * beta / pntau, dcsi2 = dcsi2_1 / (csi2_1 * csi2_1) * beta / pntau;
			graph[i]->GetFunction("func")->SetRange(0., pntau - 1.);
			
			std::fstream fstr;
			fstr.open((outname + "two_point2_fit.txt").c_str(), std::fstream::out | std::fstream::app);
			std::streambuf* sb_cout = std::cout.rdbuf();
			std::streambuf* sb_file = fstr.rdbuf();
			std::cout.rdbuf(sb_file);
			result->Print("V");
			std::cout.rdbuf(sb_cout);
			
			multigraph->Add(graph[i]);
			
			std::ostringstream string_stream2;
			string_stream2 << "part " << i << std::fixed << std::setprecision(2) << "  csi2 = " << csi2 << "+-" << dcsi2;
			std::string legend_entry = string_stream2.str();
			legend->AddEntry(graph[i], legend_entry.c_str(), "p");
		}
		
		c->SetGrid();
		// gPad->SetLogy();
		multigraph->SetTitle((folder + ";#tau;conn_corr").c_str());
		multigraph->Draw("A*");
		legend->SetTextSize(0.02);
		legend->Draw();
		c->SaveAs((outname + "two_point2_plot.pdf").c_str());
		delete c;
		delete multigraph;
		delete legend;
	}
}