#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <filesystem>
#include <cassert>
namespace fs = std::filesystem;

#include "text_io.h"
#include "cmdline_parser.h"

#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TFrame.h>
#include <TLegend.h>
#include <TNtuple.h>
#include <TF1.h>
#include <TMatrixDSym.h>
#include <TFitResult.h>

void makeGraph(const int n, double* x, double* y, double* dx, double* dy, const std::string xname, const std::string yname, const std::string title, const std::string outname) {
	TCanvas* canvas = new TCanvas(title.c_str(), title.c_str(), 600, 400);
	canvas->cd();
	canvas->SetGrid();
	TGraphErrors* graph = new TGraphErrors(n, x, y, dx, dy);
	graph->SetTitle((title + ";" + xname + ";" + yname).c_str());
	graph->Draw("AP");
	graph->SetMarkerStyle(kFullCircle);
	graph->SetMarkerSize(0.25);
	canvas->SaveAs((outname + ".pdf").c_str());
}

int main(int argc, char* argv[]) {
	std::string outfilename, inname;
	double Cxmin, Cxmax, Kxmin, Kxmax, Kbetaxmin, Kbetaxmax;
	double Cp0[2], Kp0[2], Kbetap0[3];
	
	cmdlineParser::CmdlineParser parser;
	parser.addOptParameter<std::string>("infile", &inname, "beta_results", "[std::string] Input file name prefix.");
	parser.addOptParameter<std::string>("outfile", &outfilename, "fsscaling", "[std::string] Output file name.");
	if (parser.parseAll(argc, argv) == HELP_RETURN) return 0;
    parser.kickOff(argv[0]);
	
	std::vector<double> L;
	std::vector<double> tmp, dtmp;
	// std::vector<double> Emax, dEmax, Ebetamax, dEbetamax;
	// std::vector<double> Mmax, dMmax, Mbetamax, dMbetamax;
	std::vector<double> Cmax, dCmax, Cbetamax, dCbetamax;
	std::vector<double> Kmax, dKmax, Kbetamax, dKbetamax;
	
	fs::current_path(fs::current_path() / "measures");
    // int Esize = textIo::textIn(inname + "_E_fit_results", '\t', '#', &L, &Emax, &dEmax, &tmp, &dtmp, &Ebetamax, &dEbetamax);
    // int Msize = textIo::textIn(inname + "_M_fit_results", '\t', '#', &L, &Mmax, &dMmax, &tmp, &dtmp, &Mbetamax, &dMbetamax);
    int Csize = textIo::textIn(inname + "_C_fit_results.txt", '\t', '#', &L, &Cmax, &dCmax, &tmp, &dtmp, &Cbetamax, &dCbetamax);
    int Ksize = textIo::textIn(inname + "_K_fit_results.txt", '\t', '#', &L, &Kmax, &dKmax, &tmp, &dtmp, &Kbetamax, &dKbetamax);
	
	/* std::vector<double> logL, logC, dlogC, logK, dlogK;
	for(auto& ii: L) logL.push_back(log(ii));
	for(auto& ii: Cmax) logC.push_back(log(ii));
	for(int i = 0; i < Csize; i++) dlogC.push_back(dCmax[i] / Cmax[i]);
	for(auto& ii: Kmax) logK.push_back(log(ii));
	for(int i = 0; i < Ksize; i++) dlogK.push_back(dKmax[i] / Kmax[i]); */
	
	std::vector<double> K_y, dK_y;
	double gamma_nu = 7./4.;
	for(int i = 0; i < Ksize; i++) {
		double lgn = pow(L[i], gamma_nu);
		K_y.push_back(Kmax[i] / lgn);
		dK_y.push_back(dKmax[i] / lgn);
	}
	makeGraph(Ksize, L.data(), K_y.data(), nullptr, dK_y.data(), "L", "K_{max}L^{-#frac{#gamma}{#nu}}", "K finite size scaling theory", outfilename+"_K_theory");
	
	std::vector<double> K_y_exp, dK_y_exp;
	double gamma_nu_exp = 1.709;
	for(int i = 0; i < Ksize; i++) {
		double lgn = pow(L[i], gamma_nu_exp);
		K_y_exp.push_back(Kmax[i] / lgn);
		dK_y_exp.push_back(dKmax[i] / lgn);
	}
	makeGraph(Ksize, L.data(), K_y_exp.data(), nullptr, dK_y_exp.data(), "L", "K_{max}L^{-#frac{#gamma}{#nu}}", "K finite size scaling experiment", outfilename+"_K_exp");
}