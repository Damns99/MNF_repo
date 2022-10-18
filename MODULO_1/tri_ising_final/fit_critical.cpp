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

#define BIG_DOUBLE 8192.

void fitGraphLine(const int n, double* x, double* y, double* dx, double* dy, const std::string xname, const std::string yname, const std::string title, const std::string outname, double xmin, double xmax, double* p0, const std::string m_var_name, double* m_results) {
	TCanvas* canvas = new TCanvas((outname + " canvas").c_str(), (outname + " canvas").c_str(), 600, 400);
	TGraphErrors* graph = new TGraphErrors(n, x, y, dx, dy);
	canvas->cd();
	canvas->SetGrid();
	graph->SetTitle((title + ";" + xname + ";" + yname).c_str());
	
	graph->Draw("AP");
	graph->SetMarkerStyle(kFullCircle);
	graph->SetMarkerSize(0.25);
	
	TF1* func = new TF1("func", "[0] + [1] * x", xmin, xmax);
	func->SetParameters(p0[0], p0[1]);
	func->SetParNames("q","m");
	func->SetNpx(100);
	
	func->SetLineColor(kRed);
	func->SetLineWidth(1);
	// func->DrawCopy("SAME");
	
	TFitResultPtr result = graph->Fit("func", "S+", "", xmin, xmax);
	TMatrixDSym covm = result->GetCovarianceMatrix();
	double chi2 = result->Chi2(), ndof = result->Ndf();
	double q = result->Parameter(0), m = result->Parameter(1);
	double dq = result->ParError(0), dm = result->ParError(1);
	
	std::fstream fstr;
	fstr.open((outname + "_fit.txt").c_str(), std::fstream::out);
	std::streambuf* sb_cout = std::cout.rdbuf();
	std::streambuf* sb_file = fstr.rdbuf();
	std::cout.rdbuf(sb_file);
	result->Print("V");
	std::cout.rdbuf(sb_cout);
	
	func->SetLineColor(kGreen);
	func->SetLineWidth(2);
	func->DrawCopy("SAME");
	canvas->SaveAs((outname + "_fit.pdf").c_str());
	
	std::string outfilename = "critical_fit_results.txt";
	std::ofstream outfile;
	outfile.open(outfilename, std::fstream::out | std::fstream::app);
	outfile << m_var_name << " = " << m << " +- " << dm << std::endl;
	outfile.close();
	
	m_results[0] = m;
	m_results[1] = dm;
}

void fitGraph3Variables(const int n, double* x, double* y, double* dx, double* dy, const std::string xname, const std::string yname, const std::string title, const std::string outname, double xmin, double xmax, double* p0, double* beta_c_results, double* _nu_results) {
	TCanvas* canvas = new TCanvas((outname + " canvas").c_str(), (outname + " canvas").c_str(), 600, 400);
	TGraphErrors* graph = new TGraphErrors(n, x, y, dx, dy);
	canvas->cd();
	canvas->SetGrid();
	graph->SetTitle((title + ";" + xname + ";" + yname).c_str());
	
	graph->Draw("AP");
	graph->SetMarkerStyle(kFullCircle);
	graph->SetMarkerSize(0.25);
	
	TF1* func = new TF1("func", "[0] + [1] * (x ^ (-[2]))", xmin, xmax);
	func->SetParameters(p0[0], p0[1], p0[2]);
	func->SetParNames("beta_c","x_max","_nu");
	func->SetNpx(100);
	
	func->SetLineColor(kRed);
	func->SetLineWidth(1);
	// func->DrawCopy("SAME");
	
	TFitResultPtr result = graph->Fit("func", "S+", "", xmin, xmax);
	TMatrixDSym covm = result->GetCovarianceMatrix();
	double chi2 = result->Chi2(), ndof = result->Ndf();
	double beta_c = result->Parameter(0), x_max = result->Parameter(1), _nu = result->Parameter(2);
	double dbeta_c = result->ParError(0), dx_max = result->ParError(1), d_nu = result->ParError(2);
	
	std::fstream fstr;
	fstr.open((outname + "_fit.txt").c_str(), std::fstream::out);
	std::streambuf* sb_cout = std::cout.rdbuf();
	std::streambuf* sb_file = fstr.rdbuf();
	std::cout.rdbuf(sb_file);
	result->Print("V");
	std::cout.rdbuf(sb_cout);
	
	func->SetLineColor(kGreen);
	func->SetLineWidth(2);
	func->DrawCopy("SAME");
	canvas->SaveAs((outname + "_fit.pdf").c_str());
	
	std::string outfilename = "critical_fit_results.txt";
	std::ofstream outfile;
	outfile.open(outfilename, std::fstream::out | std::fstream::app);
	outfile << "beta_c" << " = " << beta_c << " +- " << dbeta_c << std::endl;
	outfile << "_nu" << " = " << _nu << " +- " << d_nu << std::endl;
	outfile.close();
	
	beta_c_results[0] = beta_c;
	beta_c_results[1] = dbeta_c;
	_nu_results[0] = _nu;
	_nu_results[1] = d_nu;
}

int main(int argc, char* argv[]) {
	std::string outfilename, inname;
	double Cxmin, Cxmax, Kxmin, Kxmax, Kbetaxmin, Kbetaxmax;
	double Cp0[2], Kp0[2], Kbetap0[3];
	
	cmdlineParser::CmdlineParser parser;
	parser.addOptParameter<std::string>("infile", &inname, "beta_results", "[std::string] Input file name prefix.");
	parser.addOptParameter<std::string>("outfile", &outfilename, "critical_indices", "[std::string] Output file name.");
    parser.addOptParameter<double>("Cxmin", &Cxmin, -BIG_DOUBLE, "[double] Min x for C fitting.");
    parser.addOptParameter<double>("Cxmax", &Cxmax, BIG_DOUBLE, "[double] Max x for C fitting.");
    parser.addOptParameter<double>("Kxmin", &Kxmin, -BIG_DOUBLE, "[double] Min x for K fitting.");
    parser.addOptParameter<double>("Kxmax", &Kxmax, BIG_DOUBLE, "[double] Max x for K fitting.");
    parser.addOptParameter<double>("Kbetaxmin", &Kbetaxmin, -BIG_DOUBLE, "[double] Min x for Kbeta fitting.");
    parser.addOptParameter<double>("Kbetaxmax", &Kbetaxmax, BIG_DOUBLE, "[double] Max x for Kbeta fitting.");
    parser.addOptArrayParameter<double>("Cp0", Cp0, 1., "[double [2]] Init parameters for C fitting.", 2);
    parser.addOptArrayParameter<double>("Kp0", Kp0, 1., "[double [2]] Init parameters for K fitting.", 2);
    parser.addOptArrayParameter<double>("Kbetap0", Kbetap0, 1., "[double [3]] Init parameters for K fitting.", 3);
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
	
	std::vector<double> logL, logC, dlogC, logK, dlogK;
	for(auto& ii: L) logL.push_back(log(ii));
	for(auto& ii: Cmax) logC.push_back(log(ii));
	for(int i = 0; i < Csize; i++) dlogC.push_back(dCmax[i] / Cmax[i]);
	for(auto& ii: Kmax) logK.push_back(log(ii));
	for(int i = 0; i < Ksize; i++) dlogK.push_back(dKmax[i] / Kmax[i]);
	
	double m_results[2];
	fitGraphLine(Csize, logL.data(), logC.data(), nullptr, dlogC.data(), "ln(L)", "ln(#C_{max})", "fit #alpha/#nu", outfilename + "_C", Cxmin, Cxmax, Cp0, "alpha/nu", m_results);
	double alpha_nu = m_results[0], dalpha_nu = m_results[1];
	fitGraphLine(Ksize, logL.data(), logK.data(), nullptr, dlogK.data(), "ln(L)", "ln(#K_{max})", "fit #gamma/#nu", outfilename + "_K", Kxmin, Kxmax, Kp0, "gamma/nu", m_results);
	double gamma_nu = m_results[0], dgamma_nu = m_results[1];
	double beta_c_results[2], _nu_results[2];
	fitGraph3Variables(Ksize, L.data(), Kbetamax.data(), nullptr, dKbetamax.data(), "L", "#beta_{max}", "fit #beta_c e 1/#nu", outfilename + "_Kbeta", Kbetaxmin, Kbetaxmax, Kbetap0, beta_c_results, _nu_results);
	
	std::cout << std::endl << std::endl;
	std::cout << "beta_c = " << beta_c_results[0] << " +- " << beta_c_results[1] << std::endl;
	std::cout << "nu = " << 1. / _nu_results[0] << " +- " << _nu_results[1] / (_nu_results[0] * _nu_results[0]) << std::endl;
	std::cout << "alpha = " << alpha_nu / _nu_results[0] << " +- " << (dalpha_nu / alpha_nu + _nu_results[1] / _nu_results[0]) * alpha_nu / _nu_results[0] << std::endl;
	std::cout << "gamma = " << gamma_nu / _nu_results[0] << " +- " << (dgamma_nu / gamma_nu + _nu_results[1] / _nu_results[0]) * gamma_nu / _nu_results[0] << std::endl;
}