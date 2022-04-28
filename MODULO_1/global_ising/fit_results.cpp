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
#include <TF1.h>
#include <TMatrixDSym.h>
#include <TFitResult.h>

void fitGraph(const int n, double* x, double* y, double* dx, double* dy, const std::string xname, const std::string yname, const std::string title, const std::string outname, double xmin, double xmax, double* p0, int append, int length) {
	TCanvas* canvas = new TCanvas((outname + "canvas").c_str(), (outname + "canvas").c_str(), 600, 400);
	TGraphErrors* graph = new TGraphErrors(n, x, y, dx, dy);
	canvas->cd();
	canvas->SetGrid();
	graph->SetTitle((title + ";" + xname + ";" + yname).c_str());
	
	graph->Draw("A*");
	
	TF1* func = new TF1("func", "[0] + [1] * (x - [2]) * (x - [2])", xmin, xmax);
	func->SetParameters(p0[0], p0[1], p0[2]);
	func->SetParNames("y_max","c","x_max");
	func->SetNpx(100);
	
	func->SetLineColor(kRed);
	func->DrawCopy("SAME");
	
	TFitResultPtr result = graph->Fit("func", "S+", "", xmin, xmax);
	// TF1* fitted_func = new TF1(*graph->GetFunction("func"));
	TMatrixDSym covm = result->GetCovarianceMatrix();
	double chi2 = result->Chi2(), ndof = result->Ndf();
	double K_max = result->Parameter(0), c = result->Parameter(1), b_max = result->Parameter(2);
	double dK_max = result->ParError(0), dc = result->ParError(1), db_max = result->ParError(2);
	
	std::fstream fstr;
	fstr.open((outname + "_fit.txt").c_str(), std::fstream::out);
	std::streambuf* sb_cout = std::cout.rdbuf();
	std::streambuf* sb_file = fstr.rdbuf();
	std::cout.rdbuf(sb_file);
	result->Print("V");
	std::cout.rdbuf(sb_cout);
	
	func->SetLineColor(kGreen);
	func->DrawCopy("SAME");
	canvas->SaveAs((outname + "_fit.pdf").c_str());
	
	auto cpath = fs::current_path();
    fs::current_path(cpath / "..");
	std::string outfilename = outname + "_fit_results.txt";
	std::ofstream outfile;
	if (append) outfile.open(outfilename, std::fstream::out | std::fstream::app);
	else {
		outfile.open(outfilename, std::fstream::out);
		outfile << "#L \ty_max \terror \tc \terror \tx_max \terror" << std::endl;
	}
	outfile << length << "\t";
	outfile << K_max << "\t" << dK_max << "\t" << c << "\t" << dc << "\t" << b_max << "\t" << db_max << std::endl;
	outfile.close();
    fs::current_path(cpath);
}

int main(int argc, char* argv[]) {
    std::string folder, infilename, outname, xname;
    std::vector<double> x, E, dE, M, dM, C, dC, K, dK, BE, dBE, BM, dBM;
	int append, len;
	double xminC, xmaxC, xminK, xmaxK;
	
	cmdlineParser::CmdlineParser parser;
    parser.addPosParameter<std::string>("folder", &folder, "0220419144203", "[std::string] Working folder name.");
    parser.addPosParameter<std::string>("infile", &infilename, "results.txt", "[std::string] Input file for the results to plot.");
    parser.addPosParameter<int>("L", &len, 10, "[int] lattice length.");
    parser.addPosParameter<double>("xminC", &xminC, 0.4, "[double] minimum x to account for fitting C.");
    parser.addPosParameter<double>("xmaxC", &xmaxC, 0.5, "[double] maximum x to account for fitting C.");
    parser.addPosParameter<double>("xminK", &xminK, 0.4, "[double] minimum x to account for fitting K.");
    parser.addPosParameter<double>("xmaxK", &xmaxK, 0.5, "[double] maximum x to account for fitting K.");
    parser.addOptParameter<std::string>("outname", &outname, "", "[std::string] Output files name prefix.");
    parser.addOptParameter<std::string>("xname", &xname, "", "[std::string] Name for the x axis.");
    parser.addOptParameter<int>("append", &append, 0, "[int] If != 0 append means to measfile instead of overwrtiting them.");
	if (parser.parseAll(argc, argv) == HELP_RETURN) return 0;
    parser.kickOff(argv[0]);
	
    fs::current_path(fs::current_path() / "measures" / folder);
    int nmeas = textIo::textIn(infilename, '\t', '#', &x, &E, &dE, &M, &dM, &C, &dC, &K, &dK, &BE, &dBE, &BM, &dBM);
	
	double p0C[3] = {12., -25000., 0.44};
	double p0K[3] = {0.04 * len * len, -100. * len * len, 0.44};
	
	fitGraph(nmeas, x.data(), C.data(), nullptr, dC.data(), xname, "d#varepsilon/dT", folder, outname + "C", xminC, xmaxC, p0C, append, len);
	fitGraph(nmeas, x.data(), K.data(), nullptr, dK.data(), xname, "d|M|/dh", folder, outname + "K", xminK, xmaxK, p0K, append, len);
}