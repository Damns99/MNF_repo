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

constexpr int MAX_LENGTH = 65536;

void fitCorrDist(const int n, double* x, double* y, double* dx, double* dy, const std::string xname, const std::string yname, const std::string title, const std::string outname, double xmin, double xmax, double* p0, int append, int length) {
	TCanvas* canvas = new TCanvas((outname + "canvas").c_str(), (outname + "canvas").c_str(), 600, 400);
	TGraphErrors* graph = new TGraphErrors(n, x, y, dx, dy);
	canvas->cd();
	canvas->SetGrid();
	graph->SetTitle((title + ";" + xname + ";" + yname).c_str());
	
	graph->Draw("A*");
	
	TF1* func = new TF1("func", "[0] + [1] * x", xmin, xmax);
	func->SetParameters(p0[0], p0[1]);
	func->SetParNames("x0", "1/csi");
	func->SetNpx(100);
	
	func->SetLineColor(kRed);
	func->DrawCopy("SAME");
	
	TFitResultPtr result = graph->Fit("func", "S+", "", xmin, xmax);
	// TF1* fitted_func = new TF1(*graph->GetFunction("func"));
	TMatrixDSym covm = result->GetCovarianceMatrix();
	double chi2 = result->Chi2(), ndof = result->Ndf();
	double x0 = result->Parameter(0), csi = 1. / result->Parameter(1);
	double dx0 = result->ParError(0), dcsi = result->ParError(1) / (result->Parameter(1) * result->Parameter(1));
	
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
		outfile << "#L \tx0 \terror \tcsi \terror" << std::endl;
	}
	outfile << length << "\t";
	outfile << x0 << "\t" << dx0 << "\t" << csi << "\t" << dcsi << std::endl;
	outfile.close();
    fs::current_path(cpath);
}

int distsq(int i, int j, int lps) {
	switch(lps) {
		case 4:
			return (i - j) * (i - j);
		case 3:
			return 0;
		case 6:
			return 0;
		default:
			return 0;
	}
}

int main(int argc, char* argv[]) {
    std::string folder, infilename, outname;
	int spin[MAX_LENGTH];
	int len, lps, append;
	double xmin, xmax;
	
	cmdlineParser::CmdlineParser parser;
    parser.addPosParameter<std::string>("folder", &folder, "0220419144203", "[std::string] Working folder name.");
    parser.addPosParameter<std::string>("infile", &infilename, "metro_ising_conf.txt", "[std::string] Input file for the configuration to elaborate.");
    parser.addPosParameter<double>("xmin", &xmin, 0., "[double] minimum x to account for fitting.");
    parser.addPosParameter<double>("xmax", &xmax, 100., "[double] maximum x to account for fitting.");
    parser.addOptParameter<std::string>("outname", &outname, "", "[std::string] Output files name prefix.");
    parser.addOptParameter<int>("append", &append, 0, "[int] If != 0 append means to measfile instead of overwrtiting them.");
	if (parser.parseAll(argc, argv) == HELP_RETURN) return 0;
    parser.kickOff(argv[0]);
	
    fs::current_path(fs::current_path() / "measures" / folder);
    std::ifstream infile;
	infile.open(infilename, std::fstream::in);
	infile >> len;
	for(int i = 0; i < len * len; i++) infile >> spin[i];
	infile >> lps;
	infile.close();
	
	int distsq_max = 0;
	switch(lps) {
		case 4:
			distsq_max = len / 2 + 1;
			break;
		case 3:
			break;
		case 6:
			break;
		default:
			break;
	}
	
    std::vector<double> x(distsq_max, 0.), y(distsq_max, 0.), dy(distsq_max, 0.);
	std::vector<int> cy(distsq_max, 0);
	
	for(int i = 0; i < len * len; i++) for(int j = i; j < len * len; j++) {
		int index = distsq(i, j, lps);
		cy[index] += 1;
		double newvalue = 2. * spin[i] * spin[j];
		double delta1 = y[index] - newvalue;
		y[index] += delta1 / cy[index];
		double delta2 = y[index] - newvalue;
		dy[index] += delta1 * delta2;
	}
	
	for(int i = 0; i < distsq_max; i++) {
		x[i] = sqrt(i);
		y[i] = -log(y[i]);
		dy[i] = -log(sqrt(dy[i] / cy[i] / (cy[i] - 1)));
	}
	
	double p0[2] = {0., 0.};
	fitCorrDist(distsq_max, x.data(), y.data(), nullptr, dy.data(), "|i - j|", "-ln(<s_i s_j>)", folder, outname + "_corr_dist", xmin, xmax, p0, append, len);
}