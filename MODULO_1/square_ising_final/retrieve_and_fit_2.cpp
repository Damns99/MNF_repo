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

#define N_MAX_FOLDERS 64
#define BIG_DOUBLE 8192.
#define MIN_POINTS_TO_FIT 4
#define PERCENTAGE_STEP 0.05

double find_max(double* array, int length, int& index) {
	double max = array[0];
	index = 0;
	for(int i = 1; i < length; i++) {
		if(max < array[i]) {
			max = array[i];
			index = i;
		}
	}
	return max;
}

void fitGraphParable(const int n, double* x, double* y, double* dx, double* dy, const std::string xname, const std::string yname, const std::string title, const std::string outname, int append, int length, double percentage) {
	TCanvas* canvas = new TCanvas((outname + " canvas").c_str(), (outname + " canvas").c_str(), 600, 400);
	TGraphErrors* graph = new TGraphErrors(n, x, y, dx, dy);
	canvas->cd();
	canvas->SetGrid();
	graph->SetTitle((title + ";" + xname + ";" + yname).c_str());
	
	int indmax;
	double ymax = find_max(y, n, indmax);
	double xmax = x[indmax];
	int counter;
	double xleft, xright;
	double p0[3];
	do {
		counter = 0;
		for(int ii = 0; ii < n; ii++) {
			if(y[ii] >= (1. - percentage) * ymax) {
				if(counter == 0) xleft = x[ii];
				xright = x[ii];
				counter++;
			}
		}
		if(counter < MIN_POINTS_TO_FIT) {
			percentage += PERCENTAGE_STEP;
			std::cout << "! changed percentage to " << percentage << " !" << std::endl;
		}
	} while(counter < MIN_POINTS_TO_FIT && percentage <= 1.);
	p0[0] = ymax;
	p0[2] = xmax;
	p0[1] = ymax*percentage/((xleft-xmax)*(xleft-xmax));
	std::cout << "p0 = ( " << p0[0] << ", " << p0[1] << ", " << p0[2] << " )"  << std::endl;
	
	TF1* func = new TF1("func", "[0] - [1] * (x - [2]) * (x - [2])", xleft, xright);
	func->SetParameters(p0[0], p0[1], p0[2]);
	func->SetParLimits(1, 0., BIG_DOUBLE * BIG_DOUBLE);
	func->SetParLimits(2, xleft, xright);
	func->SetParNames("y_max","c","x_max");
	func->SetNpx(100);
	
	func->SetLineColor(kRed);
	func->SetLineWidth(1);
	
	TFitResultPtr result = graph->Fit("func", "SR+", "", xleft, xright);
	TMatrixDSym covm = result->GetCovarianceMatrix();
	double chi2 = result->Chi2(), ndof = result->Ndf();
	double y_max = result->Parameter(0), c = result->Parameter(1), x_max = result->Parameter(2);
	double dy_max = result->ParError(0), dc = result->ParError(1), dx_max = result->ParError(2);
	
	std::fstream fstr;
	fstr.open((outname + "_fit.txt").c_str(), std::fstream::out);
	std::streambuf* sb_cout = std::cout.rdbuf();
	std::streambuf* sb_file = fstr.rdbuf();
	std::cout.rdbuf(sb_file);
	std::cout << "xlim = [ " << xleft << " , " << xright << " ]" << std::endl;
	result->Print("V");
	std::cout.rdbuf(sb_cout);
	
	func->SetLineColor(kGreen);
	func->SetLineWidth(2);
	func->DrawCopy("SAME");
	
	graph->Draw("AP");
	graph->SetMarkerStyle(kFullCircle);
	graph->SetMarkerSize(0.25);
	
	canvas->SaveAs((outname + "_fit.pdf").c_str());
	
	std::string outfilename = outname + "_fit_results.txt";
	std::ofstream outfile;
	if (append) outfile.open(outfilename, std::fstream::out | std::fstream::app);
	else {
		outfile.open(outfilename, std::fstream::out);
		outfile << "#L \ty_max \terror \tc \terror \tx_max \terror" << std::endl;
	}
	outfile << length << "\t";
	outfile << y_max << "\t" << dy_max << "\t" << c << "\t" << dc << "\t" << x_max << "\t" << dx_max << std::endl;
	outfile.close();
}

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
	std::string outname, infilename, measfilename, folders[N_MAX_FOLDERS];
	int nfolders, length, append;
	double percC, percK;
	
	cmdlineParser::CmdlineParser parser;
    parser.addPosParameter<int>("L", &length, 0, "[int] lattice length.");
    parser.addPosParameter<int>("nf", &nfolders, N_MAX_FOLDERS, "[int] Number of folders to consider.");
    parser.addPosArrayParameter<std::string>("folders", folders, "", "[std::string []] Folder list.", N_MAX_FOLDERS);
	parser.addOptParameter<std::string>("infile", &infilename, "results.txt", "[std::string] Input file name to search.");
	parser.addOptParameter<std::string>("measfile", &measfilename, "metro_ising_meas.txt", "[std::string] Input file for beta.");
	parser.addOptParameter<std::string>("outname", &outname, "beta_results_2", "[std::string] Output file name prefix for the results.");
    parser.addOptParameter<int>("append", &append, 0, "[int] If != 0 append means to measfile instead of overwrtiting them.");
    parser.addOptParameter<double>("percC", &percC, 0.2, "[double] Percentage of the maximum y value for C to consider in fit.");
    parser.addOptParameter<double>("percK", &percK, 0.2, "[double] Percentage of the maximum y value for C to consider in fit.");
	if (parser.parseAll(argc, argv) == HELP_RETURN) return 0;
    parser.kickOff(argv[0]);
	
	assert(nfolders < N_MAX_FOLDERS);
	
	std::vector<double> beta(nfolders), E(nfolders), dE(nfolders), M(nfolders), dM(nfolders), K(nfolders), dK(nfolders), C(nfolders), dC(nfolders);
	std::vector<double> BE(nfolders), dBE(nfolders), BM(nfolders), dBM(nfolders);
	
	fs::current_path(fs::current_path() / "measures");
    auto measures_folder = fs::current_path();
	for(int i = 0; i < nfolders; i++) {
		fs::current_path(fs::current_path() / folders[i]);
		
		std::string tmp;
		
		std::ifstream tmpif;
		tmpif.open(measfilename, std::fstream::in);
		tmpif >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> beta[i];
		tmpif.close();
		
		tmpif.open(infilename, std::fstream::in);
		tmpif >> tmp >> tmp >> E[i] >> tmp >> dE[i];
		tmpif >> tmp >> tmp >> M[i] >> tmp >> dM[i];
		tmpif >> tmp >> tmp >> C[i] >> tmp >> dC[i];
		tmpif >> tmp >> tmp >> K[i] >> tmp >> dK[i];
		tmpif >> tmp >> tmp >> BE[i] >> tmp >> dBE[i];
		tmpif >> tmp >> tmp >> BM[i] >> tmp >> dBM[i];
		tmpif.close();
		
		fs::current_path(measures_folder);
	}
	
	std::string output_folder = "L" + std::to_string(length);
	fs::create_directory(output_folder);
	fs::current_path(fs::current_path() / output_folder);
	
	textIo::textOut(outname + ".txt", '\t', '#', " beta \t<E> \terror \t<|M|> \terror \td<E>/dT \terror \td<|M|>/dh \terror", nfolders, false, beta, E, dE, M, dM, C, dC, K, dK);
	
	fitGraphParable(nfolders, beta.data(), C.data(), nullptr, dC.data(), "#beta", "C(#beta)", "C(#beta) @ L = " + std::to_string(length), outname + "_C", append, length, percC);
	fitGraphParable(nfolders, beta.data(), K.data(), nullptr, dK.data(), "#beta", "K(#beta)", "K(#beta) @ L = " + std::to_string(length), outname + "_K", append, length, percK);
	
	makeGraph(nfolders, beta.data(), BE.data(), nullptr, dBE.data(), "#beta", "Binder(E)", "Binder cumulant for E", outname+"_BE");
	makeGraph(nfolders, beta.data(), BM.data(), nullptr, dBM.data(), "#beta", "Binder(M)", "Binder cumulant for M", outname+"_BM");
}
	
	
    