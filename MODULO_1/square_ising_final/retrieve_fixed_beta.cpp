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

#define N_MAX_FOLDERS 16
#define BIG_DOUBLE 8192.

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
	int nfolders, append;
	double beta;
	
	cmdlineParser::CmdlineParser parser;
    parser.addPosParameter<double>("beta", &beta, 0., "[double] Selected temperature.");
    parser.addPosParameter<int>("nf", &nfolders, N_MAX_FOLDERS, "[int] Number of folders to consider.");
    parser.addPosArrayParameter<std::string>("folders", folders, "", "[std::string []] Folder list.", N_MAX_FOLDERS);
	parser.addOptParameter<std::string>("infile", &infilename, "results.txt", "[std::string] Input file name to search.");
	parser.addOptParameter<std::string>("measfile", &measfilename, "metro_ising_meas.txt", "[std::string] Input file for beta.");
	parser.addOptParameter<std::string>("outname", &outname, "length_results_2", "[std::string] Output file name prefix for the results.");
    parser.addOptParameter<int>("append", &append, 0, "[int] If != 0 append means to measfile instead of overwrtiting them.");
	if (parser.parseAll(argc, argv) == HELP_RETURN) return 0;
    parser.kickOff(argv[0]);
	
	assert(nfolders < N_MAX_FOLDERS);
	
	std::vector<double> length(nfolders), E(nfolders), dE(nfolders), M(nfolders), dM(nfolders), K(nfolders), dK(nfolders), C(nfolders), dC(nfolders);
	std::vector<double> BE(nfolders), dBE(nfolders), BM(nfolders), dBM(nfolders);
	
	fs::current_path(fs::current_path() / "measures");
    auto measures_folder = fs::current_path();
	for(int i = 0; i < nfolders; i++) {
		fs::current_path(fs::current_path() / folders[i]);
		
		std::string tmp;
		
		std::ifstream tmpif;
		tmpif.open(measfilename, std::fstream::in);
		tmpif >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> length[i] >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp;
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
	std::vector<double> _length(nfolders);
	for(int i = 0; i < nfolders; i++) _length[i] = 1. / length[i];
	
	std::string output_folder = "b" + std::to_string(beta);
	fs::create_directory(output_folder);
	fs::current_path(fs::current_path() / output_folder);
	
	textIo::textOut(outname + ".txt", '\t', '#', " L \t<E> \terror \t<|M|> \terror \td<E>/dT \terror \td<|M|>/dh \terror", nfolders, false, length, E, dE, M, dM, C, dC, K, dK);
	
	makeGraph(nfolders, _length.data(), C.data(), nullptr, dC.data(), "1/L", "C(L)", "C(L) @ #beta = " + std::to_string(beta), outname + "_C");
	makeGraph(nfolders, _length.data(), K.data(), nullptr, dK.data(), "1/L", "K(L)", "K(L) @ #beta = " + std::to_string(beta), outname + "_K");
	
	makeGraph(nfolders, _length.data(), E.data(), nullptr, dE.data(), "1/L", "E", "Energy", outname+"_E");
	makeGraph(nfolders, _length.data(), M.data(), nullptr, dM.data(), "1/L", "|M|", "Magnetization", outname+"_M");
	
	makeGraph(nfolders, _length.data(), BE.data(), nullptr, dBE.data(), "1/L", "Binder(E)", "Binder cumulant for E", outname+"_BE");
	makeGraph(nfolders, _length.data(), BM.data(), nullptr, dBM.data(), "1/L", "Binder(|M|)", "Binder cumulant for M", outname+"_BM");
}
	
	
    