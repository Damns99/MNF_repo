#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <filesystem>
namespace fs = std::filesystem;
#include <fstream>

#include <TGraphErrors.h>
#include <TF1.h>
#include <TMatrixDSym.h>
#include <TFitResult.h>

#include "bound_cond_vecs.h"
#include "wave_plots.h"
#include "mystyle.h"
#include "text_io.h"

int main() {
	fs::current_path(fs::current_path() / "measures");
	myStyle();
	
	std::string filename = "HH_test_max.txt";
	
	std::vector<double> x, t, V;
	textIo::textIn(filename, '\t', '#', &x, &t, &V);
	x.pop_back(); t.pop_back(); V.pop_back();
	int N = x.size();
	
	std::ifstream infile;
	infile.open(filename, std::fstream::in);
	std::string tmp;
	double ddx, ddt, center, stdev, ampl;
	infile >> tmp >> tmp >> tmp >> tmp >> ddx >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> ddt;
	infile >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> center >> tmp >> stdev >> tmp >> ampl;	
	infile.close();
	std::vector<double> dx(N, ddx/2), dt(N, ddt/2);
	
	double maxdist = abs(x[N-1] - x[0]) + ddx;
	for(int i = 0; i < N; i++) {
		double tmpx = abs(x[i]-center);
		if(tmpx > maxdist / 2.) tmpx = maxdist - tmpx;
		x[i] = tmpx;
	}
	
	TCanvas* canvas = new TCanvas("canvas", "Canvas", 600, 400);
	TGraphErrors* graph = new TGraphErrors(N, x.data(), t.data(), dx.data(), dt.data());
	graph->Sort();
	graph->SetTitle("peak positions;x [cm];t [ms]");
	canvas->cd();
	canvas->SetGrid();
	// canvas->SetLogx();
	// canvas->SetLogy();
	graph->Draw("AP");
	
	/*TF1* func = new TF1("func", "[0] * (1 + (x/[1])^2)^(1./2.) + [2]", 0, 100);
	func->SetParameters(5.63, 17, 0.);
	func->SetParNames("a", "b", "c");*/
	TF1* func = new TF1("func", "[0] * x + [1]", 0, 100);
	func->SetParameters(1., 0.);
	func->SetParNames("m", "q");
	func->SetNpx(100);
	
	func->SetLineColor(kRed);
	// func->DrawCopy("LSAME");
	
	TFitResultPtr result = graph->Fit("func", "S+0", "", 60, 90);
	//TF1* fitted_func = new TF1(*graph->GetFunction("func"));
	TMatrixDSym covm = result->GetCovarianceMatrix();
	double chi2 = result->Chi2(), ndof = result->Ndf();
	double m = result->Parameter(0), q = result->Parameter(1);
	double dm = result->ParError(0), dq = result->ParError(1);
	
	result->Print("V");
	
	func->SetLineColor(kGreen);
	func->DrawCopy("LSAME");
	
	double c = graph->GetPointY(1), dc = ddt;
	TF1* func2 = new TF1("func2", "[0] * x + [1]", 0, 100);
	func2->SetParameters(0., c);
	func2->SetParNames("par0", "c");
	func2->SetNpx(100);
	
	func2->SetLineColor(kBlue);
	func2->DrawCopy("LSAME");
	
	canvas->SaveAs("HH_test_max.pdf");	
	
	/* double* xsorted = graph->GetX(), * tsorted = graph->GetY();
	std::vector<double> xsortedv(N), tsortedv(N);
	for(int iii = 0; iii < N; iii++) {
		xsortedv[iii] = xsorted[iii];
		tsortedv[iii] = tsorted[iii];
	}
	int iii = 1;
	while(iii < xsortedv.size()) {
		if(abs(xsortedv[iii] - xsortedv[iii-1]) < ddx/2) {
			xsortedv.erase(xsortedv.begin()+iii);
			tsortedv.erase(tsortedv.begin()+iii);
		}
		else {
			iii++;
		}
	}
	// for(int iii = 0; iii < xsortedv.size(); iii++) std::cout << xsortedv[iii] << " " << tsortedv[iii] << std::endl;
	
	std::vector<double> derivative(xsortedv.size(), 0.);
	for(int iii = 1; iii < xsortedv.size(); iii++) {
		derivative[iii] = (tsortedv[iii] - tsortedv[iii-1]);
	}
	TGraph* graph2 = new TGraph(xsortedv.size(), xsortedv.data(), derivative.data());
	graph2->Draw("AP");
	canvas->SetLogx();
	canvas->SetLogy();
	canvas->SaveAs("HH_test_max_2.pdf"); */
	
	std::ofstream outfile;
	outfile.open("HH_max_gaussian_fitres.txt", std::fstream::out | std::fstream::app);
	outfile << ampl << "\t" << stdev << "\t" << m << "\t" << dm << "\t" << q << "\t" << dq << "\t" << c << "\t" << dc << std::endl; 
	outfile.close();
}