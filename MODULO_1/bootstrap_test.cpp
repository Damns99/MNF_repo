#include <iostream>
#include "bootstrap.h"
#include "box_muller.h"
#include "vec_hist.h"
#include <vector>
#include <math.h>

#include <TGraph.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TFrame.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TNtuple.h>

float mean(const std::vector<float>& vec) {
	float s = 0.;
	int length = vec.size();
	for(auto& ii: vec) s += ii / length;
	return s;
}

float cumulBinder(const std::vector<float>& vec) {
	float s4 = 0., s2 = 0.;
	int length = vec.size();
	for(auto& ii: vec) {
		s4 += ii * ii * ii * ii / length;
		s2 += ii * ii / length;
	}
	return s4 / (3. * s2 * s2);
}

int main() {
	long seed = -942342;
	int n1 = 12;
	boxMuller::BoxMuller bm(seed);
	std::vector<float> vec = bm.randToVec(0., 1., pow(2, n1));
	auto hist = vecHist::makeHist(vec, 10, float(-3.), float(3.));
	
	TCanvas* canvas1 = new TCanvas("canvas1", "Canvas1", 600, 400);
	TH1F* hist_root = new TH1F("hist_root", "Hist_root", 10, -3., 3.);
	for(int ii = 0; ii < hist.size(); ii++) {
		int jj = hist[ii];
		float kk = -3. + (ii + 0.5) * 6. / 10;
		hist_root->Fill(kk, jj);
	}
	canvas1->cd();
	canvas1->SetGrid();
	hist_root->Draw("HIST");
	canvas1->SaveAs("bootstrap_test1.pdf");
	
	std::cout << bootstrap::bootstrapError(mean, vec, 1000, 1, seed) << std::endl;
	
	TCanvas* canvas2 = new TCanvas("canvas2", "Canvas2", 600, 400);
	int n2 = n1;
	double x2[n2], y2[n2];
	
	for(int ii = 0; ii < n2; ii++) {
		x2[ii] = pow(2, ii);
		y2[ii] = bootstrap::bootstrapError(mean, vec, 250, x2[ii], seed);
	}
	
	TGraph* graph2 = new TGraph(n2, x2, y2);
	canvas2->cd();
	canvas2->SetGrid();
	canvas2->SetLogx();
	graph2->Draw("AC*");
	canvas2->SaveAs("bootstrap_test2.pdf");
	
	TCanvas* canvas3 = new TCanvas("canvas3", "Canvas3", 600, 400);
	int n3 = 500;
	double x3[n3], y3[n3];
	
	for(int ii = 0; ii < n3; ii++) {
		x3[ii] = ii;
		y3[ii] = bootstrap::bootstrapError(mean, vec, x3[ii], 64, seed);
	}
	
	TGraph* graph3 = new TGraph(n3, x3, y3);
	canvas3->cd();
	canvas3->SetGrid();
	graph3->Draw("AC");
	canvas3->SaveAs("bootstrap_test3.pdf");
	
	std::cout << std::endl;
}