#include "lattice.h"
#include <iomanip>

#include <TCanvas.h>
#include <TAxis.h>
#include <TFrame.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TNtuple.h>
#include <TGraph.h>
#include <TMultiGraph.h>

double mean(const double* vec, int l) {
	double s = 0.;
	for(int ii = 0; ii < l; ii++) s += vec[ii] / l;
	return s;
}

void snapshot(const std::string out_image = "snapshot.png") {
	std::stringstream title;
	title << "beta = " << beta;
	TCanvas* c1 = new TCanvas("c1", out_image.c_str(), 2000, 1000);
	TGraph* graph[nparticles];
	TMultiGraph* multigraph = new TMultiGraph();
	TLegend* legend = new TLegend(0.75, 0.8, 0.9, 0.85);
	double time[p_length];
	for(int ii = 0; ii < p_length; ii++) time[ii] = 1. * ii / p_length * beta;
	for(int i = 0; i < nparticles; i++) {
		graph[i] = new TGraph(p_length, time, y + p_length * i);
		graph[i]->SetMarkerColor(i + 1);
		graph[i]->SetLineColor(i + 1);
		multigraph->Add(graph[i]);
		std::ostringstream string_stream2;
		string_stream2 << "part " << i << std::fixed << std::setprecision(4) << "  mean " << mean(y + p_length * i, p_length);
		std::string legend_entry = string_stream2.str();
		legend->AddEntry(graph[i], legend_entry.c_str(), "p");
	}
	c1->SetGrid();
	multigraph->SetTitle(";#tau;y");
	multigraph->Draw("AL*");
	legend->Draw();
	legend->SetTextSize(0.02);
	c1->SaveAs(out_image.c_str());
	delete c1;
	// for(int i = 0; i < nparticles; i++) delete graph[i];
	delete multigraph;
}