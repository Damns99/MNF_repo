#include "lattice.h"

#include <TCanvas.h>
#include <TAxis.h>
#include <TFrame.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TNtuple.h>
#include <TGraph.h>
#include <TMultiGraph.h>

void snapshot(const std::string out_image = "snapshot.png") {
	std::stringstream title;
	title << "beta = " << beta;
	auto c1 = new TCanvas("c1", out_image.c_str(), 1000, 1000);
	TGraph* graph[nparticles];
	TMultiGraph* multigraph = new TMultiGraph();
	TLegend* legend = new TLegend(0.75, 0.75, 0.85, 0.85);
	int time[p_length];
	for(int ii = 0; ii < p_length; ii++) time[ii] = ii;
	for(int i = 0; i < nparticles; i++) {
		graph[i] = new TGraph(p_length, time, y + p_length * i);
		graph[i]->SetMarkerColor(i + 1);
		multigraph->Add(graph[i]);
		std::ostringstream string_stream2;
		string_stream2 << "particle " << i;
		std::string legend_entry = string_stream2.str();
		legend->AddEntry(graph[i], legend_entry.c_str(), "p");
	}
	c1->SetGrid();
	c1->SaveAs(out_image.c_str());
	delete c1;
	for(int i = 0; i < nparticles; i++) delete graph[i];
	delete multigraph;
}