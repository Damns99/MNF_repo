#include "lattice.h"

#include <TCanvas.h>
#include <TAxis.h>
#include <TFrame.h>
#include <TH2I.h>
#include <TH1.h>

void Lattice2D::snapshot(const std::string out_image = "snapshot.png") {
	std::stringstream title;
	title << "beta = " << beta << "    ener = " << energy << "    magn = " << magnetization;
	auto c1 = new TCanvas("c1", out_image.c_str(), 1000, 1000);
	auto h1 = new TH2I("h1", title.str().c_str(), length, 0, length, length, 0, length);
	for(int i = 0; i < length; i++) for(int j = 0; j < length; j++) h1->Fill(i, j, spin[index2d(i, j, length)]);
	h1->SetBit(TH1::kNoStats);
	h1->SetMinimum(-1);
	h1->SetMaximum(+1);
	h1->Draw("COL");
	c1->SaveAs(out_image.c_str());
	delete c1;
	delete h1;
}