#include "lattice.h"

#include <TCanvas.h>
#include <TAxis.h>
#include <TFrame.h>
#include <TH2Poly.h>
#include <TH1.h>

void snapshot(const std::string out_image = "snapshot.png") {
	std::stringstream title;
	title << "beta = " << beta << "    ener = " << energy << "    magn = " << magnetization;
	auto c1 = new TCanvas("c1", out_image.c_str(), 1000, 1000);
	TH2Poly *h1 = new TH2Poly();
	h1->SetName("h1");
	h1->SetTitle(title.str().c_str());
	for(int i = 0; i < length; i++) for(int j = 0; j < length; j++) {
		int nt = node_type(i, j);
		switch(nt) {
			case 0:
				break;
			case 1:
				{
				double x[] = {j+1., j-0.5, j-1., j+0.5};
				double y[] = {i+1., i+0.5, i-1., i-0.5};
				h1->AddBin(4, x, y);
				}
				break;
			case 2:
				{
				double x[] = {j+1., j+0.5, j-1., j-0.5};
				double y[] = {i+0., i+0.5, i+0., i-0.5};
				h1->AddBin(4, x, y);
				}
				break;
			case 3:
				{
				double x[] = {j+0.5, j+0., j-0.5, j+0.};
				double y[] = {i+0.5, i+1., i-0.5, i-1.};
				h1->AddBin(4, x, y);
				}
				break;
		}
	}
	for(int i = 0; i < length; i++) for(int j = 0; j < length; j++) if(node_type(i, j) != 0) h1->Fill(j+1./4., i, spin[index2d(i, j, length)]);
	h1->SetBit(TH1::kNoStats);
	h1->SetMinimum(-1);
	h1->SetMaximum(+1);
	h1->Draw("COL");
	c1->SaveAs(out_image.c_str());
	delete c1;
	delete h1;
}