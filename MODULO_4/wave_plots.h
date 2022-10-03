#ifndef WAVE_PLOTS_H
#define WAVE_PLOTS_H

#include <string>
#include <vector>
#include <sstream>

#include "bound_cond_vecs.h"
#include "myfft.h"

#include <TH2D.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TFrame.h>
#include <TH1.h>
#include <TLegend.h>

#include <TError.h>

#define SURF_PLOT 101
#define CONT_PLOT 102
#define COLZ_PLOT 103

double bound(double a, double b, double c) {
	if(a > b && a < c) return a;
	else if(a >= c) return c;
	else return b;
}

namespace waveplots {
	
	using DblBcv = bound_cond_vecs::BoundCondVec<double>;
	
	void plot(const std::vector<std::vector<DblBcv>> u, double t0, double dt, int nt, double x0, double dx, int nx, const std::string file, int mode, double* mn, double* mx) {
		gErrorIgnoreLevel = kWarning;
		
		for(int k = 0; k < u[0].size(); k++) {
			std::stringstream ss;
			ss << "_" << k;
			std::string filename = file + ss.str();
			
			auto canvas = new TCanvas("canvas", filename.c_str(), 600, 400);
			auto hist = new TH2D("hist", filename.c_str(), nt, t0, t0 + nt * dt, nx, x0, x0 + nx * dx);
			for(int i = 0; i < nt; i++) for(int j = 0; j < nx; j++) hist->Fill(t0 + (i + 0.5) * dt, x0 + (j + 0.5) * dx, bound(u[i][k][j], mn[k], mx[k]));
			
			hist->SetBit(TH1::kNoStats);
			hist->SetTitle((filename + ";t;x;u(x,t)").c_str());
			hist->SetMaximum(mx[k]);
			hist->SetMinimum(mn[k]);
			
			if(mode == SURF_PLOT) hist->Draw("SURF");
			else if(mode == CONT_PLOT) hist->Draw("CONT1Z");
			else hist->Draw("COLZ0");
			
			canvas->SaveAs((filename + ".pdf").c_str());
			canvas->SaveSource((filename + ".cpp").c_str());
			
			delete canvas;
			delete hist;
			
			TCanvas* canvas2 = new TCanvas("canvas2", "Canvas2", 600, 400);
			TMultiGraph* multigraph = new TMultiGraph();
			TLegend* legend = new TLegend(0.75, 0.75, 0.85, 0.85);
			
			double tmpx[nx], first[nx], last[nx];
			for(int j = 0; j < nx; j++) {
				tmpx[j] = x0 + j * dx;
				first[j] = bound(u[0][k][j], mn[k], mx[k]);
				last[j] = bound(u[nt - 1][k][j], mn[k], mx[k]);
			}
			
			TGraph* graph1 = new TGraph(nx, tmpx, first);
			graph1->SetMarkerColor(1);
			graph1->SetLineColor(1);
			multigraph->Add(graph1);
			legend->AddEntry(graph1, "t = 0.", "p");
			TGraph* graph2 = new TGraph(nx, tmpx, last);
			graph2->SetMarkerColor(2);
			graph2->SetLineColor(2);
			multigraph->Add(graph2);
			std::stringstream ss2;
			ss2 << "t = " << (nt - 1) * dt;
			legend->AddEntry(graph2, ss2.str().c_str(), "p");
			
			multigraph->SetTitle((filename + ";x;u(x,t)").c_str());
			canvas2->cd();
			canvas2->SetGrid();
			multigraph->Draw("APL");
			legend->Draw();
			canvas2->SaveAs((filename + "_fl.pdf").c_str());
			canvas2->SaveSource((filename + "_fl.cpp").c_str());
			
			delete canvas2;
			delete multigraph;
		}
	}
	
	void plotFFT(const std::vector<std::vector<DblBcv>> u, double t0, double dt, int nt, const std::string file, int mode) {
		gErrorIgnoreLevel = kWarning;
		
		double x0 = 0., dx = 1.;
		int nx = u[0][0].len();
		
		for(int k = 0; k < u[0].size(); k++) {
			std::stringstream ss;
			ss << "_" << k;
			std::string filename = file + ss.str();
			
			auto canvas = new TCanvas("canvas", filename.c_str(), 600, 400);
			auto hist = new TH2D("hist", filename.c_str(), nt, t0, t0 + nt * dt, nx, x0, x0 + nx * dx);
			for(int i = 0; i < nt; i++) {
				double real[nx], imag[nx];
				myFFT::my_dft_r2c_1d(nx, u[i][k].data(), real, imag);
				for(int j = 1; j < nx / 2 + 1; j++) hist->Fill(t0 + (i + 0.5) * dt, x0 + (j + 0.5) * dx, real[j] * real[j] + imag[j] * imag[j]);
			}
			
			hist->SetBit(TH1::kNoStats);
			hist->SetTitle((filename + ";t;k;|u(k,t)|^2").c_str());
			canvas->SetLogz();
			canvas->SetLogy();
			
			if(mode == SURF_PLOT) hist->Draw("SURF");
			else if(mode == CONT_PLOT) hist->Draw("CONT1Z");
			else hist->Draw("COLZ");
			
			canvas->SaveAs((filename + ".pdf").c_str());
			canvas->SaveSource((filename + ".cpp").c_str());
			
			delete canvas;
			delete hist;
			
			TCanvas* canvas2 = new TCanvas("canvas2", "Canvas2", 600, 400);
			TMultiGraph* multigraph = new TMultiGraph();
			TLegend* legend = new TLegend(0.75, 0.75, 0.85, 0.85);
			
			double real0[nx], imag0[nx];
			myFFT::my_dft_r2c_1d(nx, u[0][k].data(), real0, imag0);
			double real1[nx], imag1[nx];
			myFFT::my_dft_r2c_1d(nx, u[nt - 1][k].data(), real1, imag1);
			double tmpx[nx / 2], first[nx / 2], last[nx / 2];
			for(int j = 1; j < nx / 2 + 1; j++) {
				tmpx[j - 1] = x0 + j * dx;
				first[j - 1] = real0[j] * real0[j] + imag0[j] * imag0[j];
				last[j - 1] = real1[j] * real1[j] + imag1[j] * imag1[j];
			}
			
			TGraph* graph1 = new TGraph(nx / 2, tmpx, first);
			graph1->SetMarkerColor(1);
			graph1->SetLineColor(1);
			multigraph->Add(graph1);
			legend->AddEntry(graph1, "t = 0.", "p");
			TGraph* graph2 = new TGraph(nx / 2, tmpx, last);
			graph2->SetMarkerColor(2);
			graph2->SetLineColor(2);
			multigraph->Add(graph2);
			std::stringstream ss2;
			ss2 << "t = " << (nt - 1) * dt;
			legend->AddEntry(graph2, ss2.str().c_str(), "p");
			
			multigraph->SetTitle((filename + ";k;|u(k,t)|^2").c_str());
			canvas2->SetLogx();
			canvas2->SetLogy();
			canvas2->cd();
			canvas2->SetGrid();
			multigraph->Draw("APL");
			legend->Draw();
			canvas2->SaveAs((filename + "_fl.pdf").c_str());
			canvas2->SaveSource((filename + "_fl.cpp").c_str());
			
			delete canvas2;
			delete multigraph;
		}
	}
	
}

#endif