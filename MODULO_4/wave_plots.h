#ifndef WAVE_PLOTS_H
#define WAVE_PLOTS_H

#include <string>
#include <vector>
#include <sstream>

#include "bound_cond_vecs.h"
#include "myfft.h"

#include <TH2D.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TFrame.h>
#include <TH1.h>

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
	
	void plot(const std::vector<std::vector<DblBcv>> u, double t0, double dt, int nt, double x0, double dx, int nx, const std::string file, int mode, double mn = -10., double mx = 10.) {
		gErrorIgnoreLevel = kWarning;
		
		for(int k = 0; k < u[0].size(); k++) {
			std::stringstream ss;
			ss << "_" << k;
			std::string filename = file + ss.str();
			
			auto canvas = new TCanvas("canvas", filename.c_str(), 600, 400);
			auto hist = new TH2D("hist", filename.c_str(), nt, t0, t0 + nt * dt, nx, x0, x0 + nx * dx);
			for(int i = 0; i < nt; i++) for(int j = 0; j < nx; j++) hist->Fill(t0 + (i + 0.5) * dt, x0 + (j + 0.5) * dx, bound(u[i][k][j], mn, mx));
			
			hist->SetBit(TH1::kNoStats);
			hist->SetTitle((filename + ";t;x;u(x,t)").c_str());
			hist->SetMaximum(mx);
			hist->SetMinimum(mn);
			
			if(mode == SURF_PLOT) hist->Draw("SURF");
			else if(mode == CONT_PLOT) hist->Draw("CONT1Z");
			else hist->Draw("COLZ0");
			
			canvas->SaveAs((filename + ".pdf").c_str());
			canvas->SaveSource((filename + ".cpp").c_str());
			
			delete canvas;
			delete hist;
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
			hist->SetTitle((filename + ";t;k;|u_k(t)|^2").c_str());
			canvas->SetLogz();
			canvas->SetLogy();
			
			if(mode == SURF_PLOT) hist->Draw("SURF");
			else if(mode == CONT_PLOT) hist->Draw("CONT1Z");
			else hist->Draw("COLZ");
			
			canvas->SaveAs((filename + ".pdf").c_str());
			canvas->SaveSource((filename + ".cpp").c_str());
			
			delete canvas;
			delete hist;
		}
	}
	
}

#endif