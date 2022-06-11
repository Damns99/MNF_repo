#include <iostream>
#include <string>
#include <math.h>
#include <iomanip>

#include "bound_cond_vecs.h"
#include "derivators.h"

#include <TGraph.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TFrame.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TNtuple.h>

double mse(const bound_cond_vecs::BoundCondVec<double>& a, const bound_cond_vecs::BoundCondVec<double>& b) {
	double result = 0.;
	for(int ii = 0; ii < a.len(); ii++) result += (a[ii] - b[ii]) * (a[ii] - b[ii]) / a.len();
	return result;
}

int main() {
	{
		int N = 1000;
		double dx = 1. / N;
		bound_cond_vecs::BoundCondVec<double> x(N), dx0(N);
		for(int ii = 0; ii < N; ii++) {
			double tmpx = 2. * M_PI * ii * dx;
			x[ii] = cos(tmpx);
			dx0[ii] = - 2. * M_PI * sin(tmpx);
		}
		
		bound_cond_vecs::BoundCondVec<double> dx1 = derivators::fwd_derive(x, dx);
		double mse1 = mse(dx0, dx1);
		bound_cond_vecs::BoundCondVec<double> dx2 = derivators::symm_derive(x, dx);
		double mse2 = mse(dx0, dx2);
		bound_cond_vecs::BoundCondVec<double> dx3 = derivators::fft_derive(x, dx);
		double mse3 = mse(dx0, dx3);
		
		std::cout << std::setprecision(3);
		std::cout << "function: cos(2 pi x)" << std::endl;
		std::cout << "N = " << N << "   dx = " << dx << "   L = " << N * dx << std::endl;
		std::cout << "Forward derivative: " << '\t' << mse1 << std::endl;
		std::cout << "Symmetric derivative: " << '\t' << mse2 << std::endl;
		std::cout << "FFT derivative: " << '\t' << mse3 << std::endl;
	}
	{
		int ntries = 20, N = 1000;
		double dx = 0.0005;
		double ddx = 0.75;
		bound_cond_vecs::BoundCondVec<double> x(N), dx0(N);
		double resx[ntries], resy[3][ntries];
		for(int jj = 0; jj < ntries; jj++) {
			for(int ii = 0; ii < N; ii++) {
				double tmpx = 2. * M_PI * ii * dx;
				x[ii] = cos(tmpx);
				dx0[ii] = - 2. * M_PI * sin(tmpx);
			}
			resx[jj] = dx;
			bound_cond_vecs::BoundCondVec<double> dx1 = derivators::fwd_derive(x, dx);
			resy[0][jj] = mse(dx0, dx1);
			bound_cond_vecs::BoundCondVec<double> dx2 = derivators::symm_derive(x, dx);
			resy[1][jj] = mse(dx0, dx2);
			bound_cond_vecs::BoundCondVec<double> dx3 = derivators::fft_derive(x, dx);
			resy[2][jj] = mse(dx0, dx3);
			dx = double(dx * ddx);
		}
		TCanvas* canvas = new TCanvas("canvas", "Canvas", 600, 400);
		TMultiGraph* multigraph = new TMultiGraph();
		TLegend* legend = new TLegend(0.75, 0.75, 0.85, 0.85);
		TGraph* graph[3];
		
		graph[0] = new TGraph(ntries, resx, resy[0]);
		graph[0]->SetMarkerColor(1);
		multigraph->Add(graph[0]);
		legend->AddEntry(graph[0], "forward", "p");
		graph[1] = new TGraph(ntries, resx, resy[1]);
		graph[1]->SetMarkerColor(2);
		multigraph->Add(graph[1]);
		legend->AddEntry(graph[1], "symmetric", "p");
		graph[2] = new TGraph(ntries, resx, resy[2]);
		graph[2]->SetMarkerColor(3);
		multigraph->Add(graph[2]);
		legend->AddEntry(graph[2], "fft", "p");
		
		canvas->cd();
		canvas->SetGrid();
		canvas->SetLogx();
		canvas->SetLogy();
		multigraph->Draw("A*");
		legend->Draw();
		canvas->SaveAs("derivator_scaling_dx.pdf");
	}
	{
		int ntries = 20, N = 64;
		double dx = 0.0001;
		double dN = 1.5;
		bound_cond_vecs::BoundCondVec<double> x(N), dx0(N);
		double resx[ntries], resy[3][ntries];
		for(int jj = 0; jj < ntries; jj++) {
			for(int ii = 0; ii < N; ii++) {
				double tmpx = 2. * M_PI * ii * dx;
				x[ii] = cos(tmpx);
				dx0[ii] = - 2. * M_PI * sin(tmpx);
			}
			resx[jj] = N;
			bound_cond_vecs::BoundCondVec<double> dx1 = derivators::fwd_derive(x, dx);
			resy[0][jj] = mse(dx0, dx1);
			bound_cond_vecs::BoundCondVec<double> dx2 = derivators::symm_derive(x, dx);
			resy[1][jj] = mse(dx0, dx2);
			bound_cond_vecs::BoundCondVec<double> dx3 = derivators::fft_derive(x, dx);
			resy[2][jj] = mse(dx0, dx3);
			N = int(N * dN);
		}
		TCanvas* canvas = new TCanvas("canvas1", "Canvas1", 600, 400);
		TMultiGraph* multigraph = new TMultiGraph();
		TLegend* legend = new TLegend(0.75, 0.75, 0.85, 0.85);
		TGraph* graph[3];
		
		graph[0] = new TGraph(ntries, resx, resy[0]);
		graph[0]->SetMarkerColor(1);
		multigraph->Add(graph[0]);
		legend->AddEntry(graph[0], "forward", "p");
		graph[1] = new TGraph(ntries, resx, resy[1]);
		graph[1]->SetMarkerColor(2);
		multigraph->Add(graph[1]);
		legend->AddEntry(graph[1], "symmetric", "p");
		graph[2] = new TGraph(ntries, resx, resy[2]);
		graph[2]->SetMarkerColor(3);
		multigraph->Add(graph[2]);
		legend->AddEntry(graph[2], "fft", "p");
		
		canvas->cd();
		canvas->SetGrid();
		canvas->SetLogx();
		canvas->SetLogy();
		multigraph->Draw("A*");
		legend->Draw();
		canvas->SaveAs("derivator_scaling_N.pdf");
	}
}