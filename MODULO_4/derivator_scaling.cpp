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
		bound_cond_vecs::BoundCondVec<double> tmpx(N), x(N), dx0(N);
		for(int ii = 0; ii < N; ii++) {
			tmpx[ii] = 2. * M_PI * ii * dx;
			x[ii] = cos(tmpx[ii]) + sin(2. * tmpx[ii]) - sin(3. * tmpx[ii] + 0.1);
			dx0[ii] = 2. * M_PI * (- sin(tmpx[ii]) + 2. * cos(2. * tmpx[ii]) - 3. * cos(3. * tmpx[ii] + 0.1));
		}
		
		bound_cond_vecs::BoundCondVec<double> dx1 = derivators::fwd_derive(x, dx);
		double mse1 = mse(dx0, dx1);
		bound_cond_vecs::BoundCondVec<double> dx2 = derivators::symm_derive(x, dx);
		double mse2 = mse(dx0, dx2);
		bound_cond_vecs::BoundCondVec<double> dx3 = derivators::fft_derive(x, dx);
		double mse3 = mse(dx0, dx3);
		bound_cond_vecs::BoundCondVec<double> dx4 = derivators::fwd_derive_i(x, dx);
		double mse4 = mse(dx0, dx4);
		
		std::cout << std::setprecision(3);
		std::cout << "function: cos(2 pi x)" << std::endl;
		std::cout << "N = " << N << "   dx = " << dx << "   L = " << N * dx << std::endl;
		std::cout << "Forward derivative: " << '\t' << mse1 << std::endl;
		std::cout << "Symmetric derivative: " << '\t' << mse2 << std::endl;
		std::cout << "FFT derivative: " << '\t' << mse3 << std::endl;
		std::cout << "Forward derivative i: " << '\t' << mse4 << std::endl;
		
		TCanvas* canvas = new TCanvas("canvas0", "Canvas0", 600, 400);
		TMultiGraph* multigraph = new TMultiGraph();
		TLegend* legend = new TLegend(0.75, 0.75, 0.85, 0.85);
		TGraph* graph[4];
		
		TGraph* oggraph = new TGraph(N, tmpx.data(), x.data());
		oggraph->SetMarkerColor(1);
		multigraph->Add(oggraph);
		legend->AddEntry(oggraph, "origin", "p");
		TGraph* thgraph = new TGraph(N, tmpx.data(), dx0.data());
		thgraph->SetMarkerColor(2);
		multigraph->Add(thgraph);
		legend->AddEntry(thgraph, "theory", "p");
		graph[0] = new TGraph(N, tmpx.data(), dx1.data());
		graph[0]->SetMarkerColor(1);
		multigraph->Add(graph[0]);
		legend->AddEntry(graph[0], "forward", "p");
		graph[1] = new TGraph(N, tmpx.data(), dx2.data());
		graph[1]->SetMarkerColor(2);
		multigraph->Add(graph[1]);
		legend->AddEntry(graph[1], "symmetric", "p");
		graph[2] = new TGraph(N, tmpx.data(), dx3.data());
		graph[2]->SetMarkerColor(3);
		multigraph->Add(graph[2]);
		legend->AddEntry(graph[2], "fft", "p");
		graph[3] = new TGraph(N, tmpx.data(), dx4.data());
		graph[3]->SetMarkerColor(4);
		multigraph->Add(graph[3]);
		legend->AddEntry(graph[3], "forward_i", "p");
		
		canvas->cd();
		canvas->SetGrid();
		multigraph->Draw("AP");
		legend->Draw();
		canvas->SaveAs("derivator_scaling_plot.pdf");
	}
	{
		int ntries = 12, N = 16;
		double dx = 1. / N;
		int dN = 2;
		double resx[ntries], resy[4][ntries];
		for(int jj = 0; jj < ntries; jj++) {
			bound_cond_vecs::BoundCondVec<double> tmpx(N), x(N), dx0(N);
			for(int ii = 0; ii < N; ii++) {
				tmpx[ii] = 2. * M_PI * ii * dx;
				x[ii] = cos(tmpx[ii]) + sin(2. * tmpx[ii]) - sin(3. * tmpx[ii] + 0.1);
				dx0[ii] = 2. * M_PI * (- sin(tmpx[ii]) + 2. * cos(2. * tmpx[ii]) - 3. * cos(3. * tmpx[ii] + 0.1));
			}
			resx[jj] = dx;
			bound_cond_vecs::BoundCondVec<double> dx1 = derivators::fwd_derive(x, dx);
			resy[0][jj] = mse(dx0, dx1);
			bound_cond_vecs::BoundCondVec<double> dx2 = derivators::symm_derive(x, dx);
			resy[1][jj] = mse(dx0, dx2);
			bound_cond_vecs::BoundCondVec<double> dx3 = derivators::fft_derive(x, dx);
			resy[2][jj] = mse(dx0, dx3);
			bound_cond_vecs::BoundCondVec<double> dx4 = derivators::fwd_derive_i(x, dx);
			resy[3][jj] = mse(dx0, dx4);
			N *= dN;
			dx = 1. / N;
		}
		TCanvas* canvas = new TCanvas("canvas1", "Canvas1", 600, 400);
		TMultiGraph* multigraph = new TMultiGraph();
		TLegend* legend = new TLegend(0.75, 0.25, 0.90, 0.40);
		TGraph* graph[4];
		
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
		graph[3] = new TGraph(ntries, resx, resy[3]);
		graph[3]->SetMarkerColor(4);
		multigraph->Add(graph[3]);
		legend->AddEntry(graph[3], "forward_i", "p");
		
		canvas->cd();
		canvas->SetGrid();
		canvas->SetLogx();
		canvas->SetLogy();
		multigraph->Draw("A*");
		legend->Draw();
		canvas->SaveAs("derivator_scaling_dx.pdf");
	}
	
	{
		int N = 1000;
		double dx = 1. / N;
		bound_cond_vecs::BoundCondVec<double> tmpx(N), x(N), dx0(N);
		for(int ii = 0; ii < N; ii++) {
			tmpx[ii] = 2. * M_PI * ii * dx;
			x[ii] = cos(tmpx[ii]) + sin(2. * tmpx[ii]) - sin(3. * tmpx[ii] + 0.1);
			dx0[ii] = 4. * M_PI * M_PI * (- cos(tmpx[ii]) - 4. * sin(2. * tmpx[ii]) + 9. * sin(3. * tmpx[ii] + 0.1));
		}
		
		bound_cond_vecs::BoundCondVec<double> dx1 = derivators::fwd_derive_2(x, dx);
		double mse1 = mse(dx0, dx1);
		bound_cond_vecs::BoundCondVec<double> dx2 = derivators::symm_derive_2(x, dx);
		double mse2 = mse(dx0, dx2);
		bound_cond_vecs::BoundCondVec<double> dx3 = derivators::fft_derive_2(x, dx);
		double mse3 = mse(dx0, dx3);
		
		std::cout << std::setprecision(3);
		std::cout << "function: cos(2 pi x)" << std::endl;
		std::cout << "N = " << N << "   dx = " << dx << "   L = " << N * dx << std::endl;
		std::cout << "Forward derivative 2: " << '\t' << mse1 << std::endl;
		std::cout << "Symmetric derivative 2: " << '\t' << mse2 << std::endl;
		std::cout << "FFT derivative 2: " << '\t' << mse3 << std::endl;
		
		TCanvas* canvas = new TCanvas("canvas2", "Canvas2", 600, 400);
		TMultiGraph* multigraph = new TMultiGraph();
		TLegend* legend = new TLegend(0.75, 0.75, 0.85, 0.85);
		TGraph* graph[3];
		
		TGraph* oggraph = new TGraph(N, tmpx.data(), x.data());
		oggraph->SetMarkerColor(4);
		multigraph->Add(oggraph);
		legend->AddEntry(oggraph, "origin", "p");
		TGraph* thgraph = new TGraph(N, tmpx.data(), dx0.data());
		thgraph->SetMarkerColor(5);
		multigraph->Add(thgraph);
		legend->AddEntry(thgraph, "theory", "p");
		graph[0] = new TGraph(N, tmpx.data(), dx1.data());
		graph[0]->SetMarkerColor(1);
		multigraph->Add(graph[0]);
		legend->AddEntry(graph[0], "forward", "p");
		graph[1] = new TGraph(N, tmpx.data(), dx2.data());
		graph[1]->SetMarkerColor(2);
		multigraph->Add(graph[1]);
		legend->AddEntry(graph[1], "symmetric", "p");
		graph[2] = new TGraph(N, tmpx.data(), dx3.data());
		graph[2]->SetMarkerColor(3);
		multigraph->Add(graph[2]);
		legend->AddEntry(graph[2], "fft", "p");
		
		canvas->cd();
		canvas->SetGrid();
		multigraph->Draw("AP");
		legend->Draw();
		canvas->SaveAs("derivator_scaling_plot_2.pdf");
	}
	{
		int ntries = 12, N = 16;
		double dx = 1. / N;
		int dN = 2;
		double resx[ntries], resy[3][ntries];
		for(int jj = 0; jj < ntries; jj++) {
			bound_cond_vecs::BoundCondVec<double> tmpx(N), x(N), dx0(N);
			for(int ii = 0; ii < N; ii++) {
				tmpx[ii] = 2. * M_PI * ii * dx;
				x[ii] = cos(tmpx[ii]) + sin(2. * tmpx[ii]) - sin(3. * tmpx[ii] + 0.1);
				dx0[ii] = 4. * M_PI * M_PI * (- cos(tmpx[ii]) - 4. * sin(2. * tmpx[ii]) + 9. * sin(3. * tmpx[ii] + 0.1));
			}
			resx[jj] = dx;
			bound_cond_vecs::BoundCondVec<double> dx1 = derivators::fwd_derive_2(x, dx);
			resy[0][jj] = mse(dx0, dx1);
			bound_cond_vecs::BoundCondVec<double> dx2 = derivators::symm_derive_2(x, dx);
			resy[1][jj] = mse(dx0, dx2);
			bound_cond_vecs::BoundCondVec<double> dx3 = derivators::fft_derive_2(x, dx);
			resy[2][jj] = mse(dx0, dx3);
			N *= dN;
			dx = 1. / N;
		}
		TCanvas* canvas = new TCanvas("canvas3", "Canvas3", 600, 400);
		TMultiGraph* multigraph = new TMultiGraph();
		TLegend* legend = new TLegend(0.75, 0.25, 0.90, 0.40);
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
		canvas->SaveAs("derivator_scaling_dx_2.pdf");
	}
}