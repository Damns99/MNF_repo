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
	
	std::string filename = "HH_max_gaussian_fitres.txt";
	
	std::vector<double> ampl, stdev, m, dm, q, dq, c, dc;
	textIo::textIn(filename, '\t', '#', &ampl, &stdev, &m, &dm, &q, &dq, &c, &dc);
	ampl.pop_back(); stdev.pop_back(); m.pop_back(); dm.pop_back(); q.pop_back(); dq.pop_back(); c.pop_back(); dc.pop_back();
	int N = ampl.size();
	
	{
	TCanvas* canvas = new TCanvas("canvas", "Canvas", 600, 400);
	TGraphErrors* graph = new TGraphErrors(N, ampl.data(), m.data(), nullptr, dm.data());
	graph->Sort();
	graph->SetTitle(";ampl [mV];m [s/dam]");
	canvas->cd();
	canvas->SetGrid();
	// canvas->SetLogx();
	// canvas->SetLogy();
	graph->Draw("AP");
	
	canvas->SaveAs("HH_test_max_fitres_ampl_m.pdf"); 
	}
	
	{
	TCanvas* canvas = new TCanvas("canvas", "Canvas", 600, 400);
	TGraphErrors* graph = new TGraphErrors(N, ampl.data(), q.data(), nullptr, dq.data());
	graph->Sort();
	graph->SetTitle(";ampl [mV];q [ms]");
	canvas->cd();
	canvas->SetGrid();
	// canvas->SetLogx();
	// canvas->SetLogy();
	graph->Draw("AP");
	
	canvas->SaveAs("HH_test_max_fitres_ampl_q.pdf"); 
	}
	
	{
	TCanvas* canvas = new TCanvas("canvas", "Canvas", 600, 400);
	TGraphErrors* graph = new TGraphErrors(N, ampl.data(), c.data(), nullptr, dc.data());
	graph->Sort();
	graph->SetTitle(";ampl [mV];c [ms]");
	canvas->cd();
	canvas->SetGrid();
	// canvas->SetLogx();
	// canvas->SetLogy();
	graph->Draw("AP");
	
	canvas->SaveAs("HH_test_max_fitres_ampl_c.pdf"); 
	}
	
	{
	TCanvas* canvas = new TCanvas("canvas", "Canvas", 600, 400);
	TGraphErrors* graph = new TGraphErrors(N, stdev.data(), m.data(), nullptr, dm.data());
	graph->Sort();
	graph->SetTitle(";stdev [mV];m [s/dam]");
	canvas->cd();
	canvas->SetGrid();
	// canvas->SetLogx();
	// canvas->SetLogy();
	graph->Draw("AP");
	
	canvas->SaveAs("HH_test_max_fitres_stdev_m.pdf"); 
	}
	
	{
	TCanvas* canvas = new TCanvas("canvas", "Canvas", 600, 400);
	TGraphErrors* graph = new TGraphErrors(N, stdev.data(), q.data(), nullptr, dq.data());
	graph->Sort();
	graph->SetTitle(";stdev [mV];q [ms]");
	canvas->cd();
	canvas->SetGrid();
	// canvas->SetLogx();
	// canvas->SetLogy();
	graph->Draw("AP");
	
	canvas->SaveAs("HH_test_max_fitres_stdev_q.pdf"); 
	}
	
	{
	TCanvas* canvas = new TCanvas("canvas", "Canvas", 600, 400);
	TGraphErrors* graph = new TGraphErrors(N, stdev.data(), c.data(), nullptr, dc.data());
	graph->Sort();
	graph->SetTitle(";stdev [mV];c [ms]");
	canvas->cd();
	canvas->SetGrid();
	// canvas->SetLogx();
	// canvas->SetLogy();
	graph->Draw("AP");
	
	canvas->SaveAs("HH_test_max_fitres_stdev_c.pdf"); 
	}
}