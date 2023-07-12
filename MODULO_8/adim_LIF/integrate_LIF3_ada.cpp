#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <math.h>
#include <chrono>
#include <ctime>
#include <filesystem>
namespace fs = std::filesystem;

#include <TGraph.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TFrame.h>
#include <TLegend.h>

#include "mystyle.h"
#include "integrate_LIF.h"
#include "text_io.h"

void addToMultigraph(TMultiGraph* multigraph, TLegend* legend, std::vector<double>& x, std::vector<double>& y, int n, int color, const char* name, const char* labeltype = "p", const char* drawoption = "") {
	TGraph* graph = new TGraph(n, x.data(), y.data());
	graph->SetMarkerColor(color);
	graph->SetLineColor(color);
	multigraph->Add(graph, drawoption);
	legend->AddEntry(graph, name, labeltype);
}

void addToMultigraph(TMultiGraph* multigraph, TLegend* legend, std::vector<int>& x, std::vector<double>& y, int n, double simtime, int color, const char* name, const char* labeltype = "p", const char* drawoption = "") {
	std::vector<double> xx(n);
	for(int ii = 0; ii < n; ii++) xx[ii] = simtime/(x[ii]-1);
	TGraph* graph = new TGraph(n, xx.data(), y.data());
	graph->SetMarkerColor(color);
	graph->SetLineColor(color);
	multigraph->Add(graph, drawoption);
	legend->AddEntry(graph, name, labeltype);
}

/* double exact_sol(double x, double offsetz, double periodz, double az) {
	return 2. - ((1. - offsetz) * (1 - exp(-x)) + exp(-x) / (4*M_PI*M_PI/periodz/periodz + 1) * az * (-exp(x)*sin(2*M_PI/periodz*x) + exp(x)*2*M_PI/periodz*cos(2*M_PI/periodz*x) - 2*M_PI/periodz) + 1.*exp(-x)); 
} */
/* double exact_sol(double x, double startz, double durz, double az) {
	if(x<startz) return 1.;
	//else if(x<=startz+durz) return ((az+1)*exp(-(x-startz))*(exp((x-startz))+1) + 1.*exp(-(x-startz)));
	else if(x<startz+durz) return 1. + az - az*exp(-(x-startz));
	else return (exact_sol((startz+durz-1e-10),startz,durz,az)-1)*exp(-(x-startz-durz)) + 1;
} */

double func1(double x, double simtime, double az) {
	if(x<simtime/4) return az*x;
	else if(x<simtime*3/4) return az*(simtime/2-x);
	else return az*(-simtime+x);
}
double exact_sol(double x, double simtime, double az) {
	if(x<simtime/4) return (az-1)*exp(-x)+az*(x-1)+1+1*exp(-x);
	else if(x<simtime*3/4) return exact_sol(simtime/4-1e-10,simtime,az)*exp(simtime/4-x) - (az*simtime/4+az+1)*exp(simtime/4-x) + az*(2*simtime/4-x+1) +1;
	else return exact_sol(simtime*3/4-1e-10,simtime,az)*exp(simtime*3/4-x) + (az*simtime/4+az-1)*exp(simtime*3/4-x) + az*(-simtime+x-1) +1;
}

int main() {
	// Voltages in [Vrest = -70e-3 V], Currents in [g*Vrest = -7e-10 A], Times in [Cm/g = 1e-2 s]
	double params[2] = {
		0.7, // yth
		1.,  // yreset
	};
	double yth = params[0], yreset = params[1];
	double simtime = 1.e1;
	int N0 = 1000;
	double h0 = simtime/(N0-1);
	double x0 = 0.;
	double y0 = 1.;
	
	// Currents
	double az = -1;
	// pulse_train
	double startz = 2, durz = 5, intz = 20;
	int numz = 25;
	// sine_wave
	double periodz = 1., phasez = 0., offsetz = -0.5;
	/* double meanz = -0.10, variancez = 0.0001;*/
	
	double epsilon = 1e-5;
	double h;
	
	// Heun
	std::vector<double> x1, y1;
	x1.push_back(x0); y1.push_back(y0);
	int n1 = 0;
	h = h0;
	while(x1[n1] < simtime) {
		int iter = 0;
		while(true) {
			double curr_0 = az*sin(2*M_PI*x1[n1]/periodz+phasez) + offsetz;
			double curr_2 = az*sin(2*M_PI*(x1[n1]+h)/periodz+phasez) + offsetz;
			double yh;
			double a = h - h*h/2;
			yh = (1-a) * y1[n1] + (a-h/2) * curr_0 + (h/2) * curr_2 + a;
			if(yh <= yth) yh = yreset;
			
			double h2 = h/2;
			double curr_1 = az*sin(2*M_PI*(x1[n1]+h2)/periodz+phasez) + offsetz;
			double yh2;
			double a2 = h2 - h2*h2/2;
			yh2 = (1-a2) * y1[n1] + (a2-h2/2) * curr_0 + (h2/2) * curr_1 + a2;
			if(yh2 <= yth) yh2 = yreset;
			yh2 = (1-a2) * yh2 + (a2-h2/2) * curr_1 + (h2/2) * curr_2 + a2;
			if(yh2 <= yth) yh2 = yreset;
			
			double trunc = abs(1./3. * (yh2 - yh));
			if((0.5*epsilon*h/2 > trunc || 2.*epsilon*h/2 < trunc) && trunc > 0. && iter < 100) {
				std::cout << "changing h from " << h << " to " << h*pow(epsilon*h/2/trunc, 1./3.) << " at step " << n1 << std::endl;
				h = h*pow(epsilon*h/2/trunc, 1./3.);
				iter++;
			}
			else {
				y1.push_back(yh2);
				x1.push_back(x1[n1]+h);
				break;
			}
		}
		n1++;
	}
	
	std::cout << std::endl;
	
	epsilon = 1e-6;
	
	// Heun
	std::vector<double> x11, y11;
	x11.push_back(x0); y11.push_back(y0);
	int n11 = 0;
	h = h0;
	while(x11[n11] < simtime) {
		int iter = 0;
		while(true) {
			double curr_0 = az*sin(2*M_PI*x11[n11]/periodz+phasez) + offsetz;
			double curr_2 = az*sin(2*M_PI*(x11[n11]+h)/periodz+phasez) + offsetz;
			double yh;
			double a = h - h*h/2;
			yh = (1-a) * y11[n11] + (a-h/2) * curr_0 + (h/2) * curr_2 + a;
			if(yh <= yth) yh = yreset;
			
			double h2 = h/2;
			double curr_1 = az*sin(2*M_PI*(x11[n11]+h2)/periodz+phasez) + offsetz;
			double yh2;
			double a2 = h2 - h2*h2/2;
			yh2 = (1-a2) * y11[n11] + (a2-h2/2) * curr_0 + (h2/2) * curr_1 + a2;
			if(yh2 <= yth) yh2 = yreset;
			yh2 = (1-a2) * yh2 + (a2-h2/2) * curr_1 + (h2/2) * curr_2 + a2;
			if(yh2 <= yth) yh2 = yreset;
			
			double trunc = abs(1./3. * (yh2 - yh));
			if((0.5*epsilon*h/2 > trunc || 2.*epsilon*h/2 < trunc) && trunc > 0. && iter < 100) {
				std::cout << "changing h from " << h << " to " << h*pow(epsilon*h/2/trunc, 1./3.) << " at step " << n11 << std::endl;
				h = h*pow(epsilon*h/2/trunc, 1./3.);
				iter ++;
			}
			else {
				y11.push_back(yh2);
				x11.push_back(x11[n11]+h);
				break;
			}
		}
		n11++;
	}
	
	std::cout << std::endl;
	
	/*epsilon = 1e-10;
	// RK4
	std::vector<double> x2, y2;
	x2.push_back(x0); y2.push_back(y0);
	int n2 = 0;
	h = h0;
	while(x2[n2] < simtime) {
		while(true) {
			double curr_0 = az*sin(2*M_PI*x1[n1]/periodz+phasez) + offsetz;
			double curr_1 = az*sin(2*M_PI*(x1[n1]+h/2)/periodz+phasez) + offsetz;
			double curr_2 = az*sin(2*M_PI*(x1[n1]+h)/periodz+phasez) + offsetz;
			double yh;
			double a1 = -h*h*h*h/24. + h*h*h/6. - h*h/2. + h, a2 = h*h*h/12. - h*h/3. + 2.*h/3.;
			yh = (1-a1) * y1[n1] + (a1-a2-h/6) * curr_0 + (a2) * curr_1 + (h/6) * curr_2 + a1;
			if(yh <= yth) yh = yreset;
			
			double h2 = h/2;
			double curr_12 = az*sin(2*M_PI*(x1[n1]+h2/2)/periodz+phasez) + offsetz;
			double curr_32 = az*sin(2*M_PI*(x1[n1]+h2/2*3)/periodz+phasez) + offsetz;
			double yh2;
			double a12 = -h2*h2*h2*h2/24. + h2*h2*h2/6. - h2*h2/2. + h2, a22 = h2*h2*h2/12. - h2*h2/3. + 2.*h2/3.;
			yh2 = (1-a12) * y1[n1] + (a12-a22-h2/6) * curr_0 + (a22) * curr_12 + (h2/6) * curr_1 + a12;
			if(yh2 <= yth) yh2 = yreset;
			yh2 = (1-a12) * yh2 + (a12-a22-h2/6) * curr_1 + (a22) * curr_32 + (h2/6) * curr_2 + a12;
			if(yh2 <= yth) yh2 = yreset;
			
			double trunc = abs(1./15. * (yh2 - yh));
			if(0.5*epsilon*h/2 > trunc || 2.*epsilon*h/2 < trunc) {
				std::cout << "changing h from " << h << " to " << h*pow(epsilon*h/2/trunc, 1./5.) << " at step " << n2 << std::endl;
				h = h*pow(epsilon*h/2/trunc, 1./5.);
			}
			else {
				y2.push_back(yh2);
				x2.push_back(x2[n2]+h);
				break;
			}
		}
		n2++;
	}*/
	
	std::vector<double> y3, x3, spkt3;
	std::vector<double> z3;
	for(int i = 0; i < N0; i++) z3.push_back(az*sin(2*M_PI*(x0+h0*i)/periodz+phasez) + offsetz);
	y3 = int_lif::Heun(y0, h0, N0, z3, params);
	x3 = int_lif::utils::linspace(x0, x0 + simtime, N0);
	
	myStyle();
	
	TCanvas* canvas = new TCanvas("canvas0", "Canvas0", 600, 400);
	TMultiGraph* multigraph = new TMultiGraph();
	TLegend* legend = new TLegend(0.75, 0.75, 0.85, 0.85);
	
	addToMultigraph(multigraph, legend, x1, y1, n1, 1, "adaHeun #epsilon = 1e-5", "p", "");
	addToMultigraph(multigraph, legend, x11, y11, n11, 2, "adaHeun #epsilon = 1e-6", "p", "");
	//addToMultigraph(multigraph, legend, x2, y2, n2, 2, "adaRK4", "p", "");
	addToMultigraph(multigraph, legend, x3, y3, N0, 3, "Heun", "p", "");
	
	canvas->cd();
	canvas->SetGrid();
	multigraph->Draw("AP");
	multigraph->SetTitle("LIF integrators comparison;t [Cm/g];V [Vrest]");
	legend->Draw();
	canvas->SaveAs("integration_comparison_ada.pdf");
	
	
	
	// Richardson error
	/*// Heun
	std::vector<double> x12, y12, error12;
	x12.push_back(x0); y12.push_back(y0); error12.push_back(0.);
	h = h0 * 2;
	N0 = N0 / 2;
	for(int n = 1; n < N0; n++) {
		double curr_0 = az*sin(2*M_PI*x12[n-1]/periodz+phasez) + offsetz;
		double curr_2 = az*sin(2*M_PI*(x12[n-1]+h)/periodz+phasez) + offsetz;
		double yh;
		double a = h - h*h/2;
		yh = (1-a) * y12[n-1] + (a-h/2) * curr_0 + (h/2) * curr_2 + a;
		if(yh <= yth) yh = yreset;
		
		double h2 = h/2;
		double curr_1 = az*sin(2*M_PI*(x12[n-1]+h2)/periodz+phasez) + offsetz;
		double yh2;
		double a2 = h2 - h2*h2/2;
		yh2 = (1-a2) * y12[n-1] + (a2-h2/2) * curr_0 + (h2/2) * curr_1 + a2;
		if(yh2 <= yth) yh2 = yreset;
		yh2 = (1-a2) * yh2 + (a2-h2/2) * curr_1 + (h2/2) * curr_2 + a2;
		if(yh2 <= yth) yh2 = yreset;
		
		double trunc = abs(1./3. * (yh2 - yh));
		error12.push_back(trunc);
		y12.push_back(yh2);
		x12.push_back(x12[n-1]+h);
	}
	
	TCanvas* canvas2 = new TCanvas("canvas2", "Canvas2", 600, 400);
	TMultiGraph* multigraph2 = new TMultiGraph();
	TLegend* legend2 = new TLegend(0.75, 0.75, 0.85, 0.85);
	
	addToMultigraph(multigraph2, legend2, x12, error12, N0, 3, "Heun", "p", "");
	
	canvas2->cd();
	canvas2->SetGrid();
	multigraph2->Draw("AP");
	multigraph2->SetTitle("LIF integrators Richardson error;t [Cm/g];error [Vrest]");
	canvas2->SetLogy();
	legend2->Draw();
	canvas2->SaveAs("integration_comparison_ada_richardson.pdf");*/
	
	simtime = 1e0;
	int nNvec = 400;
	double Nvec0 = 2., Nvec1 = 7.;
	std::vector<int> Nvec; // occhio che sono double ma sono interi
	for(int i = 0; i < nNvec; i++) {
		int newN = floor(pow(10, Nvec0+((Nvec1-Nvec0)*i)/nNvec));
		if(i == 0 || newN != Nvec[Nvec.size()-1]) {
			Nvec.push_back(newN);
			//std::cout << Nvec[Nvec.size()-1] << " ";
		}
	}
	
	std::vector<double> richvec1;
	int percent1 = 0, ii1 = -1;
	for(int N : Nvec) {
		int_lif::utils::printPercent(ii1++, percent1, Nvec.size(), "");
		h = simtime/(N-1);
		double rich = 0.;
		
		std::vector<double> x12, y12;
		x12.push_back(x0); y12.push_back(y0);
		for(int n = 1; n < N; n++) {
			double curr_0 = az*sin(2*M_PI*x12[n-1]/periodz+phasez) + offsetz;
			double curr_2 = az*sin(2*M_PI*(x12[n-1]+h)/periodz+phasez) + offsetz;
			double yh;
			double a = h - h*h/2;
			yh = (1-a) * y12[n-1] + (a-h/2) * curr_0 + (h/2) * curr_2 + a;
			if(yh <= yth) yh = yreset;
			
			double h2 = h/2;
			double curr_1 = az*sin(2*M_PI*(x12[n-1]+h2)/periodz+phasez) + offsetz;
			double yh2;
			double a2 = h2 - h2*h2/2;
			yh2 = (1-a2) * y12[n-1] + (a2-h2/2) * curr_0 + (h2/2) * curr_1 + a2;
			if(yh2 <= yth) yh2 = yreset;
			yh2 = (1-a2) * yh2 + (a2-h2/2) * curr_1 + (h2/2) * curr_2 + a2;
			if(yh2 <= yth) yh2 = yreset;
			
			double trunc = abs(1./3. * (yh2 - yh));
			y12.push_back(yh2);
			x12.push_back(x12[n-1]+h);
			rich += trunc / nNvec;
		}
		richvec1.push_back(rich);
	}
	
	TCanvas* canvas10 = new TCanvas("canvas10", "Canvas10", 600, 400);
	TMultiGraph* multigraph10 = new TMultiGraph();
	TLegend* legend10 = new TLegend(0.15, 0.75, 0.25, 0.85);
	
	addToMultigraph(multigraph10, legend10, Nvec, richvec1, Nvec.size(), simtime, 3, "Heun", "p", "");
	
	canvas10->cd();
	canvas10->SetGrid();
	multigraph10->Draw("AP");
	multigraph10->SetTitle("Richardson Error;h [];error [V]");
	legend10->Draw();
	canvas10->SaveAs("integration_comparison_ada_richardson.pdf");
	multigraph10->GetXaxis()->SetLimits(0.5e-7, 2e-2);
	canvas10->SetLogx();
	canvas10->SaveAs("integration_comparison_ada_richardson_lx.pdf");
	multigraph10->SetMaximum(pow(10, ceil(log10(multigraph10->GetHistogram()->GetMaximum()))));
	multigraph10->SetMinimum(pow(10, -15));
	canvas10->SetLogy();
	canvas10->SaveAs("integration_comparison_ada_richardson_lxly.pdf");
	
}