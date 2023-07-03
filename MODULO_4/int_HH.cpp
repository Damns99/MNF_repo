#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <filesystem>
namespace fs = std::filesystem;
#include <fstream>

#include "bound_cond_vecs.h"
#include "derivators.h"
#include "integrators.h"
#include "wave_plots.h"
#include "tridiag.h"
#include "mystyle.h"
#include "cmdline_parser.h"

void inline printPercent(int ii, int& percent, const int& N, std::string prefix = "") {
	if(100 * ii / N > percent) {
		percent = 100 * ii / N;
		std::cout << "\r" << prefix << percent << "%";
		std::flush(std::cout);
	}
}

std::vector<bound_cond_vecs::BoundCondVec<double>> HH_FTCS(double t0, double dt, double nsteps, const bound_cond_vecs::BoundCondVec<double> & x, const bound_cond_vecs::BoundCondVec<double> & u0, const bound_cond_vecs::BoundCondVec<double> & m0, const bound_cond_vecs::BoundCondVec<double> & h0, const bound_cond_vecs::BoundCondVec<double> & n0, double pars[8], std::vector<bound_cond_vecs::BoundCondVec<double>> & m, std::vector<bound_cond_vecs::BoundCondVec<double>> & h, std::vector<bound_cond_vecs::BoundCondVec<double>> & n) {
	double cm = pars[0], D_2_R = pars[1], g_na = pars[2], e_na = pars[3], g_k = pars[4], e_k = pars[5], g_i = pars[6], e_i = pars[7];
	int N = x.len();
	double dx = x[1] - x[0];
	std::vector<bound_cond_vecs::BoundCondVec<double>> u;
	u.reserve(nsteps + 1);
	u.push_back(u0);
	std::vector<double> t;
	t.reserve(nsteps + 1);
	t.push_back(t0);
	m.reserve(nsteps + 1);
	m.push_back(n0);
	h.reserve(nsteps + 1);
	h.push_back(h0);
	n.reserve(nsteps + 1);
	n.push_back(n0);
	
	int percent = 0;
	for(int ii = 1; ii <= nsteps; ii++) {
		printPercent(ii, percent, nsteps, "");
		bound_cond_vecs::BoundCondVec<double> new_u(N), new_m(N), new_h(N), new_n(N);
		
		for(int jj = 0; jj < N; jj++) {
			double V_j = u[ii-1][jj];
			
			double a_m, a_h, a_n, b_m, b_h, b_n;
			a_m = 0.1 * (V_j + 40) / (1 - exp(-0.1 * (V_j + 40))); // [mV/ms] V è in mV!!!!!!!!!
			b_m = 4. * exp(-0.0556 * (V_j + 65));
			a_h = 0.07 * exp(-0.05 * (V_j + 65));
			b_h = 1. / (1. + exp(-0.1 * (V_j + 35)));
			a_n = 0.01 * (V_j + 55) / (1 - exp(-0.1 * (V_j + 55)));
			b_n = 0.125 * exp(-0.0125 * (V_j + 65));
			
			double t_m, t_h, t_n, l_m, l_h, l_n;
			t_m = 1. / (a_m + b_m);
			l_m = a_m * t_m;
			t_h = 1. / (a_h + b_h);
			l_h = a_h * t_h;
			t_n = 1. / (a_n + b_n);
			l_n = a_n * t_n;
			
			double mm, hh, nn;
			mm = (1 - dt/t_m) * m[ii-1][jj] + dt*l_m/t_m;
			hh = (1 - dt/t_h) * h[ii-1][jj] + dt*l_h/t_h;
			nn = (1 - dt/t_n) * n[ii-1][jj] + dt*l_n/t_n;
			
			new_m[jj] = mm;
			new_h[jj] = hh;
			new_n[jj] = nn;
			
			double eta_1 = dt/cm, eta_2 = D_2_R/dx/dx;
			
			double F_j = g_na*mm*mm*mm*hh*(e_na-V_j) + g_k*nn*nn*nn*nn*(e_k-V_j) + g_i*(e_i-V_j);
			
			new_u[jj] = eta_1 * eta_2 * (F_j / eta_2 + u[ii-1][jj+1] + u[ii-1][jj-1] + (1/eta_1/eta_2 - 2) * V_j);
			
		}
		
		u.push_back(new_u);
		m.push_back(new_m);
		h.push_back(new_h);
		n.push_back(new_n);
		t.push_back(t0 + ii * dt);
	}
	
	return u;
}

std::vector<bound_cond_vecs::BoundCondVec<double>> HH_FTCS_implicit(double t0, double dt, double nsteps, const bound_cond_vecs::BoundCondVec<double> & x, const bound_cond_vecs::BoundCondVec<double> & u0, const bound_cond_vecs::BoundCondVec<double> & m0, const bound_cond_vecs::BoundCondVec<double> & h0, const bound_cond_vecs::BoundCondVec<double> & n0, double pars[8], std::vector<bound_cond_vecs::BoundCondVec<double>> & m, std::vector<bound_cond_vecs::BoundCondVec<double>> & h, std::vector<bound_cond_vecs::BoundCondVec<double>> & n) {
	double cm = pars[0], D_2_R = pars[1], g_na = pars[2], e_na = pars[3], g_k = pars[4], e_k = pars[5], g_i = pars[6], e_i = pars[7];
	int N = x.len();
	double dx = x[1] - x[0];
	std::vector<bound_cond_vecs::BoundCondVec<double>> u;
	u.reserve(nsteps + 1);
	u.push_back(u0);
	std::vector<double> t;
	t.reserve(nsteps + 1);
	t.push_back(t0);
	m.reserve(nsteps + 1);
	m.push_back(n0);
	h.reserve(nsteps + 1);
	h.push_back(h0);
	n.reserve(nsteps + 1);
	n.push_back(n0);
	
	int percent = 0;
	for(int ii = 1; ii <= nsteps; ii++) {
		printPercent(ii, percent, nsteps, "");
		bound_cond_vecs::BoundCondVec<double> new_u(N), new_m(N), new_h(N), new_n(N);
		double a[N], b[N], c[N], d[N], uu[N];
		
		for(int jj = 0; jj < N; jj++) {
			double V_j = u[ii-1][jj];
			
			double a_m, a_h, a_n, b_m, b_h, b_n;
			a_m = 0.1 * (V_j + 40) / (1 - exp(-0.1 * (V_j + 40))); // [mV/ms] V è in mV!!!!!!!!!
			b_m = 4. * exp(-0.0556 * (V_j + 65));
			a_h = 0.07 * exp(-0.05 * (V_j + 65));
			b_h = 1. / (1. + exp(-0.1 * (V_j + 35)));
			a_n = 0.01 * (V_j + 55) / (1 - exp(-0.1 * (V_j + 55)));
			b_n = 0.125 * exp(-0.0125 * (V_j + 65));
			
			double t_m, t_h, t_n, l_m, l_h, l_n;
			t_m = 1. / (a_m + b_m);
			l_m = a_m * t_m;
			t_h = 1. / (a_h + b_h);
			l_h = a_h * t_h;
			t_n = 1. / (a_n + b_n);
			l_n = a_n * t_n;
			
			double mm, hh, nn;
			mm = (1 - dt/t_m) * m[ii-1][jj] + dt*l_m/t_m;
			hh = (1 - dt/t_h) * h[ii-1][jj] + dt*l_h/t_h;
			nn = (1 - dt/t_n) * n[ii-1][jj] + dt*l_n/t_n;
			
			new_m[jj] = mm;
			new_h[jj] = hh;
			new_n[jj] = nn;
			
			double eta_1 = dt/cm, eta_2 = D_2_R/dx/dx;
			
			double F_j = g_na*mm*mm*mm*hh*(e_na-V_j) + g_k*nn*nn*nn*nn*(e_k-V_j) + g_i*(e_i-V_j);
			
			a[jj] = - eta_1 * eta_2;
			b[jj] = 1 + 2 * eta_1 * eta_2;
			c[jj] = - eta_1 * eta_2;
			d[jj] = eta_1 * F_j + V_j;
			
		}
		
		triDiag::cTriDiagSolve(N, a, b, c, d, uu);
		for(int jj = 0; jj < N; jj++) new_u[jj] = uu[jj];
		//if(ii == 1) for(int jj = 0; jj < N; jj++) std::cout << "d[" << jj << "] = " << uu[jj] << std::endl;
		
		u.push_back(new_u);
		m.push_back(new_m);
		h.push_back(new_h);
		n.push_back(new_n);
		t.push_back(t0 + ii * dt);
	}
	
	return u;
}

std::vector<bound_cond_vecs::BoundCondVec<double>> HH_CN(double t0, double dt, double nsteps, const bound_cond_vecs::BoundCondVec<double> & x, const bound_cond_vecs::BoundCondVec<double> & u0, const bound_cond_vecs::BoundCondVec<double> & m0, const bound_cond_vecs::BoundCondVec<double> & h0, const bound_cond_vecs::BoundCondVec<double> & n0, double pars[8], std::vector<bound_cond_vecs::BoundCondVec<double>> & m, std::vector<bound_cond_vecs::BoundCondVec<double>> & h, std::vector<bound_cond_vecs::BoundCondVec<double>> & n) {
	double cm = pars[0], D_2_R = pars[1], g_na = pars[2], e_na = pars[3], g_k = pars[4], e_k = pars[5], g_i = pars[6], e_i = pars[7];
	int N = x.len();
	double dx = x[1] - x[0];
	std::vector<bound_cond_vecs::BoundCondVec<double>> u;
	u.reserve(nsteps + 1);
	u.push_back(u0);
	std::vector<double> t;
	t.reserve(nsteps + 1);
	t.push_back(t0);
	m.reserve(nsteps + 1);
	m.push_back(m0);
	h.reserve(nsteps + 1);
	h.push_back(h0);
	n.reserve(nsteps + 1);
	n.push_back(n0);
	
	int percent = 0;
	for(int ii = 1; ii <= nsteps; ii++) {
		printPercent(ii, percent, nsteps, "");
		bound_cond_vecs::BoundCondVec<double> new_u(N, x.getMode(), u0.getPlaceholder()), new_m(N, x.getMode(), m0.getPlaceholder()), new_h(N, x.getMode(), h0.getPlaceholder()), new_n(N, x.getMode(), n0.getPlaceholder());
		double a[N], b[N], c[N], d[N], uu[N];
		
		for(int jj = 0; jj < N; jj++) {
			double V_j = u[ii-1][jj];
			
			double a_m, a_h, a_n, b_m, b_h, b_n;
			a_m = 0.1 * (V_j + 40) / (1 - exp(-0.1 * (V_j + 40))); // [mV/ms] V è in mV!!!!!!!!!
			b_m = 4. * exp(-0.0556 * (V_j + 65));
			a_h = 0.07 * exp(-0.05 * (V_j + 65));
			b_h = 1. / (1. + exp(-0.1 * (V_j + 35)));
			a_n = 0.01 * (V_j + 55) / (1 - exp(-0.1 * (V_j + 55)));
			b_n = 0.125 * exp(-0.0125 * (V_j + 65));
			
			double t_m, t_h, t_n, l_m, l_h, l_n;
			t_m = 1. / (a_m + b_m);
			l_m = a_m * t_m;
			t_h = 1. / (a_h + b_h);
			l_h = a_h * t_h;
			t_n = 1. / (a_n + b_n);
			l_n = a_n * t_n;
			
			double mm, hh, nn;
			mm = (1 - dt/t_m) * m[ii-1][jj] + dt*l_m/t_m;
			hh = (1 - dt/t_h) * h[ii-1][jj] + dt*l_h/t_h;
			nn = (1 - dt/t_n) * n[ii-1][jj] + dt*l_n/t_n;
			
			new_m[jj] = mm;
			new_h[jj] = hh;
			new_n[jj] = nn;
			
			double eta_1 = dt/cm, eta_2 = D_2_R/dx/dx;
			
			double F_j = g_na*mm*mm*mm*hh*(e_na-V_j) + g_k*nn*nn*nn*nn*(e_k-V_j) + g_i*(e_i-V_j);
			
			a[jj] = - eta_1 * eta_2 / 2;
			if(x.getMode() == ABSORBANT_BC && jj == 0) a[jj] = 0;
			if(x.getMode() == REFLECTIVE_BC && jj == 0) a[jj] = 0; if(x.getMode() == REFLECTIVE_BC && jj == N-1) a[jj] = 2*a[jj];
			b[jj] = 1 + eta_1 * eta_2;
			c[jj] = - eta_1 * eta_2 / 2;
			if(x.getMode() == ABSORBANT_BC && jj == N-1) c[jj] = 0;
			if(x.getMode() == REFLECTIVE_BC && jj == 0) c[jj] = 2*c[jj]; if(x.getMode() == REFLECTIVE_BC && jj == N-1) c[jj] = 0;
			d[jj] = eta_1 * F_j + V_j * (1 - eta_1 * eta_2) + eta_1 * eta_2 / 2 * (u[ii-1][jj-1] + u[ii-1][jj+1]);
			
		}
		
		triDiag::cTriDiagSolve(N, a, b, c, d, uu);
		for(int jj = 0; jj < N; jj++) new_u[jj] = uu[jj];
		
		for(int jj = 0; jj < N; jj++) {
			double V_j = new_u[jj];
			
			double a_m, a_h, a_n, b_m, b_h, b_n;
			a_m = 0.1 * (V_j + 40) / (1 - exp(-0.1 * (V_j + 40))); // [mV/ms] V è in mV!!!!!!!!!
			b_m = 4. * exp(-0.0556 * (V_j + 65));
			a_h = 0.07 * exp(-0.05 * (V_j + 65));
			b_h = 1. / (1. + exp(-0.1 * (V_j + 35)));
			a_n = 0.01 * (V_j + 55) / (1 - exp(-0.1 * (V_j + 55)));
			b_n = 0.125 * exp(-0.0125 * (V_j + 65));
			
			double t_m, t_h, t_n, l_m, l_h, l_n;
			t_m = 1. / (a_m + b_m);
			l_m = a_m * t_m;
			t_h = 1. / (a_h + b_h);
			l_h = a_h * t_h;
			t_n = 1. / (a_n + b_n);
			l_n = a_n * t_n;
			
			double mm, hh, nn;
			mm = m[ii-1][jj]/2 + (1./2.-dt/2/t_m) * new_m[jj] + dt*l_m/2/t_m;
			hh = h[ii-1][jj]/2 + (1./2.-dt/2/t_h) * new_h[jj] + dt*l_h/2/t_h;
			nn = n[ii-1][jj]/2 + (1./2.-dt/2/t_n) * new_n[jj] + dt*l_n/2/t_n;
			
			new_m[jj] = mm;
			new_h[jj] = hh;
			new_n[jj] = nn;
			
			double eta_1 = dt/cm, eta_2 = D_2_R/dx/dx;
			
			double F_j = g_na*mm*mm*mm*hh*(e_na-V_j) + g_k*nn*nn*nn*nn*(e_k-V_j) + g_i*(e_i-V_j);
			
			a[jj] = - eta_1 * eta_2 / 2;
			if(x.getMode() == ABSORBANT_BC && jj == 0) a[jj] = 0;
			if(x.getMode() == REFLECTIVE_BC && jj == 0) a[jj] = 0; if(x.getMode() == REFLECTIVE_BC && jj == N-1) a[jj] = 2*a[jj];
			b[jj] = 1 + eta_1 * eta_2;
			c[jj] = - eta_1 * eta_2 / 2;
			if(x.getMode() == ABSORBANT_BC && jj == N-1) c[jj] = 0;
			if(x.getMode() == REFLECTIVE_BC && jj == 0) c[jj] = 2*c[jj]; if(x.getMode() == REFLECTIVE_BC && jj == N-1) c[jj] = 0;
			d[jj] = eta_1 * F_j + V_j * (1 - eta_1 * eta_2) + eta_1 * eta_2 / 2 * (u[ii-1][jj-1] + u[ii-1][jj+1]);
			
		}
		
		triDiag::cTriDiagSolve(N, a, b, c, d, uu);
		for(int jj = 0; jj < N; jj++) new_u[jj] = uu[jj];
		
		u.push_back(new_u);
		m.push_back(new_m);
		h.push_back(new_h);
		n.push_back(new_n);
		t.push_back(t0 + ii * dt);
	}
	return u;
}

int main(int argc, char* argv[]) {
	double sigma, ampl;
	cmdlineParser::CmdlineParser parser;
    parser.addOptParameter<double>("sigma", &sigma, 0.35*100, "[double] StdDev of the gaussian starting condition.");
    parser.addOptParameter<double>("ampl", &ampl, -5500, "[double] Amplitude of the gaussian starting condition.");
	if (parser.parseAll(argc, argv) == HELP_RETURN) return 0;
    parser.kickOff(argv[0]);
	
	double t0 = 0., dt = 0.01;
	int nsteps = 5000, nx = 301;
	double x0 = -1.*100., dx = 2.*100. / (nx-1);
	bound_cond_vecs::BoundCondVec<double> x = integrators::linspace(x0, x0 + (nx-1) * dx, nx, PERIODIC_BC);
	
	double pars[8] = {
		100e-6, // Cm [mF/cm = ms/ohm/cm]
		0.0238/140/2, // D/2/R [cm^2/ohm]
		12000e-6, // g_Na [S/cm = 1/ohm/cm]
		50, // e_Na [mV]
		3600e-6, // g_K [S/cm = 1/ohm/cm]
		-77, // e_K [mV]
		30e-6, // g_i [S/cm = 1/ohm/cm]
		-54.402 // e_i [mV]		
	};
	
	bound_cond_vecs::BoundCondVec<double> u0(nx, x.getMode(), -65), m0(nx, x.getMode()), h0(nx, x.getMode()), n0(nx, x.getMode());
	std::vector<bound_cond_vecs::BoundCondVec<double>> m, h, n;
	for(int ii = 0; ii < nx; ii++) {
		u0[ii] = -65;
		m0[ii] = 0; // 0.724143431;
		h0[ii] = 0; // 0.065423336;
		n0[ii] = 0; // 0.905660902;
	}
	
	bool startFromEquilibrium = true;
	int nskip = int(50./dt);
	if(startFromEquilibrium) {
		std::vector<bound_cond_vecs::BoundCondVec<double>> m1, h1, n1;
		auto u1 = HH_CN(t0, dt, nskip, x, u0, m0, h0, n0, pars, m1, h1, n1);
		for(int ii = 0; ii < nx; ii++) {
			u0[ii] = u1[nskip][ii];
			m0[ii] = m1[nskip][ii];
			h0[ii] = h1[nskip][ii];
			n0[ii] = n1[nskip][ii];
		}
		u0.setPlaceholder(u0[0]);
		m0.setPlaceholder(m0[0]);
		h0.setPlaceholder(h0[0]);
		n0.setPlaceholder(n0[0]);
	}
	double mu = 0.2*100;
	for(int ii = 0; ii < nx; ii++) {
		// u0[ii] += (x[ii]-mu) * (x[ii]-mu) < sigma * sigma ? -30 : 0; // constant
		u0[ii] += exp(- (x[ii]-mu) * (x[ii]-mu) / (2. * sigma * sigma)) / (sigma * 2 * M_PI) * ampl; // gaussian
		// u0[ii] += sin(2. * M_PI * 1. * x[ii]/50.)*-8; // cosine
	}
	
	// auto u1 = HH_FTCS_implicit(t0, dt, nsteps, x, u0, m0, h0, n0, pars, m, h, n);
	auto u1 = HH_CN(t0, dt, nsteps, x, u0, m0, h0, n0, pars, m, h, n);
	std::cout << std::endl;
	
	fs::current_path(fs::current_path() / "measures");
	double minu = -100, maxu = 50;
	
	// waveplots::plot(u1, t0, dt, nsteps+1, x0, dx, nx, "HH_test_SURF", SURF_PLOT, minu, maxu);
	waveplots::plot(u1, t0, dt, nsteps+1, x0, dx, nx, "HH_test_CONT", CONT_PLOT, minu, maxu);
	waveplots::plot(u1, t0, dt, nsteps+1, x0, dx, nx, "HH_test_COLZ", COLZ_PLOT, minu, maxu);
	// waveplots::plotFFT(u1, t0, dt, nsteps, "HH_test_FFT_SURF", SURF_PLOT);
	
	double u1_slice[nsteps], m_slice[nsteps], h_slice[nsteps], n_slice[nsteps], time[nsteps];
	int j_slice = nx*1/2;
	for(int ii = 0; ii < nsteps; ii++) {
		u1_slice[ii] = u1[ii][j_slice];
		m_slice[ii] = m[ii][j_slice];
		h_slice[ii] = h[ii][j_slice];
		n_slice[ii] = n[ii][j_slice];
		time[ii] = t0 + ii*dt;
	}
	int nshift = 1;
	double u1_diff[nsteps-nshift], t_u1_max, u1_max = u1[0][j_slice];
	for(int ii = 0; ii < nsteps-nshift; ii++) {
		u1_diff[ii] = u1[ii+nshift][j_slice] - u1[ii][j_slice];
		if(ii > 0) if(u1_diff[ii-1] > 0 && u1_diff[ii] <= 0 && u1_max < u1[ii][j_slice]) {
			t_u1_max = t0 + ii*dt;
			u1_max = u1[ii][j_slice];
		}
	}
	std::cout << "(t_max [ms], V_max [mV]) = (" << t_u1_max << ", " << u1_max << ")" << std::endl;
	
	myStyle();
	
	TCanvas* canvas1 = new TCanvas("canvas1", "Canvas1", 600, 400);
	TGraph* graphu1 = new TGraph(nsteps, time, u1_slice);
	graphu1->SetMarkerColor(1);
	graphu1->SetLineColor(1);
	graphu1->SetTitle(("x = "+std::to_string(x[j_slice])+" cm;time [ms];V [mV]").c_str());
	canvas1->cd();
	canvas1->SetGrid();
	graphu1->Draw("AP");
	canvas1->SaveAs("HH_test_u1_slice.pdf");
	
	TCanvas* canvas2 = new TCanvas("canvas2", "Canvas2", 600, 400);
	TMultiGraph* multigraph = new TMultiGraph();
	TLegend* legend = new TLegend(0.75, 0.75, 0.85, 0.85);
	TGraph* graphm = new TGraph(nsteps, time, m_slice);
	graphm->SetMarkerColor(2);
	graphm->SetLineColor(2);
	multigraph->Add(graphm);
	legend->AddEntry(graphm, "m", "p");
	TGraph* graphh = new TGraph(nsteps, time, h_slice);
	graphh->SetMarkerColor(3);
	graphh->SetLineColor(3);
	multigraph->Add(graphh);
	legend->AddEntry(graphh, "h", "p");
	TGraph* graphn = new TGraph(nsteps, time, n_slice);
	graphn->SetMarkerColor(4);
	graphn->SetLineColor(4);
	multigraph->Add(graphn);
	legend->AddEntry(graphn, "n", "p");
	multigraph->SetTitle(("x = "+std::to_string(x[j_slice])+" cm;time [ms];[mV]").c_str());
	canvas2->cd();
	canvas2->SetGrid();
	multigraph->Draw("AP");
	legend->Draw();
	canvas2->SaveAs("HH_test_mhn_slice.pdf");
	
	TCanvas* canvas3 = new TCanvas("canvas3", "Canvas3", 600, 400);
	TMultiGraph* multigraph2 = new TMultiGraph();
	TLegend* legend2 = new TLegend(0.75, 0.75, 0.85, 0.85);
	TGraph* graphm2 = new TGraph(nsteps, u1_slice, m_slice);
	graphm2->SetMarkerColor(2);
	graphm2->SetLineColor(2);
	multigraph2->Add(graphm2);
	legend2->AddEntry(graphm2, "m", "p");
	TGraph* graphh2 = new TGraph(nsteps, u1_slice, h_slice);
	graphh2->SetMarkerColor(3);
	graphh2->SetLineColor(3);
	multigraph2->Add(graphh2);
	legend2->AddEntry(graphh2, "h", "p");
	TGraph* graphn2 = new TGraph(nsteps, u1_slice, n_slice);
	graphn2->SetMarkerColor(4);
	graphn2->SetLineColor(4);
	multigraph2->Add(graphn2);
	legend2->AddEntry(graphn2, "n", "p");
	multigraph2->SetTitle(("x = "+std::to_string(x[j_slice])+" cm;u1 [mV];[mV]").c_str());
	canvas3->cd();
	canvas3->SetGrid();
	multigraph2->Draw("AP");
	legend2->Draw();
	canvas3->SaveAs("HH_test_u1_vs_mhn_slice.pdf");
	
	std::ofstream outfile;
	outfile.open("HH_test_u1.txt", std::fstream::out);
	outfile << "# x0= " << x0 << " dx= " << dx << " nx= " << nx << std::endl;
	outfile << "# t0= " << t0 << " dt= " << dt << " nt= " << nsteps << std::endl;
	outfile << "# u0 gaussian mean= " << mu << " stdev= " << sigma << " ampl= " << ampl << std::endl;
	for(int ii = 0; ii < nsteps; ii++) {
		for(int jj = 0; jj < nx; jj++) outfile << u1[ii][jj] << " ";
		outfile << std::endl;
	}
	outfile.close();
	outfile.open("HH_test_m.txt", std::fstream::out);
	outfile << "# x0= " << x0 << " dx= " << dx << " nx= " << nx << std::endl;
	outfile << "# t0= " << t0 << " dt= " << dt << " nt= " << nsteps << std::endl;
	outfile << "# u0 gaussian mean= " << mu << " stdev= " << sigma << " ampl= " << ampl << std::endl;
	for(int ii = 0; ii < nsteps; ii++) {
		for(int jj = 0; jj < nx; jj++) outfile << m[ii][jj] << " ";
		outfile << std::endl;
	}
	outfile.close();
	outfile.open("HH_test_h.txt", std::fstream::out);
	outfile << "# x0= " << x0 << " dx= " << dx << " nx= " << nx << std::endl;
	outfile << "# t0= " << t0 << " dt= " << dt << " nt= " << nsteps << std::endl;
	outfile << "# u0 gaussian mean= " << mu << " stdev= " << sigma << " ampl= " << ampl << std::endl;
	for(int ii = 0; ii < nsteps; ii++) {
		for(int jj = 0; jj < nx; jj++) outfile << h[ii][jj] << " ";
		outfile << std::endl;
	}
	outfile.close();
	outfile.open("HH_test_n.txt", std::fstream::out);
	outfile << "# x0= " << x0 << " dx= " << dx << " nx= " << nx << std::endl;
	outfile << "# t0= " << t0 << " dt= " << dt << " nt= " << nsteps << std::endl;
	outfile << "# u0 gaussian mean= " << mu << " stdev= " << sigma << " ampl= " << ampl << std::endl;
	for(int ii = 0; ii < nsteps; ii++) {
		for(int jj = 0; jj < nx; jj++) outfile << n[ii][jj] << " ";
		outfile << std::endl;
	}
	outfile.close();
	
	outfile.open("HH_test_max.txt", std::fstream::out);
	outfile << "# x0= " << x0 << " dx= " << dx << " nx= " << nx << std::endl;
	outfile << "# t0= " << t0 << " dt= " << dt << " nt= " << nsteps << std::endl;
	outfile << "# u0 gaussian mean= " << mu << " stdev= " << sigma << " ampl= " << ampl << std::endl;
	outfile << "# xmax [cm]\ttmax[ms] \tVmax [mV]" << std::endl;
	for(int jj = 0; jj < nx; jj++) {
		double u1_diff, t_u1_max, u1_max = u1[0][jj];
		for(int ii = 0; ii < nsteps-1; ii++) {
			double u1_diff_tmp = u1[ii+1][jj] - u1[ii][jj];
			if(ii > 0) if(u1_diff > 0 && u1_diff_tmp <= 0 && u1_max < u1[ii][jj]) {
				t_u1_max = t0 + ii*dt;
				u1_max = u1[ii][jj];
			}
			u1_diff = u1_diff_tmp;
		}
		outfile << x[jj] << "\t" << t_u1_max << "\t" << u1_max << std::endl;
	}
	outfile.close();
}