#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <filesystem>
namespace fs = std::filesystem;

#include "bound_cond_vecs.h"
#include "derivators.h"
#include "integrators.h"
#include "wave_plots.h"
#include "tridiag.h"

void inline printPercent(int ii, int& percent, const int& N, std::string prefix = "") {
	if(100 * ii / N > percent) {
		percent = 100 * ii / N;
		std::cout << "\r" << prefix << percent << "%";
		std::flush(std::cout);
	}
}

std::vector<bound_cond_vecs::BoundCondVec<double>> HH_FTCS(double t0, double dt, double nsteps, const bound_cond_vecs::BoundCondVec<double> & x, const bound_cond_vecs::BoundCondVec<double> & u0, const bound_cond_vecs::BoundCondVec<double> & m0, const bound_cond_vecs::BoundCondVec<double> & h0, const bound_cond_vecs::BoundCondVec<double> & n0, double pars[8]) {
	double cm = pars[0], D_2_R = pars[1], g_na = pars[2], e_na = pars[3], g_k = pars[4], e_k = pars[5], g_i = pars[6], e_i = pars[7];
	int N = x.len();
	double dx = x[1] - x[0];
	std::vector<bound_cond_vecs::BoundCondVec<double>> u;
	u.reserve(nsteps + 1);
	u.push_back(u0);
	std::vector<double> t;
	t.reserve(nsteps + 1);
	t.push_back(t0);
	std::vector<bound_cond_vecs::BoundCondVec<double>> m;
	m.reserve(nsteps + 1);
	m.push_back(n0);
	std::vector<bound_cond_vecs::BoundCondVec<double>> h;
	h.reserve(nsteps + 1);
	h.push_back(h0);
	std::vector<bound_cond_vecs::BoundCondVec<double>> n;
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

std::vector<bound_cond_vecs::BoundCondVec<double>> HH_FTCS_implicit(double t0, double dt, double nsteps, const bound_cond_vecs::BoundCondVec<double> & x, const bound_cond_vecs::BoundCondVec<double> & u0, const bound_cond_vecs::BoundCondVec<double> & m0, const bound_cond_vecs::BoundCondVec<double> & h0, const bound_cond_vecs::BoundCondVec<double> & n0, double pars[8]) {
	double cm = pars[0], D_2_R = pars[1], g_na = pars[2], e_na = pars[3], g_k = pars[4], e_k = pars[5], g_i = pars[6], e_i = pars[7];
	int N = x.len();
	double dx = x[1] - x[0];
	std::vector<bound_cond_vecs::BoundCondVec<double>> u;
	u.reserve(nsteps + 1);
	u.push_back(u0);
	std::vector<double> t;
	t.reserve(nsteps + 1);
	t.push_back(t0);
	std::vector<bound_cond_vecs::BoundCondVec<double>> m;
	m.reserve(nsteps + 1);
	m.push_back(n0);
	std::vector<bound_cond_vecs::BoundCondVec<double>> h;
	h.reserve(nsteps + 1);
	h.push_back(h0);
	std::vector<bound_cond_vecs::BoundCondVec<double>> n;
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
		bound_cond_vecs::BoundCondVec<double> new_u(N, x.getMode()), new_m(N, x.getMode()), new_h(N, x.getMode()), new_n(N, x.getMode());
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
			
			a[jj] = - eta_1 * eta_2 / 2; if(x.getMode() == ABSORBANT_BC && jj == 0) a[jj] = 0; if(x.getMode() == REFLECTIVE_BC && jj == 0) a[jj] = -a[jj];
			b[jj] = 1 + eta_1 * eta_2;
			c[jj] = - eta_1 * eta_2 / 2; if(x.getMode() == ABSORBANT_BC && jj == N-1) c[jj] = 0; if(x.getMode() == REFLECTIVE_BC && jj == N-1) c[jj] = -c[jj];
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
			
			a[jj] = - eta_1 * eta_2 / 2; if(x.getMode() == ABSORBANT_BC && jj == 0) a[jj] = 0; if(x.getMode() == REFLECTIVE_BC && jj == 0) a[jj] = -a[jj];
			b[jj] = 1 + eta_1 * eta_2;
			c[jj] = - eta_1 * eta_2 / 2; if(x.getMode() == ABSORBANT_BC && jj == N-1) c[jj] = 0; if(x.getMode() == REFLECTIVE_BC && jj == N-1) c[jj] = -c[jj];
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

int main() {
	double t0 = 0., dt = 0.01;
	int nsteps = 10000, nx = 65;
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
	
	bound_cond_vecs::BoundCondVec<double> u0(nx, x.getMode()), m0(nx, x.getMode()), h0(nx, x.getMode()), n0(nx, x.getMode());
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
	}
	for(int ii = 0; ii < nx; ii++) {
		u0[ii] += exp(- (x[ii]-0.5*100) * (x[ii]-0.5*100) / (2. * 0.1*100 * 0.1*100))*-30; // gaussian
		// u0[ii] += sin(2. * M_PI * 1. * x[ii]/100.)*10; // cosine
	}
	
	// auto u1 = HH_FTCS_implicit(t0, dt, nsteps, x, u0, m0, h0, n0, pars);
	auto u1 = HH_CN(t0, dt, nsteps, x, u0, m0, h0, n0, pars, m, h, n);
	std::cout << std::endl;
	
	fs::current_path(fs::current_path() / "measures");
	double minu = -100, maxu = 50;
	
	waveplots::plot(u1, t0, dt, nsteps+1, x0, dx, nx, "HH_test_SURF", SURF_PLOT, minu, maxu);
	waveplots::plot(u1, t0, dt, nsteps+1, x0, dx, nx, "HH_test_CONT", CONT_PLOT, minu, maxu);
	waveplots::plot(u1, t0, dt, nsteps+1, x0, dx, nx, "HH_test_COLZ", COLZ_PLOT, minu, maxu);
	// waveplots::plotFFT(u1, t0, dt, nsteps, "HH_test_FFT_SURF", SURF_PLOT);
}