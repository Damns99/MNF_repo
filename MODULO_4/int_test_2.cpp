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

// simple equation dy/dt = v dy/dx + w d2y/dx2

std::vector<bound_cond_vecs::BoundCondVec<double>> HH_FTCS(double t0, double dt, double nsteps, const bound_cond_vecs::BoundCondVec<double> & x, const bound_cond_vecs::BoundCondVec<double> & u0, double pars[1]) {
	double w = pars[0];
	int N = x.len();
	double dx = x[1] - x[0];
	std::vector<bound_cond_vecs::BoundCondVec<double>> u;
	u.reserve(nsteps + 1);
	u.push_back(u0);
	std::vector<double> t;
	t.reserve(nsteps + 1);
	t.push_back(t0);
	
	int percent = 0;
	for(int ii = 1; ii <= nsteps; ii++) {
		printPercent(ii, percent, nsteps, "");
		bound_cond_vecs::BoundCondVec<double> new_u(N);
		
		for(int jj = 0; jj < N; jj++) {
			double V_j = u[ii-1][jj];
			
			double eta_1 = dt, eta_2 = w/dx/dx;
			
			new_u[jj] = eta_1 * eta_2 * (u[ii-1][jj+1] + u[ii-1][jj-1] + (1/eta_1/eta_2 - 2) * V_j);
			
		}
		
		u.push_back(new_u);
		t.push_back(t0 + ii * dt);
	}
	
	return u;
}

std::vector<bound_cond_vecs::BoundCondVec<double>> HH_FTCS_implicit(double t0, double dt, double nsteps, const bound_cond_vecs::BoundCondVec<double> & x, const bound_cond_vecs::BoundCondVec<double> & u0, double pars[1]) {
	double w = pars[0];
	int N = x.len();
	double dx = x[1] - x[0];
	std::vector<bound_cond_vecs::BoundCondVec<double>> u;
	u.reserve(nsteps + 1);
	u.push_back(u0);
	std::vector<double> t;
	t.reserve(nsteps + 1);
	t.push_back(t0);
	
	int percent = 0;
	for(int ii = 1; ii <= nsteps; ii++) {
		printPercent(ii, percent, nsteps, "");
		bound_cond_vecs::BoundCondVec<double> new_u(N);
		double a[N], b[N], c[N], d[N], uu[N];
		
		for(int jj = 0; jj < N; jj++) {
			double V_j = u[ii-1][jj];
			
			double eta_1 = dt, eta_2 = w/dx/dx;
			
			a[jj] = - eta_1 * eta_2;
			if(x.getMode() == ABSORBANT_BC && jj == 0) a[jj] = 0;
			b[jj] = 1 + 2 * eta_1 * eta_2;
			c[jj] = - eta_1 * eta_2;
			if(x.getMode() == ABSORBANT_BC && jj == N-1) c[jj] = 0;
			d[jj] = V_j;
			
		}
		
		triDiag::cTriDiagSolve(N, a, b, c, d, uu);
		for(int jj = 0; jj < N; jj++) new_u[jj] = uu[jj];
		//if(ii == 1) for(int jj = 0; jj < N; jj++) std::cout << "d[" << jj << "] = " << uu[jj] << std::endl;
		
		u.push_back(new_u);
		t.push_back(t0 + ii * dt);
	}
	
	return u;
}

std::vector<bound_cond_vecs::BoundCondVec<double>> HH_CN(double t0, double dt, double nsteps, const bound_cond_vecs::BoundCondVec<double> & x, const bound_cond_vecs::BoundCondVec<double> & u0, double pars[1]) {
	double w = pars[0];
	int N = x.len();
	double dx = x[1] - x[0];
	std::vector<bound_cond_vecs::BoundCondVec<double>> u;
	u.reserve(nsteps + 1);
	u.push_back(u0);
	std::vector<double> t;
	t.reserve(nsteps + 1);
	t.push_back(t0);
	
	int percent = 0;
	for(int ii = 1; ii <= nsteps; ii++) {
		printPercent(ii, percent, nsteps, "");
		bound_cond_vecs::BoundCondVec<double> new_u(N, x.getMode(), u0.getPlaceholder());
		double a[N], b[N], c[N], d[N], uu[N];
		
		for(int jj = 0; jj < N; jj++) {
			double V_j = u[ii-1][jj];
			
			double eta_1 = dt, eta_2 = w/dx/dx;
			
			a[jj] = - eta_1 * eta_2 / 2;
			if(x.getMode() == ABSORBANT_BC && jj == 0) a[jj] = 0;
			if(x.getMode() == REFLECTIVE_BC && jj == 0) a[jj] = 0; if(x.getMode() == REFLECTIVE_BC && jj == N-1) a[jj] = 2*a[jj];
			b[jj] = 1 + eta_1 * eta_2;
			c[jj] = - eta_1 * eta_2 / 2;
			if(x.getMode() == ABSORBANT_BC && jj == N-1) c[jj] = 0;
			if(x.getMode() == REFLECTIVE_BC && jj == 0) c[jj] = 2*c[jj]; if(x.getMode() == REFLECTIVE_BC && jj == N-1) c[jj] = 0;
			d[jj] = V_j * (1 - eta_1 * eta_2) + eta_1 * eta_2 / 2 * (u[ii-1][jj-1] + u[ii-1][jj+1]);
			
		}
		
		triDiag::cTriDiagSolve(N, a, b, c, d, uu);
		for(int jj = 0; jj < N; jj++) new_u[jj] = uu[jj];
		
		u.push_back(new_u);
		t.push_back(t0 + ii * dt);
	}
	return u;
}

int main() {
	double t0 = 0., dt = 1;
	int nsteps = 500, nx = 301;
	double x0 = -1.*100., dx = 2.*100. / (nx-1);
	bound_cond_vecs::BoundCondVec<double> x = integrators::linspace(x0, x0 + (nx-1) * dx, nx, PERIODIC_BC);
	
	double pars[8] = {
		3.85 // D []	
	};
	
	bound_cond_vecs::BoundCondVec<double> u0(nx, x.getMode(), 0.);
	for(int ii = 0; ii < nx; ii++) {
		u0[ii] = 0.;
	}
	
	for(int ii = 0; ii < nx; ii++) {
		u0[ii] += exp(- (x[ii]-0.5*100) * (x[ii]-0.5*100) / (2. * 0.1*100 * 0.1*100))*30; // gaussian
		// u0[ii] += sin(2. * M_PI * 1. * x[ii]/100.)*10; // cosine
	}
	
	//auto u1 = HH_FTCS(t0, dt, nsteps, x, u0, pars);
	auto u1 = HH_FTCS_implicit(t0, dt, nsteps, x, u0, pars);
	// auto u1 = HH_CN(t0, dt, nsteps, x, u0, pars);
	std::cout << std::endl;
	
	fs::current_path(fs::current_path() / "measures");
	double minu = 0, maxu = 30;
	
	// waveplots::plot(u1, t0, dt, nsteps+1, x0, dx, nx, "HH_test_SURF", SURF_PLOT, minu, maxu);
	waveplots::plot(u1, t0, dt, nsteps+1, x0, dx, nx, "int_test_2_CONT", CONT_PLOT, minu, maxu);
	waveplots::plot(u1, t0, dt, nsteps+1, x0, dx, nx, "int_test_2_COLZ", COLZ_PLOT, minu, maxu);
	// waveplots::plotFFT(u1, t0, dt, nsteps, "HH_test_FFT_SURF", SURF_PLOT);
}