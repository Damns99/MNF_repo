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


std::vector<std::vector<bound_cond_vecs::BoundCondVec<double>>> CrankNicolson(double t0, double dt, double nsteps, const bound_cond_vecs::BoundCondVec<double> & x, const std::vector<bound_cond_vecs::BoundCondVec<double>> & u0, double v, double w, double q) {
	int N = x.len(), ndim = u0.size();
	double dx = x[1] - x[0];
	std::vector<std::vector<bound_cond_vecs::BoundCondVec<double>>> u;
	u.reserve(nsteps + 1);
	u.push_back(u0);
	std::vector<double> t;
	t.reserve(nsteps + 1);
	t.push_back(t0);
	
	for(int ii = 1; ii <= nsteps; ii++) {
		std::vector<bound_cond_vecs::BoundCondVec<double>> new_u;
		new_u.reserve(ndim);
		
		for(int jj = 0; jj < ndim; jj++) {
			double a[N], b[N], c[N], d[N], uu[N];
			
			for(int kk = 0; kk < N; kk++) {
				double alpha = dt / dx * v, beta = dt / dx / dx * w;
				a[kk] = - alpha / 4 - beta / 2;
				b[kk] = 1 + beta;
				c[kk] = + alpha / 4 - beta / 2;
				d[kk] = (alpha / 4 + beta / 2) * u[ii - 1][jj][kk - 1] + (1 - beta) * u[ii - 1][jj][kk] + (- alpha / 4 + beta / 2) * u[ii - 1][jj][kk + 1] + dt * q;
			}
			
			triDiag::triDiagSolve(N, a, b, c, d, uu);
			
			bound_cond_vecs::BoundCondVec<double> tmp_u(N, uu);
			new_u.push_back(tmp_u);
		}
		
		u.push_back(new_u);
		t.push_back(t0 + ii * dt);
	}
	
	return u;
}


int main() {
	double t0 = 0., dt = 0.0001;
	int nsteps = 500, nx = 200;
	double x0 = -0.5, dx = 1. / nx;
	bound_cond_vecs::BoundCondVec<double> x = integrators::linspace(x0, x0 + nx * dx, nx, PERIODIC_BC);
	std::vector<bound_cond_vecs::BoundCondVec<double>> u0;
	bound_cond_vecs::BoundCondVec<double> u0_tmp(nx, x.getMode());
	for(int ii = 0; ii < nx; ii++) {
		
		u0_tmp[ii] = sin(2. * M_PI * 1. * x[ii]); // cosine
		// u0_tmp[ii] = exp(- (x[ii]) * (x[ii]) / (2. * 0.01 * 0.01)); // gaussian
		
	}
	u0.push_back(u0_tmp);
	
	double v = 4., w = 0.5, q = 0.;
	
	auto u1 = CrankNicolson(t0, dt, nsteps, x, u0, v, w, q);
	
	fs::current_path(fs::current_path() / "measures");
	double minu[1] = {-2.}, maxu[2] = {2.};
	
	waveplots::plot(u1, t0, dt, nsteps, x0, dx, nx, "CN_test_SURF", SURF_PLOT, minu, maxu);
	waveplots::plot(u1, t0, dt, nsteps, x0, dx, nx, "CN_test_CONT", CONT_PLOT, minu, maxu);
	waveplots::plot(u1, t0, dt, nsteps, x0, dx, nx, "CN_test_COLZ", COLZ_PLOT, minu, maxu);
	waveplots::plotFFT(u1, t0, dt, nsteps, "CN_test_FFT_SURF", SURF_PLOT);
}