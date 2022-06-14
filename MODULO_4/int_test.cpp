#include <iostream>
#include <vector>
#include <string>
#include <math.h>

#include "bound_cond_vecs.h"
#include "derivators.h"
#include "integrators.h"
#include "wave_plots.h"

std::vector<bound_cond_vecs::BoundCondVec<double>> f(std::vector<bound_cond_vecs::BoundCondVec<double>> u, double t, bound_cond_vecs::BoundCondVec<double> x) {
	double v = 1.;
	
	std::vector<bound_cond_vecs::BoundCondVec<double>> fu;
	for(int ii = 0; ii < u.size(); ii++) {
		bound_cond_vecs::BoundCondVec<double> fu_tmp(x.len(), x.getMode());
		
		for(int jj = 0; jj < x.len(); jj++) fu_tmp[jj] = - v * u[ii][jj];
		
		fu.push_back(fu_tmp);
	}
	return fu;
}

std::vector<bound_cond_vecs::BoundCondVec<double>> g(std::vector<bound_cond_vecs::BoundCondVec<double>> u, double t, bound_cond_vecs::BoundCondVec<double> x) {
	double nu = 0.01;
	
	std::vector<bound_cond_vecs::BoundCondVec<double>> gu;
	for(int ii = 0; ii < u.size(); ii++) {
		bound_cond_vecs::BoundCondVec<double> gu_tmp(x.len(), x.getMode());
		
		for(int jj = 0; jj < x.len(); jj++) gu_tmp[jj] = nu * u[ii][jj];
		
		gu.push_back(gu_tmp);
	}
	return gu;
}

int main() {
	double t0 = 0., dt = 0.01;
	int nsteps = 50, nx = 100;
	double x0 = 0., dx = 1. / nx;
	bound_cond_vecs::BoundCondVec<double> x = integrators::linspace(x0, x0 + nx * dx, nx);
	std::vector<bound_cond_vecs::BoundCondVec<double>> u0;
	bound_cond_vecs::BoundCondVec<double> u0_tmp(nx, x.getMode());
	for(int ii = 0; ii < nx; ii++) {
		u0_tmp[ii] = cos(2. * M_PI * x[ii]);
	}
	u0.push_back(u0_tmp);
	
	auto u = integrators::FTCS(t0, dt, nsteps, x, u0, f, derivators::symm_derive, g, derivators::symm_derive_2);
	
	double maxu = 2.;
	waveplots::plot(u, t0, dt, nsteps, x0, dx, nx, "first_FTCS_test_SURF", SURF_PLOT, maxu);
	waveplots::plot(u, t0, dt, nsteps, x0, dx, nx, "first_FTCS_test_CONT", CONT_PLOT, maxu);
	waveplots::plot(u, t0, dt, nsteps, x0, dx, nx, "first_FTCS_test_COLZ", COLZ_PLOT, maxu);
}