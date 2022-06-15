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

// equation: du/dt = df(u)/dx + d2g(u)/dx2

std::vector<bound_cond_vecs::BoundCondVec<double>> f(std::vector<bound_cond_vecs::BoundCondVec<double>> u, double t, bound_cond_vecs::BoundCondVec<double> x) {
	double v = 0.99;
	
	std::vector<bound_cond_vecs::BoundCondVec<double>> fu;
	for(int ii = 0; ii < u.size(); ii++) {
		bound_cond_vecs::BoundCondVec<double> fu_tmp(x.len(), x.getMode());
		
		for(int jj = 0; jj < x.len(); jj++) fu_tmp[jj] = - v * u[ii][jj];
		
		fu.push_back(fu_tmp);
	}
	return fu;
}

std::vector<bound_cond_vecs::BoundCondVec<double>> g(std::vector<bound_cond_vecs::BoundCondVec<double>> u, double t, bound_cond_vecs::BoundCondVec<double> x) {
	double nu = 0.001;
	
	std::vector<bound_cond_vecs::BoundCondVec<double>> gu;
	for(int ii = 0; ii < u.size(); ii++) {
		bound_cond_vecs::BoundCondVec<double> gu_tmp(x.len(), x.getMode());
		
		for(int jj = 0; jj < x.len(); jj++) gu_tmp[jj] = nu * u[ii][jj];
		
		gu.push_back(gu_tmp);
	}
	return gu;
}

int main() {
	double t0 = 0., dt = 0.02;
	int nsteps = 100, nx = 50;
	double x0 = 0., dx = 1. / nx;
	bound_cond_vecs::BoundCondVec<double> x = integrators::linspace(x0, x0 + nx * dx, nx);
	std::vector<bound_cond_vecs::BoundCondVec<double>> u0;
	bound_cond_vecs::BoundCondVec<double> u0_tmp(nx, x.getMode());
	for(int ii = 0; ii < nx; ii++) {
		u0_tmp[ii] = cos(2. * M_PI * 4. * x[ii]) + 0.5 * sin(2. * M_PI * 2. * x[ii]);
	}
	u0.push_back(u0_tmp);
	
	auto u1 = integrators::FTCS(t0, dt, nsteps, x, u0, f, derivators::symm_derive, g, derivators::symm_derive_2);
	auto u2 = integrators::Lax(t0, dt, nsteps, x, u0, f, derivators::symm_derive, g, derivators::symm_derive_2);	
	auto u3 = integrators::LeapFrog(t0, dt, nsteps, x, u0, f, derivators::symm_derive, g, derivators::symm_derive_2);
	auto u4 = integrators::LaxWendroff(t0, dt, nsteps, x, u0, f, derivators::fwd_derive, g, derivators::symm_derive_2);
	
	fs::current_path(fs::current_path() / "measures");
	double maxu = 2.;
	waveplots::plot(u1, t0, dt, nsteps, x0, dx, nx, "first_FTCS_test_SURF", SURF_PLOT, maxu);
	waveplots::plot(u1, t0, dt, nsteps, x0, dx, nx, "first_FTCS_test_CONT", CONT_PLOT, maxu);
	waveplots::plot(u1, t0, dt, nsteps, x0, dx, nx, "first_FTCS_test_COLZ", COLZ_PLOT, maxu);
	waveplots::plot(u2, t0, dt, nsteps, x0, dx, nx, "first_Lax_test_SURF", SURF_PLOT, maxu);
	waveplots::plot(u2, t0, dt, nsteps, x0, dx, nx, "first_Lax_test_CONT", CONT_PLOT, maxu);
	waveplots::plot(u2, t0, dt, nsteps, x0, dx, nx, "first_Lax_test_COLZ", COLZ_PLOT, maxu);
	waveplots::plot(u3, t0, dt, nsteps, x0, dx, nx, "first_LeapFrog_test_SURF", SURF_PLOT, maxu);
	waveplots::plot(u3, t0, dt, nsteps, x0, dx, nx, "first_LeapFrog_CONT", CONT_PLOT, maxu);
	waveplots::plot(u3, t0, dt, nsteps, x0, dx, nx, "first_LeapFrog_COLZ", COLZ_PLOT, maxu);
	waveplots::plot(u4, t0, dt, nsteps, x0, dx, nx, "first_LaxWendroff_test_SURF", SURF_PLOT, maxu);
	waveplots::plot(u4, t0, dt, nsteps, x0, dx, nx, "first_LaxWendroff_test_CONT", CONT_PLOT, maxu);
	waveplots::plot(u4, t0, dt, nsteps, x0, dx, nx, "first_LaxWendroff_test_COLZ", COLZ_PLOT, maxu);
}