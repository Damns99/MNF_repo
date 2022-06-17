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

// equation: du/dt = v * df(u)/dx + w * d2g(u)/dx2 + q

constexpr double M_R = 0.21;
constexpr double D_R = - 1.;

std::vector<bound_cond_vecs::BoundCondVec<double>> v(std::vector<bound_cond_vecs::BoundCondVec<double>> u, double t, bound_cond_vecs::BoundCondVec<double> x) {	
	std::vector<bound_cond_vecs::BoundCondVec<double>> vu;
	bound_cond_vecs::BoundCondVec<double> vu0_tmp(x.len(), x.getMode()), vu1_tmp(x.len(), x.getMode());
	
	for(int jj = 0; jj < x.len(); jj++) {
		vu0_tmp[jj] = 1.;
		vu1_tmp[jj] = 1.;
	}
	
	vu.push_back(vu0_tmp);
	vu.push_back(vu1_tmp);
	return vu;
}

std::vector<bound_cond_vecs::BoundCondVec<double>> w(std::vector<bound_cond_vecs::BoundCondVec<double>> u, double t, bound_cond_vecs::BoundCondVec<double> x) {
	std::vector<bound_cond_vecs::BoundCondVec<double>> wu;
	bound_cond_vecs::BoundCondVec<double> wu0_tmp(x.len(), x.getMode()), wu1_tmp(x.len(), x.getMode());
	
	for(int jj = 0; jj < x.len(); jj++) {
		wu0_tmp[jj] = 0.;
		wu1_tmp[jj] = 0.;
	}
	
	wu.push_back(wu0_tmp);
	wu.push_back(wu1_tmp);
	return wu;
}

std::vector<bound_cond_vecs::BoundCondVec<double>> q(std::vector<bound_cond_vecs::BoundCondVec<double>> u, double t, bound_cond_vecs::BoundCondVec<double> x) {
	std::vector<bound_cond_vecs::BoundCondVec<double>> qu;
	bound_cond_vecs::BoundCondVec<double> qu0_tmp(x.len(), x.getMode()), qu1_tmp(x.len(), x.getMode());
	
	for(int jj = 0; jj < x.len(); jj++) {
		qu0_tmp[jj] = u[1][jj];
		qu1_tmp[jj] = u[0][jj] * (2 * M_R / (x[jj] - D_R) - 1) * 2 * M_R / ((x[jj] - D_R) * (x[jj] - D_R) * (x[jj] - D_R));
	}
	
	qu.push_back(qu0_tmp);
	qu.push_back(qu1_tmp);
	return qu;
}

std::vector<bound_cond_vecs::BoundCondVec<double>> f(std::vector<bound_cond_vecs::BoundCondVec<double>> u, double t, bound_cond_vecs::BoundCondVec<double> x) {
	std::vector<bound_cond_vecs::BoundCondVec<double>> fu;
	bound_cond_vecs::BoundCondVec<double> fu0_tmp(x.len(), x.getMode()), fu1_tmp(x.len(), x.getMode());
	
	for(int jj = 0; jj < x.len(); jj++) {
		fu0_tmp[jj] = - u[0][jj];
		fu1_tmp[jj] = u[1][jj];
	}
	
	fu.push_back(fu0_tmp);
	fu.push_back(fu1_tmp);
	return fu;
}

std::vector<bound_cond_vecs::BoundCondVec<double>> g(std::vector<bound_cond_vecs::BoundCondVec<double>> u, double t, bound_cond_vecs::BoundCondVec<double> x) {
	std::vector<bound_cond_vecs::BoundCondVec<double>> gu;
	bound_cond_vecs::BoundCondVec<double> gu0_tmp(x.len(), x.getMode()), gu1_tmp(x.len(), x.getMode());
	
	for(int jj = 0; jj < x.len(); jj++) {
		gu0_tmp[jj] = 0.;
		gu1_tmp[jj] = 0.;
	}
	
	gu.push_back(gu0_tmp);
	gu.push_back(gu1_tmp);
	return gu;
}

int main() {
	double t0 = 0., dt = 0.002;
	int nsteps = 500, nx = 500;
	double x0 = -0.6, dx = 1. / nx;
	bound_cond_vecs::BoundCondVec<double> x = integrators::linspace(x0, x0 + nx * dx, nx, PERIODIC_BC);
	std::vector<bound_cond_vecs::BoundCondVec<double>> u0;
	bound_cond_vecs::BoundCondVec<double> u0_tmp1(nx, x.getMode()), u0_tmp2(nx, x.getMode());
	double sigma = 0.05, vel = 0.;
	for(int ii = 0; ii < nx; ii++) {
		u0_tmp1[ii] = exp(- x[ii] * x[ii] / (2. * sigma * sigma)); // gaussian
		u0_tmp2[ii] = - (1. + vel) * x[ii] / (sigma * sigma) * exp(- x[ii] * x[ii] / (2. * sigma * sigma)); // gaussian
	}
	u0.push_back(u0_tmp1);
	u0.push_back(u0_tmp2);
	
	auto u2 = integrators::Lax(t0, dt, nsteps, x, u0, f, derivators::symm_derive, g, derivators::symm_derive_2, v, w, q);	
	auto u4 = integrators::LaxWendroff(t0, dt, nsteps, x, u0, f, derivators::fwd_derive, g, derivators::fwd_derive_2, v, w, q);
	
	fs::current_path(fs::current_path() / "measures");
	double minu[2] = {-0.1, -25.}, maxu[2] = {1., 25.};
	waveplots::plot(u2, t0, dt, nsteps, x0, dx, nx, "grav_Lax_SURF", SURF_PLOT, minu, maxu);
	waveplots::plot(u2, t0, dt, nsteps, x0, dx, nx, "grav_Lax_CONT", CONT_PLOT, minu, maxu);
	waveplots::plot(u2, t0, dt, nsteps, x0, dx, nx, "grav_Lax_COLZ", COLZ_PLOT, minu, maxu);
	waveplots::plotFFT(u2, t0, dt, nsteps, "grav_Lax_FFT_SURF", SURF_PLOT);
	waveplots::plot(u4, t0, dt, nsteps, x0, dx, nx, "grav_LaxWendroff_SURF", SURF_PLOT, minu, maxu);
	waveplots::plot(u4, t0, dt, nsteps, x0, dx, nx, "grav_LaxWendroff_CONT", CONT_PLOT, minu, maxu);
	waveplots::plot(u4, t0, dt, nsteps, x0, dx, nx, "grav_LaxWendroff_COLZ", COLZ_PLOT, minu, maxu);
	waveplots::plotFFT(u4, t0, dt, nsteps, "grav_LaxWendroff_test_FFT_SURF", SURF_PLOT);
}