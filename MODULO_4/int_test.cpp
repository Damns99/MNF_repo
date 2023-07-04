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
#include "mystyle.h"

constexpr double vv = 0.5;
constexpr double ww = 0.001;
constexpr double qq = 0.;

// equation: du/dt = v * df(u)/dx + w * d2g(u)/dx2 + q

std::vector<bound_cond_vecs::BoundCondVec<double>> v(std::vector<bound_cond_vecs::BoundCondVec<double>> u, double t, bound_cond_vecs::BoundCondVec<double> x) {	
	std::vector<bound_cond_vecs::BoundCondVec<double>> vu;
	for(int ii = 0; ii < u.size(); ii++) {
		bound_cond_vecs::BoundCondVec<double> vu_tmp(x.len(), x.getMode());
		
		for(int jj = 0; jj < x.len(); jj++) vu_tmp[jj] = vv * u[ii][jj];
		
		vu.push_back(vu_tmp);
	}
	return vu;
}

std::vector<bound_cond_vecs::BoundCondVec<double>> w(std::vector<bound_cond_vecs::BoundCondVec<double>> u, double t, bound_cond_vecs::BoundCondVec<double> x) {
	std::vector<bound_cond_vecs::BoundCondVec<double>> wu;
	for(int ii = 0; ii < u.size(); ii++) {
		bound_cond_vecs::BoundCondVec<double> wu_tmp(x.len(), x.getMode());
		
		for(int jj = 0; jj < x.len(); jj++) wu_tmp[jj] = ww;
		
		wu.push_back(wu_tmp);
	}
	return wu;
}

std::vector<bound_cond_vecs::BoundCondVec<double>> q(std::vector<bound_cond_vecs::BoundCondVec<double>> u, double t, bound_cond_vecs::BoundCondVec<double> x) {
	std::vector<bound_cond_vecs::BoundCondVec<double>> qu;
	for(int ii = 0; ii < u.size(); ii++) {
		bound_cond_vecs::BoundCondVec<double> qu_tmp(x.len(), x.getMode());
		
		for(int jj = 0; jj < x.len(); jj++) qu_tmp[jj] = qq;
		
		qu.push_back(qu_tmp);
	}
	return qu;
}

std::vector<bound_cond_vecs::BoundCondVec<double>> f(std::vector<bound_cond_vecs::BoundCondVec<double>> u, double t, bound_cond_vecs::BoundCondVec<double> x) {
	std::vector<bound_cond_vecs::BoundCondVec<double>> fu;
	for(int ii = 0; ii < u.size(); ii++) {
		bound_cond_vecs::BoundCondVec<double> fu_tmp(x.len(), x.getMode());
		
		for(int jj = 0; jj < x.len(); jj++) fu_tmp[jj] = u[ii][jj];
		
		fu.push_back(fu_tmp);
	}
	return fu;
}

std::vector<bound_cond_vecs::BoundCondVec<double>> g(std::vector<bound_cond_vecs::BoundCondVec<double>> u, double t, bound_cond_vecs::BoundCondVec<double> x) {
	std::vector<bound_cond_vecs::BoundCondVec<double>> gu;
	for(int ii = 0; ii < u.size(); ii++) {
		bound_cond_vecs::BoundCondVec<double> gu_tmp(x.len(), x.getMode());
		
		for(int jj = 0; jj < x.len(); jj++) gu_tmp[jj] = u[ii][jj];
		
		gu.push_back(gu_tmp);
	}
	return gu;
}

int main() {
	myStyle();
	
	double t0 = 0., dt = 0.001;
	int nsteps = 300, nx = 500;
	double x0 = -0.5, dx = 1. / nx;
	bound_cond_vecs::BoundCondVec<double> x = integrators::linspace(x0, x0 + nx * dx, nx, PERIODIC_BC);
	std::vector<bound_cond_vecs::BoundCondVec<double>> u0;
	bound_cond_vecs::BoundCondVec<double> u0_tmp(nx, x.getMode());
	for(int ii = 0; ii < nx; ii++) {
		
		// u0_tmp[ii] = sin(2. * M_PI * 1. * x[ii]); // cosine
		u0_tmp[ii] = exp(- (x[ii]) * (x[ii]) / (2. * 0.01 * 0.01)); // gaussian
		
	}
	u0.push_back(u0_tmp);
	
	std::cout << "v*dt/dx = " << vv * dt / dx << std::endl;
	std::cout << "w*dt/dx^2 = " << ww * dt / (dx * dx) << std::endl;
	
	auto u1 = integrators::FTCS(t0, dt, nsteps, x, u0, f, derivators::symm_derive, g, derivators::symm_derive_2, v, w, q);
	auto u2 = integrators::Lax(t0, dt, nsteps, x, u0, f, derivators::symm_derive, g, derivators::symm_derive_2, v, w, q);	
	auto u3 = integrators::LeapFrog(t0, dt, nsteps, x, u0, f, derivators::symm_derive, g, derivators::symm_derive_2, v, w, q);
	auto u4 = integrators::LaxWendroff(t0, dt, nsteps, x, u0, f, derivators::fwd_derive, g, derivators::symm_derive_2, v, w, q);
	
	fs::current_path(fs::current_path() / "measures");
	double minu[1] = {-2.}, maxu[2] = {2.};
	// double minu[2] = {0.}, maxu[2] = {1.};
	// waveplots::plot(u1, t0, dt, nsteps, x0, dx, nx, "first_FTCS_test_SURF", SURF_PLOT, minu, maxu);
	// waveplots::plot(u1, t0, dt, nsteps, x0, dx, nx, "first_FTCS_test_CONT", CONT_PLOT, minu, maxu);
	// waveplots::plot(u1, t0, dt, nsteps, x0, dx, nx, "first_FTCS_test_COLZ", COLZ_PLOT, minu, maxu);
	// waveplots::plot(u2, t0, dt, nsteps, x0, dx, nx, "first_Lax_test_SURF", SURF_PLOT, minu, maxu);
	// waveplots::plot(u2, t0, dt, nsteps, x0, dx, nx, "first_Lax_test_CONT", CONT_PLOT, minu, maxu);
	// waveplots::plot(u2, t0, dt, nsteps, x0, dx, nx, "first_Lax_test_COLZ", COLZ_PLOT, minu, maxu);
	// waveplots::plotFFT(u2, t0, dt, nsteps, "first_Lax_test_FFT_SURF", SURF_PLOT);
	// waveplots::plot(u3, t0, dt, nsteps, x0, dx, nx, "first_LeapFrog_test_SURF", SURF_PLOT, minu, maxu);
	// waveplots::plot(u3, t0, dt, nsteps, x0, dx, nx, "first_LeapFrog_CONT", CONT_PLOT, minu, maxu);
	// waveplots::plot(u3, t0, dt, nsteps, x0, dx, nx, "first_LeapFrog_COLZ", COLZ_PLOT, minu, maxu);
	// waveplots::plotFFT(u3, t0, dt, nsteps, "first_LeapFrog_test_FFT_SURF", SURF_PLOT);
	waveplots::plot(u4, t0, dt, nsteps, x0, dx, nx, "first_LaxWendroff_test_SURF", SURF_PLOT, minu, maxu);
	// waveplots::plot(u4, t0, dt, nsteps, x0, dx, nx, "first_LaxWendroff_test_CONT", CONT_PLOT, minu, maxu);
	waveplots::plot(u4, t0, dt, nsteps, x0, dx, nx, "first_LaxWendroff_test_COLZ", COLZ_PLOT, minu, maxu);
	waveplots::plotFFT(u4, t0, dt, nsteps, "first_LaxWendroff_test_FFT_SURF", SURF_PLOT);
}