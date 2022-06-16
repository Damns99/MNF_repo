#ifndef INTEGRATORS_H
#define INTEGRATORS_H

#include <math.h>
#include <vector>

#include "bound_cond_vecs.h"
#include "derivators.h"

namespace integrators {
	
	using DblBcv = bound_cond_vecs::BoundCondVec<double>;
	
	DblBcv linspace(double x0, double x1, double n = 100, int m = PERIODIC_BC) {
		DblBcv res(n, m);
		double t = (x1 - x0) / n;
		for(int ii = 0; ii < n; ii++) res[ii] = x0 + t * ii;
		return res;
	}
	
	typedef std::vector<std::vector<DblBcv>> (*Integrator)(double, double, double, const DblBcv &, const std::vector<DblBcv> &, 
		std::vector<DblBcv> (*)(std::vector<DblBcv>, double, DblBcv), derivators::Derivator, 
		std::vector<DblBcv> (*)(std::vector<DblBcv>, double, DblBcv), derivators::Derivator,
		std::vector<DblBcv> (*)(std::vector<DblBcv>, double, DblBcv), std::vector<DblBcv> (*)(std::vector<DblBcv>, double, DblBcv),
		std::vector<DblBcv> (*)(std::vector<DblBcv>, double, DblBcv));
	
	namespace utilities {
		// equation: du/dt = v * df(u)/dx + w * d2g(u)/dx2 + q
		
		std::vector<std::vector<DblBcv>> FTCS_function(
			double t0, double dt, double nsteps, const DblBcv & x, const std::vector<DblBcv> & u0, 
			std::vector<DblBcv> (*f)(std::vector<DblBcv>, double, DblBcv), derivators::Derivator d1, 
			std::vector<DblBcv> (*g)(std::vector<DblBcv>, double, DblBcv), derivators::Derivator d2,
			std::vector<DblBcv> (*v)(std::vector<DblBcv>, double, DblBcv), std::vector<DblBcv> (*w)(std::vector<DblBcv>, double, DblBcv),
			std::vector<DblBcv> (*q)(std::vector<DblBcv>, double, DblBcv)) 
		{
			double dx = x[1] - x[0];
			int ndim = u0.size();
			std::vector<std::vector<DblBcv>> u;
			u.reserve(nsteps + 1);
			u.push_back(u0);
			std::vector<double> t;
			t.reserve(nsteps + 1);
			t.push_back(t0);
			
			for(int ii = 1; ii <= nsteps; ii++) {
				std::vector<DblBcv> new_u;
				new_u.reserve(ndim);
				std::vector<DblBcv> fu = f(u[ii - 1], t0 + (ii - 1) * dt, x);
				std::vector<DblBcv> gu = g(u[ii - 1], t0 + (ii - 1) * dt, x);
				std::vector<DblBcv> vu = v(u[ii - 1], t0 + (ii - 1) * dt, x);
				std::vector<DblBcv> wu = w(u[ii - 1], t0 + (ii - 1) * dt, x);
				std::vector<DblBcv> qu = q(u[ii - 1], t0 + (ii - 1) * dt, x);
				for(int jj = 0; jj < ndim; jj++) {
					DblBcv d1fu = d1(fu[jj], dx), d2gu = d2(gu[jj], dx), tmp_u(x.len(), x.getMode());
					for(int kk = 0; kk < x.len(); kk++) tmp_u[kk] = u[ii- 1][jj][kk] + dt * (vu[jj][kk] * d1fu[kk] + wu[jj][kk] * d2gu[kk] + qu[jj][kk]);
					new_u.push_back(tmp_u);
				}
				u.push_back(new_u);
				t.push_back(t0 + ii * dt);
			}
			
			return u;
		}
		
		std::vector<std::vector<DblBcv>> Lax_function(
			double t0, double dt, double nsteps, const DblBcv & x, const std::vector<DblBcv> & u0, 
			std::vector<DblBcv> (*f)(std::vector<DblBcv>, double, DblBcv), derivators::Derivator d1, 
			std::vector<DblBcv> (*g)(std::vector<DblBcv>, double, DblBcv), derivators::Derivator d2,
			std::vector<DblBcv> (*v)(std::vector<DblBcv>, double, DblBcv), std::vector<DblBcv> (*w)(std::vector<DblBcv>, double, DblBcv),
			std::vector<DblBcv> (*q)(std::vector<DblBcv>, double, DblBcv)) 
		{
			double dx = x[1] - x[0];
			int ndim = u0.size();
			std::vector<std::vector<DblBcv>> u;
			u.reserve(nsteps + 1);
			u.push_back(u0);
			std::vector<double> t;
			t.reserve(nsteps + 1);
			t.push_back(t0);
			
			for(int ii = 1; ii <= nsteps; ii++) {
				std::vector<DblBcv> new_u;
				new_u.reserve(ndim);
				std::vector<DblBcv> fu = f(u[ii - 1], t0 + (ii - 1) * dt, x);
				std::vector<DblBcv> gu = g(u[ii - 1], t0 + (ii - 1) * dt, x);
				std::vector<DblBcv> vu = v(u[ii - 1], t0 + (ii - 1) * dt, x);
				std::vector<DblBcv> wu = w(u[ii - 1], t0 + (ii - 1) * dt, x);
				std::vector<DblBcv> qu = q(u[ii - 1], t0 + (ii - 1) * dt, x);
				for(int jj = 0; jj < ndim; jj++) {
					DblBcv d1fu = d1(fu[jj], dx), d2gu = d2(gu[jj], dx), tmp_u(x.len(), x.getMode());
					for(int kk = 0; kk < x.len(); kk++) tmp_u[kk] = (u[ii- 1][jj][kk - 1] + u[ii- 1][jj][kk + 1]) / 2 + dt * (vu[jj][kk] * d1fu[kk] + wu[jj][kk] * d2gu[kk] + qu[jj][kk]);
					new_u.push_back(tmp_u);
				}
				u.push_back(new_u);
				t.push_back(t0 + ii * dt);
			}
			
			return u;
		}
		
		std::vector<std::vector<DblBcv>> LeapFrog_function(
			double t0, double dt, double nsteps, const DblBcv & x, const std::vector<DblBcv> & u0, 
			std::vector<DblBcv> (*f)(std::vector<DblBcv>, double, DblBcv), derivators::Derivator d1, 
			std::vector<DblBcv> (*g)(std::vector<DblBcv>, double, DblBcv), derivators::Derivator d2,
			std::vector<DblBcv> (*v)(std::vector<DblBcv>, double, DblBcv), std::vector<DblBcv> (*w)(std::vector<DblBcv>, double, DblBcv),
			std::vector<DblBcv> (*q)(std::vector<DblBcv>, double, DblBcv)) 
		{
			double dx = x[1] - x[0];
			int ndim = u0.size();
			std::vector<std::vector<DblBcv>> u;
			u.reserve(nsteps + 1);
			u.push_back(u0);
			std::vector<double> t;
			t.reserve(nsteps + 1);
			t.push_back(t0);
			auto u1 = FTCS_function(t0, dt, 1, x, u0, f, d1, g, d2, v, w, q);
			u.push_back(u1[1]);
			t.push_back(t0 + dt);
			
			for(int ii = 2; ii <= nsteps; ii++) {
				std::vector<DblBcv> new_u;
				new_u.reserve(ndim);
				std::vector<DblBcv> fu = f(u[ii - 1], t0 + (ii - 1) * dt, x);
				std::vector<DblBcv> gu = g(u[ii - 1], t0 + (ii - 1) * dt, x);
				std::vector<DblBcv> vu = v(u[ii - 1], t0 + (ii - 1) * dt, x);
				std::vector<DblBcv> wu = w(u[ii - 1], t0 + (ii - 1) * dt, x);
				std::vector<DblBcv> qu = q(u[ii - 1], t0 + (ii - 1) * dt, x);
				for(int jj = 0; jj < ndim; jj++) {
					DblBcv d1fu = d1(fu[jj], dx), d2gu = d2(gu[jj], dx), tmp_u(x.len(), x.getMode());
					for(int kk = 0; kk < x.len(); kk++) tmp_u[kk] = u[ii- 2][jj][kk] + 2. * dt * (vu[jj][kk] * d1fu[kk] + wu[jj][kk] * d2gu[kk] + qu[jj][kk]);
					new_u.push_back(tmp_u);
				}
				u.push_back(new_u);
				t.push_back(t0 + ii * dt);
			}
			
			return u;
		}
		
		std::vector<std::vector<DblBcv>> LaxWendroff_function(
			double t0, double dt, double nsteps, const DblBcv & x, const std::vector<DblBcv> & u0, 
			std::vector<DblBcv> (*f)(std::vector<DblBcv>, double, DblBcv), derivators::Derivator d1, 
			std::vector<DblBcv> (*g)(std::vector<DblBcv>, double, DblBcv), derivators::Derivator d2,
			std::vector<DblBcv> (*v)(std::vector<DblBcv>, double, DblBcv), std::vector<DblBcv> (*w)(std::vector<DblBcv>, double, DblBcv),
			std::vector<DblBcv> (*q)(std::vector<DblBcv>, double, DblBcv)) 
		{
			double dx = x[1] - x[0];
			int ndim = u0.size();
			std::vector<std::vector<DblBcv>> u;
			u.reserve(nsteps + 1);
			u.push_back(u0);
			std::vector<double> t;
			t.reserve(nsteps + 1);
			t.push_back(t0);
			DblBcv x2 = x;
			for(int ii = 0; ii < x2.len(); ii++) x2[ii] += 0.5 * dx;
			
			for(int ii = 1; ii <= nsteps; ii++) {
				std::vector<DblBcv> new_u, new_u2;
				new_u.reserve(ndim);
				new_u2.reserve(ndim);
				std::vector<DblBcv> fu = f(u[ii - 1], t0 + (ii - 1) * dt, x);
				std::vector<DblBcv> gu = g(u[ii - 1], t0 + (ii - 1) * dt, x);
				std::vector<DblBcv> vu = v(u[ii - 1], t0 + (ii - 1) * dt, x);
				std::vector<DblBcv> wu = w(u[ii - 1], t0 + (ii - 1) * dt, x);
				std::vector<DblBcv> qu = q(u[ii - 1], t0 + (ii - 1) * dt, x);
				for(int jj = 0; jj < ndim; jj++) {
					DblBcv d1fu = d1(fu[jj], dx), d2gu = d2(gu[jj], dx), tmp_u(x.len(), x.getMode());
					for(int kk = 0; kk < x.len(); kk++) tmp_u[kk] = (u[ii- 1][jj][kk] + u[ii- 1][jj][kk + 1]) / 2 + dt / 2 * (vu[jj][kk] * d1fu[kk] + wu[jj][kk] * d2gu[kk] + qu[jj][kk]);
					new_u.push_back(tmp_u);
				}
				std::vector<DblBcv> fu2 = f(new_u, t0 + (ii - 0.5) * dt, x2);
				std::vector<DblBcv> gu2 = g(new_u, t0 + (ii - 0.5) * dt, x2);
				std::vector<DblBcv> vu2 = v(new_u, t0 + (ii - 0.5) * dt, x2);
				std::vector<DblBcv> wu2 = w(new_u, t0 + (ii - 0.5) * dt, x2);
				std::vector<DblBcv> qu2 = q(new_u, t0 + (ii - 0.5) * dt, x2);
				for(int jj = 0; jj < ndim; jj++) {
					DblBcv d1fu = d1(fu2[jj], dx), d2gu = d2(gu2[jj], dx), tmp_u(x.len(), x.getMode());
					for(int kk = 0; kk < x.len(); kk++) tmp_u[kk] = u[ii- 1][jj][kk] + dt * (vu2[jj][kk - 1] * d1fu[kk - 1] + wu2[jj][kk - 1] * d2gu[kk - 1] + qu2[jj][kk - 1]);
					new_u2.push_back(tmp_u);
				}
				
				u.push_back(new_u2);
				t.push_back(t0 + ii * dt);
			}
			
			return u;
		}
		
	}
	
	Integrator FTCS = &(utilities::FTCS_function);
	Integrator Lax = &(utilities::Lax_function);
	Integrator LeapFrog = &(utilities::LeapFrog_function);
	Integrator LaxWendroff = &(utilities::LaxWendroff_function);
	
}

#endif