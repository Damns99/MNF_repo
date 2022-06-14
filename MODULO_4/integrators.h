#ifndef INTEGRATORS_H
#define INTEGRATORS_H

#include <math.h>
#include <vector>
#include <iostream>

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
		std::vector<DblBcv> (*)(std::vector<DblBcv>, double, DblBcv), derivators::Derivator);
	
	namespace utilities {
		// equation: du/dt = df(u)/dx + d2g(u)/dx2
		
		std::vector<std::vector<DblBcv>> FTCS_function(
			double t0, double dt, double nsteps, const DblBcv & x, const std::vector<DblBcv> & u0, 
			std::vector<DblBcv> (*f)(std::vector<DblBcv>, double, DblBcv), derivators::Derivator d1, 
			std::vector<DblBcv> (*g)(std::vector<DblBcv>, double, DblBcv), derivators::Derivator d2) 
		{
			double dx = x[1] - x[0];
			int ndim = u0.size();
			std::vector<std::vector<DblBcv>> u;
			u.reserve(nsteps + 1);
			u.push_back(u0);
			std::vector<double> t;
			t.reserve(nsteps + 1);
			t.push_back(t0);
			
			// for(int qq = 0; qq < x.len(); qq++) std::cout << "x[" << qq << "] = " << x[qq] << std::endl;
			
			for(int ii = 1; ii <= nsteps; ii++) {
				std::vector<DblBcv> new_u;
				new_u.reserve(ndim);
				std::vector<DblBcv> fu = f(u[ii - 1], t0 + (ii - 1) * dt, x);
				std::vector<DblBcv> gu = g(u[ii - 1], t0 + (ii - 1) * dt, x);
				for(int jj = 0; jj < ndim; jj++) {
					DblBcv d1fu = d1(fu[jj], dx), d2gu = d2(gu[jj], dx), tmp_u(x.len(), x.getMode());
					for(int kk = 0; kk < x.len(); kk++) tmp_u[kk] = u[ii- 1][jj][kk] + dt * (d1fu[kk] + d2gu[kk]);
					new_u.push_back(tmp_u);
				}					
				u.push_back(new_u);
				t.push_back(t0 + ii * dt);
			}
			
			return u;
		}
		
	}
	
	Integrator FTCS = &(utilities::FTCS_function);
	
}

#endif