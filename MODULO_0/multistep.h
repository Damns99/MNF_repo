#ifndef MULTISTEP_H
#define MULTISTEP_H

#include <vector>
#include <tuple>
#include <iostream>
#include <cassert>

#include "runge_kutta.h"
#include "text_io.h"

#define MAX(X, Y) X >= Y ? X : Y

#define MS_ASSERT_EPS 5e-16

namespace multistep {
	
	template <typename... Parameters>
	class ExplicitMultistep {
		private:
			std::vector<double> a;
			std::vector<double> b;
			int n_a;
			int n_b;
			int p;
			int n_max;
			int n_dim;
			std::vector<double> (*eq)(std::vector<double>, double, Parameters...);
			std::tuple<Parameters...> pars;
			rungekutta::ExplicitRungeKutta<Parameters...> rk = rungekutta::RK4<Parameters...>();
			
		public:
			int order;
			std::vector<double> x;
			std::vector<std::vector<double>> y;
			std::vector<std::vector<double>> y1;
			
			ExplicitMultistep(std::vector<double> & aa, std::vector<double> & bb, int o) {
				for(auto ai: aa) a.push_back(ai);
				for(auto bi: bb) b.push_back(bi);
				n_a = a.size();
				n_b = b.size();
				p = MAX(n_a, n_b) - 1;
				order = o;
				checkSetupConsistency();
			}
			
			void checkSetupConsistency() {
				double asum = -1.;
				for(auto& ii: a) asum += ii;
				assert(asum * asum <= MS_ASSERT_EPS * MS_ASSERT_EPS);
				double bsum = 0., kasum = 0.;
				for(auto& ii: b) bsum += ii;
				for(int ii = 0; ii < n_a; ii++) kasum += (p - ii) * a[ii];
				double temp = bsum + kasum - p - 1;
				assert(temp * temp <= MS_ASSERT_EPS * MS_ASSERT_EPS);
			}
			
			void print() {
				std::cout << "   ";
				for(int i = 0; i <= p; i++) std::cout << " \t" << i;
				std::cout << std::endl;
				std::cout << "a :";
				for(auto& ii: a) std::cout << " \t" << ii;
				std::cout << std::endl;
				std::cout << "b :";
				for(auto& ii: b) std::cout << " \t" << ii;
				std::cout << std::endl << std::endl;
			}
			
			void clear() {
				n_max = 0;
				n_dim = 0;
				x.clear();
				y.clear();
				y1.clear();
				rk.clear();
			}
			
			void setEquation(std::vector<double> (*equation)(std::vector<double>, double, Parameters...), Parameters... parameters) {
				eq = equation;
				pars = std::tie(parameters...);
				rk.setEquation(equation, parameters...);
			}
			
			void setInit(double init_x, std::vector<double> & init_y) {
				n_max = 1;
				n_dim = init_y.size();
				x.reserve(p + 1);
				x.push_back(init_x);
				y.reserve(p + 1);
				y.push_back(init_y);
				y1.reserve(p + 1);
				y1.push_back(std::apply(eq, std::tuple_cat(std::make_tuple(init_y, init_x), pars)));
				rk.setInit(init_x, init_y);
			}
			
			void step(double h) {
				std::vector<double> sum_(n_dim, 0.);
				for(int i = 0; i < n_dim; i++) {
					for(int j = 0; j < n_a; j++) {
						sum_[i] += a[j] * y[n_max - 1 - j][i];
					}
					for(int j = 0; j < n_b; j++) sum_[i] += h * b[j] * y1[n_max - 1 - j][i];
				}
				x.push_back(x[n_max - 1] + h);
				y.push_back(sum_);
				y1.push_back(std::apply(eq, std::tuple_cat(std::make_tuple(y[n_max], x[n_max]), pars)));
				n_max += 1;
			}
			
			std::vector<std::vector<double>> run(int nsteps, double h) {
				std::cout << "Calculating first " << p << " points with RK4" << std::endl;
				auto y_temp = rk.run(p, h);
				auto x_temp = rk.x;
				bool skip = true;
				for(auto& ii: x_temp) {
					if(skip) skip = false;
					else x.push_back(ii);
				}
				skip = true;
				for(auto& ii: y_temp) {
					if(skip) skip = false;
					else {
						y.push_back(ii);
						y1.push_back(std::apply(eq, std::tuple_cat(std::make_tuple(y[n_max], x[n_max]), pars)));
						n_max++;
					}
				}
				nsteps -= p;
				
				std::cout << "Multistep: " << nsteps << " points with step " << h << std::endl;
				x.reserve(nsteps + n_max);
				y.reserve(nsteps + n_max);
				y1.reserve(nsteps + n_max);
				for(int i = 0; i < nsteps; i++) step(h);
				return y;
			}
			
			void printToFile(std::string filename = "ms.txt") {
				textIo::textOut(filename, '\t', '#', "x \t y \t y1", n_max, false, x, y, y1);
			}
	};
	
	template <typename... Parameters>
	ExplicitMultistep<Parameters...> AB2() {
		std::vector<double> a = {1.};
		std::vector<double> b = {3./2, -1./2};
		return ExplicitMultistep<Parameters...>(a, b, 2);
	}
	
	template <typename... Parameters>
	void runWithRichardsonError(ExplicitMultistep<Parameters...> ms, double init_x, std::vector<double> & init_y, int nsteps, double h, std::vector<double>& x, std::vector<std::vector<double>>& y, std::vector<std::vector<double>>& y1, std::vector<std::vector<double>>& dy) {
		std::cout << "Running with Richardson Error Estimation. It may take a while..." << std::endl;
		ms.clear();
		ms.setInit(init_x, init_y);
		auto yh2 = ms.run(2 * nsteps, h / 2);
		ms.clear();
		ms.setInit(init_x, init_y);
		y = ms.run(nsteps, h);
		x = ms.x;
		y1 = ms.y1;
		std::vector<std::vector<double>> yhh;
		for(int i = 0; i < nsteps + 1; i++) {
			std::vector<double> tmp_rich;
			for(int j = 0; j < y[0].size(); j++) tmp_rich.push_back((yh2[2 * i][j] - y[i][j]) / (pow(2, ms.order) - 1));
			dy.push_back(tmp_rich);
		}
	}
}

#endif