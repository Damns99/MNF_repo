#ifndef RUNGE_KUTTA_H
#define RUNGE_KUTTA_H

#include <vector>
#include <tuple>
#include "text_io.h"
#include <iostream>

namespace rungekutta {
	
	template <typename... Parameters>
	class ExplicitRungeKutta {
		private:
			std::vector<std::vector<double>> a;
			std::vector<double> b;
			std::vector<double> c;
			int s;
			int n_max;
			int n_dim;
			std::vector<double> (*eq)(std::vector<double>, double, Parameters...);
			std::tuple<Parameters...> pars;
			
		public:
			std::vector<double> x;
			std::vector<std::vector<double>> y;
			
			ExplicitRungeKutta(std::vector<std::vector<double>> & aa, std::vector<double> & bb, std::vector<double> & cc) {
				for(auto ai: aa) a.push_back(ai);
				for(auto bi: bb) b.push_back(bi);
				for(auto ci: cc) c.push_back(ci);
				s = b.size();
			}
			
			void printTableau() {
				for(int i = 0; i < s; i++) {
					std::cout << c[i] << " \t|";
					for(int j = 0; j < i; j++) std::cout << " \t" << a[i-1][j];
					std::cout << std::endl;
				}
				std::cout << " \t|";
				for(int i = 0; i < s; i++) std::cout << " \t" << b[i];
				std::cout << std::endl;
			}
			
			void clear() {
				n_max = 0;
				n_dim = 0;
				x.clear();
				y.clear();
			}
			
			void setEquation(std::vector<double> (*equation)(std::vector<double>, double, Parameters...), Parameters... parameters) {
				eq = equation;
				pars = std::tie(parameters...);
			}
			
			int setInit(double init_x, std::vector<double> & init_y) {
				n_max = 1;
				n_dim = init_y.size();
				x.reserve(n_max);
				x.push_back(init_x);
				y.reserve(n_max);
				y.push_back(init_y);
				return 0;
			}
			
			void step(double h) {
				std::vector<double> sum_;
				std::vector<std::vector<double>> k;
				for(int j = 0; j < s; j++) {
					std::vector<double> tmpy_ = y[n_max - 1];
					for(int jj = 0; jj < j; jj++) for(int i = 0; i < n_dim; i++) {
						tmpy_[i] += h * a[j-1][jj] * k[jj][i];
					}
					k.push_back(std::apply(eq, std::tuple_cat(std::make_tuple(tmpy_, x[n_max - 1] + c[j] * h), pars)));
				}
				for(int i = 0; i < n_dim; i++) {
					double tmp_ = y[n_max - 1][i];
					for(int j = 0; j < s; j++) tmp_ += b[j] * k[j][i] * h;
					sum_.push_back(tmp_);
				}
				x.push_back(x[n_max - 1] + h);
				y.push_back(sum_);
				n_max += 1;
			}
			
			std::vector<std::vector<double>> run(int nsteps, double h) {
				x.reserve(nsteps + n_max);
				y.reserve(nsteps + n_max);
				for(int i = 0; i < nsteps; i++) step(h);
				return y;
			}
			
			void printToFile(std::string filename = "rk.txt") {
				textIo::textOut(filename, '\t', '#', "x \t y", n_max, false, x, y);
			}
	};
	
	template <typename... Parameters>
	ExplicitRungeKutta<Parameters...> RK4() {
		std::vector<std::vector<double>> a = {{0.5}, {0., 0.5}, {0., 0., 1.}};
		std::vector<double> b = {1./6, 1./3, 1./3, 1./6};
		std::vector<double> c = {0., 0.5, 0.5, 1.};
		return ExplicitRungeKutta<Parameters...>(a, b, c);
	}
	
	template <typename... Parameters>
	void runWithRichardsonError(ExplicitRungeKutta<Parameters...> rk, double init_x, std::vector<double> & init_y, int nsteps, double h, int order, std::vector<double>& x, std::vector<std::vector<double>>& y, std::vector<std::vector<double>>& dy) {
		std::cout << "Running with Richardson Error Estimation." << std::endl;
		std::cout << "It may take a while..." << std::endl;
		rk.clear();
		rk.setInit(init_x, init_y);
		auto yh2 = rk.run(2 * nsteps, h / 2);
		rk.clear();
		rk.setInit(init_x, init_y);
		y = rk.run(nsteps, h);
		x = rk.x;
		std::vector<std::vector<double>> yhh;
		for(int i = 0; i < nsteps + 1; i++) {
			std::vector<double> tmp_rich;
			for(int j = 0; j < y[0].size(); j++) tmp_rich.push_back((yh2[2 * i][j] - y[i][j]) / (pow(2, order) - 1));
			dy.push_back(tmp_rich);
		}
	}
	
}

#endif