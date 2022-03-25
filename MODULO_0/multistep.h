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
	
	// BaseMultistep ---------------------------------------------------
	
	class BaseMultistep {
		protected:
			int n_a;
			int n_b;
			int p;
			int n_dim;
		
		public:
			std::vector<double> a;
			std::vector<double> b;
			int order;
			int n_max;
			std::vector<double> x;
			std::vector<std::vector<double>> y;
			std::vector<std::vector<double>> y1;
			
			BaseMultistep() {;}
			
			BaseMultistep(std::vector<double> & aa, std::vector<double> & bb, int o) {
				for(auto ai: aa) a.push_back(ai);
				for(auto bi: bb) b.push_back(bi);
				n_a = a.size();
				n_b = b.size();
				p = MAX(n_a, n_b) - 1;
				order = o;
			}
			
			virtual void checkSetupConsistency() = 0;
			
			virtual void print() = 0;
			
			void clear() {
				n_max = 0;
				n_dim = 0;
				x.clear();
				y.clear();
				y1.clear();
			}
			
			void setInit(double init_x, std::vector<double> & init_y) {
				n_max = 1;
				n_dim = init_y.size();
				x.reserve(p + 1);
				x.push_back(init_x);
				y.reserve(p + 1);
				y.push_back(init_y);
				y1.reserve(p + 1);
			}
			
			virtual void step(double h) = 0;
			
			std::vector<std::vector<double>> run(int nsteps, double h) {
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
	
	// ExplicitMultistep ---------------------------------------------------
	
	template <typename... Parameters>
	class ExplicitMultistep : public BaseMultistep{
		protected:
			std::vector<double> (*eq)(std::vector<double>, double, Parameters...);
			std::tuple<Parameters...> pars;
			rungekutta::ExplicitRungeKutta<Parameters...> rk;
			
		public:
			
			ExplicitMultistep() {;}
			
			ExplicitMultistep(std::vector<double> & aa, std::vector<double> & bb, int o, rungekutta::ExplicitRungeKutta<Parameters...> ex_rk) : BaseMultistep(aa, bb, o) {
				rk = ex_rk;
				checkSetupConsistency();
			}
			
			void checkSetupConsistency() override {
				double asum = -1.;
				for(auto& ii: a) asum += ii;
				assert(asum * asum <= MS_ASSERT_EPS * MS_ASSERT_EPS);
				double bsum = 0., kasum = 0.;
				for(auto& ii: b) bsum += ii;
				for(int ii = 0; ii < n_a; ii++) kasum += (p - ii) * a[ii];
				double temp = bsum + kasum - p - 1;
				assert(temp * temp <= MS_ASSERT_EPS * MS_ASSERT_EPS);
			}
			
			void print() override {
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
				BaseMultistep::clear();
				rk.clear();
			}
			
			void setEquation(std::vector<double> (*equation)(std::vector<double>, double, Parameters...), Parameters... parameters) {
				eq = equation;
				pars = std::tie(parameters...);
				rk.setEquation(equation, parameters...);
			}
			
			void setInit(double init_x, std::vector<double> & init_y) {
				y1.push_back(std::apply(eq, std::tuple_cat(std::make_tuple(init_y, init_x), pars)));
				BaseMultistep::setInit(init_x, init_y);
				rk.setInit(init_x, init_y);
			}
			
			void step(double h) override {
				std::vector<double> sum_(n_dim, 0.);
				for(int i = 0; i < n_dim; i++) {
					for(int j = 0; j < n_a; j++) sum_[i] += a[j] * y[n_max - 1 - j][i];
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
				
				return BaseMultistep::run(nsteps, h);
			}
	};
	
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
	
	template <typename... Parameters>
	ExplicitMultistep<Parameters...> AB1() {
		std::vector<double> a = {1.};
		std::vector<double> b = {1.};
		return ExplicitMultistep<Parameters...>(a, b, 1, rungekutta::RK4<Parameters...>());
	}
	
	template <typename... Parameters>
	ExplicitMultistep<Parameters...> AB2() {
		std::vector<double> a = {1.};
		std::vector<double> b = {3./2, -1./2};
		return ExplicitMultistep<Parameters...>(a, b, 2, rungekutta::RK4<Parameters...>());
	}
	
	template <typename... Parameters>
	ExplicitMultistep<Parameters...> AB3() {
		std::vector<double> a = {1.};
		std::vector<double> b = {23./12, -16./12, 5./12};
		return ExplicitMultistep<Parameters...>(a, b, 3, rungekutta::RK4<Parameters...>());
	}
	
	template <typename... Parameters>
	ExplicitMultistep<Parameters...> AB4() {
		std::vector<double> a = {1.};
		std::vector<double> b = {55./24, -59./24, 37./24, -9./24};
		return ExplicitMultistep<Parameters...>(a, b, 4, rungekutta::RK4<Parameters...>());
	}
	
	template <typename... Parameters>
	ExplicitMultistep<Parameters...> Midpoint() {
		std::vector<double> a = {0., 1.};
		std::vector<double> b = {2.};
		return ExplicitMultistep<Parameters...>(a, b, 2, rungekutta::RK4<Parameters...>());
	}
	
	template <typename... Parameters>
	ExplicitMultistep<Parameters...> Euler() {
		return AB1<Parameters...>();
	}
	
	// ImplicitMultistep ---------------------------------------------------
	
	template <typename... Parameters>
	class ImplicitMultistep: public BaseMultistep {
		protected:
			int n_a_ex;
			int n_b_ex;
			std::vector<double> (*eq)(std::vector<double>, double, Parameters...);
			std::tuple<Parameters...> pars;
			rungekutta::ExplicitRungeKutta<Parameters...> rk;
			
		public:
			double b_1;
			std::vector<double> a_ex;
			std::vector<double> b_ex;
			
			ImplicitMultistep() {;}
		
			ImplicitMultistep(std::vector<double> & aa, std::vector<double> & bb, int o, rungekutta::ExplicitRungeKutta<Parameters...> ex_rk, double bb_1, ExplicitMultistep<Parameters...> ex_ms): BaseMultistep(aa, bb, o) {
				rk = ex_rk;
				b_1 = bb_1;
				for(auto ai: ex_ms.a) a_ex.push_back(ai);
				for(auto bi: ex_ms.b) b_ex.push_back(bi);
				n_a_ex = a_ex.size();
				n_b_ex = b_ex.size();
				checkSetupConsistency();
			}
			
			void checkSetupConsistency() override {
				double asum = -1.;
				for(auto& ii: a) asum += ii;
				assert(asum * asum <= MS_ASSERT_EPS * MS_ASSERT_EPS);
				double bsum = b_1, kasum = 0.;
				for(auto& ii: b) bsum += ii;
				for(int ii = 0; ii < n_a; ii++) kasum += (p - ii) * a[ii];
				double temp = bsum + kasum - p - 1;
				assert(temp * temp <= MS_ASSERT_EPS * MS_ASSERT_EPS);
			}
			
			void print() override {
				std::cout << "   ";
				for(int i = -1; i <= p; i++) std::cout << " \t" << i;
				std::cout << std::endl;
				std::cout << "a :" << "   \t";
				for(auto& ii: a) std::cout << " \t" << ii;
				std::cout << std::endl;
				std::cout << "b :" << " \t" << b_1;
				for(auto& ii: b) std::cout << " \t" << ii;
				std::cout << std::endl << std::endl;
			}
			
			void clear() {
				BaseMultistep::clear();
				rk.clear();
			}
			
			void setEquation(std::vector<double> (*equation)(std::vector<double>, double, Parameters...), Parameters... parameters) {
				eq = equation;
				pars = std::tie(parameters...);
				rk.setEquation(equation, parameters...);
			}
			
			void setInit(double init_x, std::vector<double> & init_y) {
				y1.push_back(std::apply(eq, std::tuple_cat(std::make_tuple(init_y, init_x), pars)));
				BaseMultistep::setInit(init_x, init_y);
				rk.setInit(init_x, init_y);
			}
			
			void step(double h) override {
				double pred_x = x[n_max - 1] + h;
				std::vector<double> pred_y(n_dim, 0.);
				for(int i = 0; i < n_dim; i++) {
					for(int j = 0; j < n_a_ex; j++) pred_y[i] += a_ex[j] * y[n_max - 1 - j][i];
					for(int j = 0; j < n_b_ex; j++) pred_y[i] += h * b_ex[j] * y1[n_max - 1 - j][i];
				}
				std::vector<double> pred_y1 = std::apply(eq, std::tuple_cat(std::make_tuple(pred_y, pred_x), pars));
				std::vector<double> sum_(n_dim, 0.);
				for(int i = 0; i < n_dim; i++) {
					sum_[i] += h * b_1 * pred_y1[i];
					for(int j = 0; j < n_a; j++) sum_[i] += a[j] * y[n_max - 1 - j][i];
					for(int j = 0; j < n_b; j++) sum_[i] += h * b[j] * y1[n_max - 1 - j][i];
				}
				auto tmp_y1 = std::apply(eq, std::tuple_cat(std::make_tuple(sum_, pred_x), pars));
				x.push_back(pred_x);
				y.push_back(sum_);
				y1.push_back(tmp_y1);
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
				
				return BaseMultistep::run(nsteps, h);
			}
	};
	
	template <typename... Parameters>
	ImplicitMultistep<Parameters...> AM0() {
		std::vector<double> a = {1.};
		std::vector<double> b = {};
		double b_1 = 1.;
		return ImplicitMultistep<Parameters...>(a, b, 1, rungekutta::RK4<Parameters...>(), b_1, multistep::AB2<Parameters...>());
	}
	
	template <typename... Parameters>
	ImplicitMultistep<Parameters...> AM1() {
		std::vector<double> a = {1.};
		std::vector<double> b = {1./2};
		double b_1 = 1./2;
		return ImplicitMultistep<Parameters...>(a, b, 2, rungekutta::RK4<Parameters...>(), b_1, multistep::AB2<Parameters...>());
	}
	
	template <typename... Parameters>
	ImplicitMultistep<Parameters...> AM2() {
		std::vector<double> a = {1.};
		std::vector<double> b = {8./12, -1./12};
		double b_1 = 5./12;
		return ImplicitMultistep<Parameters...>(a, b, 3, rungekutta::RK4<Parameters...>(), b_1, multistep::AB2<Parameters...>());
	}
	
	template <typename... Parameters>
	ImplicitMultistep<Parameters...> AM3() {
		std::vector<double> a = {1.};
		std::vector<double> b = {19./24, -5./24, 1./24};
		double b_1 = 9./24;
		return ImplicitMultistep<Parameters...>(a, b, 4, rungekutta::RK4<Parameters...>(), b_1, multistep::AB2<Parameters...>());
	}
	
	template <typename... Parameters>
	ImplicitMultistep<Parameters...> Trapezioidal() {
		return AM1<Parameters...>();
	}
	
	template <typename... Parameters>
	void runWithRichardsonError(ImplicitMultistep<Parameters...> ms, double init_x, std::vector<double> & init_y, int nsteps, double h, std::vector<double>& x, std::vector<std::vector<double>>& y, std::vector<std::vector<double>>& y1, std::vector<std::vector<double>>& dy) {
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