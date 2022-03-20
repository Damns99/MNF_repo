#ifndef ROOT_FINDER_H
#define ROOT_FINDER_H

#define RF_EPS 5e-16
#define RF_MAX_LOOP 128
#define abs(X) X >= 0 ? X : -X

namespace rootFinder {
	
	class RootResult {
		public:
			double value;
			double error;
			RootResult(double v, double e) {
				value = v;
				error = e;
			}
	};
	
	template <typename... Parameters>
	double dummyDerivator(double f(double, Parameters...), double x, Parameters... pars) {
		double eps = 1e-8;
		return (f(x + eps, pars...) - f(x - eps, pars...)) / (2. * eps);
	}
	
	template <typename... Parameters>
	RootResult bisection(double f(double, Parameters...), double a, double b, Parameters... pars) {
		double fa = f(a, pars...);
		if (fa * fa < RF_EPS * RF_EPS) return RootResult(a, 0.);
		double fb = f(b, pars...);
		if (fb * fb < RF_EPS * RF_EPS) return RootResult(b, 0.);
		double c = (a + b) / 2.;
		double fc = f(c);
		while (fc * fc >= RF_EPS * RF_EPS) {
			if (fa * fc > 0) a = c;
			else b = c;
			c = (a + b) / 2.;
			fc = f(c);
		}
		return RootResult(c, (b - a) / 2.);
	}
	
	template <typename... Parameters>
	RootResult newton(double f(double, Parameters...), double f1(double, Parameters...), double x, Parameters... pars) {
		double fx = f(x, pars...);
		if (fx * fx < RF_EPS * RF_EPS) return RootResult(x, 0.);
		double delta = fx / f1(x, pars...);
		unsigned niter = 0;
		while (delta >= RF_EPS * RF_EPS && niter < RF_MAX_LOOP) {
			x = x - delta;
			fx = f(x, pars...);
			delta = fx / f1(x, pars...);
			niter++;
		}
		return RootResult(x, abs(delta));
	}
	
	template <typename... Parameters>
	RootResult newton(double f(double, Parameters...), double f1(double (*)(double, Parameters...), Parameters...), double x0, Parameters... pars) {
		double fx = f(x, pars...);
		if (fx * fx < RF_EPS * RF_EPS) return RootResult(x, 0.);
		double delta = fx / f1(f, x, pars...);
		unsigned niter = 0;
		while (delta >= RF_EPS * RF_EPS && niter < RF_MAX_LOOP) {
			x = x - delta;
			fx = f(x, pars...);
			delta = fx / f1(f, x, pars...);
			niter++;
		}
		return RootResult(x, abs(delta));
	}
	
	template <typename... Parameters>
	RootResult mixed(double f(double, Parameters...), double f1(double, Parameters...), double a, double b, Parameters... pars) {
		double fa = f(a, pars...);
		if (fa * fa < RF_EPS * RF_EPS) return RootResult(a, 0.);
		double fb = f(b, pars...);
		if (fb * fb < RF_EPS * RF_EPS) return RootResult(b, 0.);
		double c = (a + b) / 2.;
		double fc = f(c);
		for(ii = 0; ii < 4; ii++) {
			if (fa * fc > 0) a = c;
			else b = c;
			c = (a + b) / 2.;
			fc = f(c);
		}
		return newton(f, f1, c, pars...);
	}
	
}

#endif