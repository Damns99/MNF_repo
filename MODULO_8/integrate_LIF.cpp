#include "integrate_LIF.h"

// single neuron Leaky Integrate and Fire model with external input
// {
//    Cm * dV/dt = -g * (V - Vrest) + I(t)
//    if V >= Vth => V -> V_reset and I(t:t+tref) = 0
// }

// Forward Euler method
// dV(n)/dt = (V(n+1) - V(n)) / h
std::vector<double> int_lif::fwdEuler(double V0, double h, int N, std::vector<double>& I, double params[6]) {
	double Cm = params[0], g = params[1], Vrest = params[2], Vth = params[3], Vreset = params[4], tauref = params[5];
	double a = g/Cm*h;
	std::vector<double> V(N);
	int refr = 0;
	V[0] = V0;
	double new_V = V0;
	for(int n = 1; n < N; n++) {
		if(refr <= 0) new_V = (1-a) * new_V + a * Vrest + h/Cm * I[n-1];
		else {
			new_V = (1-a) * new_V + a * Vrest;
			refr--;
		}
		if(new_V >= Vth) {
			V[n] = Vreset;
			refr = int(tauref / h);
		}
		else V[n] = new_V;
	}
	return V;
}

// Backward Euler method
// dV(n)/dt = (V(n) - V(n-1)) / h
std::vector<double> int_lif::bwdEuler(double V0, double h, int N, std::vector<double>& I, double params[6]) {
	double Cm = params[0], g = params[1], Vrest = params[2], Vth = params[3], Vreset = params[4], tauref = params[5];
	double a = g/Cm*h;
	std::vector<double> V(N);
	int refr = 0;
	V[0] = V0;
	double new_V = V0;
	for(int n = 1; n < N; n++) {
		if(refr <= 0) new_V = (new_V + a * Vrest +  h/Cm * I[n]) / (1+a);
		else {
			new_V = (new_V + a * Vrest) / (1+a);
			refr--;
		}
		if(new_V >= Vth) {
			V[n] = Vreset;
			refr = int(tauref / h);
		}
		else V[n] = new_V;
	}
	return V;
}

// Heun method
// 
std::vector<double> int_lif::Heun(double V0, double h, int N, std::vector<double>& I, double params[6]) {
	double Cm = params[0], g = params[1], Vrest = params[2], Vth = params[3], Vreset = params[4], tauref = params[5];
	double a = g/Cm*h;
	std::vector<double> V(N);
	int refr = 0;
	V[0] = V0;
	double new_V = V0;
	for(int n = 1; n < N; n++) {
		if(refr <= 0) new_V = (a*a/2.-a+1) * new_V - (a*a/2.-a) * Vrest + (1-a)/2. * h/Cm * I[n-1] + 1/2. * h/Cm * I[n]; // (2-a)/2. * h/Cm * I[n-1];
		else {
			new_V = (a*a/2.-a+1) * new_V - (a*a/2.-a) * Vrest;
			refr--;
		}
		if(new_V >= Vth) {
			V[n] = Vreset;
			refr = int(tauref / h);
		}
		else V[n] = new_V;
	}
	return V;
}

// RungeKutta4 method
// I has 2N elements, one every h/2 seconds
std::vector<double> int_lif::RK4(double V0, double h, int N, std::vector<double>& I, double params[6]) {
	double Cm = params[0], g = params[1], Vrest = params[2], Vth = params[3], Vreset = params[4], tauref = params[5];
	double a = g/Cm*h;
	double c1 = (a*a*a*a/24.-a*a*a/6.+a*a/2.-a+1), c2 = 1-c1, c3 = (-a*a*a/24.+a*a/12.-a/6.+1/6.) * h/Cm, c4 = (a*a/12.-a/3.+2/3.) * h/Cm, c5 = 1/6. * h/Cm;
	std::vector<double> V(N);
	int refr = 0;
	V[0] = V0;
	double new_V = V0;
	for(int n = 1; n < N; n++) {
		if(refr <= 0) new_V = c1 * new_V + c2 * Vrest + c3 * I[2*n-2] + c4 * I[2*n-1] + c5 * I[2*n];
		else {
			new_V = c1 * new_V + c2 * Vrest;
			refr--;
		}
		if(new_V >= Vth) {
			V[n] = Vreset;
			refr = int(tauref / h);
		}
		else V[n] = new_V;
	}
	return V;
}

// currents

std::vector<double> int_lif::currents::square_wave(int n, int period, double amplitude, int phase, double offset) {
	std::vector<double> res(n);
	for(int ii = 0; ii < n; ii++) {
		if((ii + phase + period) % period < period/2) res[ii] = offset + amplitude;
		else res[ii] = offset - amplitude;
	}
	return res;
}

std::vector<double> int_lif::currents::constant_pulse(int n, int start, int duration, double amplitude) {
	std::vector<double> res(n, 0.);
	for(int ii = start; ii < start+duration && ii < n; ii++) {
		res[ii] = amplitude;
	}
	return res;
}

std::vector<double> int_lif::currents::pulse_train(int n, int start, int duration, double amplitude, int interval, int numpulses) {
	std::vector<double> res(n, 0.);
	for(int pp = 0; pp < numpulses; pp++) {
		for(int ii = start+pp*(duration+interval); ii < start+pp*(duration+interval)+duration && ii < n; ii++) {
			res[ii] = amplitude;
		}
	}
	return res;
}

// utils

std::vector<double> int_lif::utils::linspace(double x0, double x1, double n = 100) {
	std::vector<double> res(n);
	double dx = (x1 - x0) / n;
	for(int ii = 0; ii < n; ii++) res[ii] = x0 + dx * ii;
	return res;
}

std::vector<double> int_lif::utils::logspace(double e0, double e1, double n = 100) {
	std::vector<double> res(n);
	double de = (e1 - e0) / n;
	for(int ii = 0; ii < n; ii++) res[ii] = pow(10, e0 + de * ii);
	return res;
}

double int_lif::utils::mean(const std::vector<double>& vec) {
	double s = 0.;
	int length = vec.size();
	for(auto& ii: vec) s += ii / length;
	return s;
}

double int_lif::utils::stdev(const std::vector<double>& vec) {
	double s1 = 0., s2 = 0.;
	int length = vec.size();
	for(auto& ii: vec) {
		s1 += ii / length;
		s2 += (ii * ii) / length;
	}
	return sqrt(s2 - s1 * s1);
}

void int_lif::utils::statsCalcOnePass(const std::vector<double>& vec, double& mean, double& variance, double& skewness, double& kurtosis) {
	double xbar = 0., M2 = 0., M3 = 0., M4 = 0.;
	int N = 0;
	for(auto& ii:vec) {
		N++;
		double delta = (ii - xbar) / N;
		M4 += N*(N-1)*(N*N-3*N+3)*delta*delta*delta*delta + 6*M2*delta*delta - 4*M3*delta;
		M3 += N*(N-1)*(N-2)*delta*delta*delta - 3*M2*delta;
		M2 += N*(N-1)*delta*delta;
		xbar += delta;
	}
	double m2, m3, m4;
	m2 = M2 / (N-1);
	m3 = M3 * N/(N-1)/(N-2);
	m4 = M4 * (N-1)*(N*N-3*N+3)/N/N/N/N + M2*M2 * 3*(2*N-3)*(N-1)/N/N/N/N/N;
	mean = xbar;
	variance = m2;
	skewness = m3 / pow(m2,1.5);
	kurtosis = m4 / m2/m2;
}

void int_lif::utils::statsCalcTwoPass(const std::vector<double>& vec, double& mean, double& variance, double& skewness, double& kurtosis) {
	int N = vec.size();
	mean = 0.;
	for(auto& ii: vec) mean += ii;
	mean /= N;
	double m2 = 0., m3 = 0., m4 = 0.;
	for(auto& ii:vec) {
		double delta = ii - mean;
		m2 += delta*delta;
		m3 += delta*delta*delta;
		m4 += delta*delta*delta*delta;
	}
	m2 /= N;
	m3 /= N;
	m4 /= N;
	// m4 = m4 * (N-1)*(N*N-3*N+3)/N/N/N/N + m2*m2 * 3*(2*N-3)*(N-1)/N/N/N/N/N;
	// m3 = m3 * N/(N-1)/(N-2);
	// m2 = m2 / (N-1);
	variance = m2;
	skewness = m3 / pow(m2,1.5);
	kurtosis = m4 / m2/m2;
}

double int_lif::utils::JarqueBera(const std::vector<double>& vec) {
	double mean, variance, skewness, kurtosis;
	statsCalcTwoPass(vec, mean, variance, skewness, kurtosis);
	int N = vec.size();
	// std::cout << "mean = " << mean << " variance = " << variance << " skewness = " << skewness << " kurtosis = " << kurtosis << std::endl;
	return N/6*(skewness*skewness + (kurtosis-3)*(kurtosis-3)/4);
}

double int_lif::utils::mse(const std::vector<double>& vec1, const std::vector<double>& vec2) {
	double res = 0.;
	for(auto it1 = vec1.begin(), it2 = vec2.begin(); it1 != vec1.end() && it2 != vec2.end(); it1++, it2++) res += pow((*it1)-(*it2), 2);
	return res / vec1.size();
}

double int_lif::utils::mae(const std::vector<double>& vec1, const std::vector<double>& vec2) {
	double res = 0.;
	for(auto it1 = vec1.begin(), it2 = vec2.begin(); it1 != vec1.end() && it2 != vec2.end(); it1++, it2++) res += abs((*it1)-(*it2));
	return res / vec1.size();
}