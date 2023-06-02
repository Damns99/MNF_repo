#include "integrate_LIF.h"

// single neuron Leaky Integrate and Fire model with external input
// ADIMENSIONAL
// {
//    dy(x)/dx = -(y(x) - 1) + z(x)
//    if y(x) >= yth => y(x) -> yreset
// }

// Forward Euler method
// dy(n)/dx = (y(n+1) - y(n)) / h
std::vector<double> int_lif::fwdEuler(double y0, double h, int N, std::vector<double>& z, double params[2], std::vector<double>* spiketimes) {
	double yth = params[0], yreset = params[1];
	std::vector<double> y(N);
	y[0] = y0;
	double newy = y0;
	for(int n = 1; n < N; n++) {
		newy = (1-h) * newy + h * (z[n-1] + 1);
		if(newy <= yth) {
			y[n] = yreset;
			newy = yreset;
			if(spiketimes != nullptr) spiketimes->push_back(n*h);
		}
		else y[n] = newy;
	}
	return y;
}

double int_lif::fwdEulerLocError(std::vector<double>& y_ref, double h_ref, double h, int Npoints, std::vector<double>& z_ref, double params[2]) {
	double yth = params[0], yreset = params[1];
	int N_ref = y_ref.size();
	int methodSteps = 1;
	int init_ii = ceil(methodSteps * h/h_ref);
	std::vector<double> locerrorvec;
	for(int ii = init_ii; ii < N_ref; ii += int((N_ref-init_ii)/Npoints)) {
		double newy = (1-h) * y_ref[int(ii-h/h_ref)] + h * (z_ref[int(ii-h/h_ref)] + 1);
		if(newy <= yth) newy = yreset;
		locerrorvec.push_back(abs(newy - y_ref[ii]));
	}
	return utils::mean(locerrorvec);
}

// Backward Euler method
// dy(n)/dx = (y(n) - y(n-1)) / h
std::vector<double> int_lif::bwdEuler(double y0, double h, int N, std::vector<double>& z, double params[2], std::vector<double>* spiketimes) {
	double yth = params[0], yreset = params[1];
	std::vector<double> y(N);
	y[0] = y0;
	double newy = y0;
	for(int n = 1; n < N; n++) {
		newy = (newy + h * (z[n] + 1)) / (1+h);
		if(newy <= yth) {
			y[n] = yreset;
			newy = yreset;
			if(spiketimes != nullptr) spiketimes->push_back(n*h);
		}
		else y[n] = newy;
	}
	return y;
}

double int_lif::bwdEulerLocError(std::vector<double>& y_ref, double h_ref, double h, int Npoints, std::vector<double>& z_ref, double params[2]) {
	double yth = params[0], yreset = params[1];
	int N_ref = y_ref.size();
	int methodSteps = 1;
	int init_ii = ceil(methodSteps * h/h_ref);
	std::vector<double> locerrorvec;
	for(int ii = init_ii; ii < N_ref; ii += int((N_ref-init_ii)/Npoints)) {
		double newy = (y_ref[int(ii-h/h_ref)] + h * (z_ref[int(ii)] + 1)) / (1+h);
		if(newy <= yth) newy = yreset;
		locerrorvec.push_back(abs(newy - y_ref[ii]));
	}
	return utils::mean(locerrorvec);
}

// Heun method
// 2 step RK with equal derivative weights and no intermediate steps
std::vector<double> int_lif::Heun(double y0, double h, int N, std::vector<double>& z, double params[2], std::vector<double>* spiketimes) {
	double yth = params[0], yreset = params[1];
	double a = h - h*h/2;
	std::vector<double> y(N);
	y[0] = y0;
	double newy = y0;
	for(int n = 1; n < N; n++) {
		newy = (1-a) * newy + (a-h/2) * z[n-1] + (h/2) * z[n] + a;
		if(newy <= yth) {
			y[n] = yreset;
			newy = yreset;
			if(spiketimes != nullptr) spiketimes->push_back(n*h);
		}
		else y[n] = newy;
	}
	return y;
}

double int_lif::HeunLocError(std::vector<double>& y_ref, double h_ref, double h, int Npoints, std::vector<double>& z_ref, double params[2]) {
	double yth = params[0], yreset = params[1];
	double a = h - h*h/2;
	int N_ref = y_ref.size();
	int methodSteps = 2;
	int init_ii = ceil(methodSteps * h/h_ref);
	std::vector<double> locerrorvec;
	for(int ii = init_ii; ii < N_ref; ii += int((N_ref-init_ii)/Npoints)) {
		double newy = (1-a) * y_ref[int(ii-h/h_ref)] + (a-h/2) * z_ref[int(ii-h/h_ref)] + (h/2) * z_ref[int(ii)] + a;
		if(newy <= yth) newy = yreset;
		locerrorvec.push_back(abs(newy - y_ref[ii]));
	}
	return utils::mean(locerrorvec);
}

// RungeKutta4 method
// I has 2N elements, one every h/2 seconds
std::vector<double> int_lif::RK4(double y0, double h, int N, std::vector<double>& z, double params[2], std::vector<double>* spiketimes) {
	double yth = params[0], yreset = params[1];
	double a1 = -h*h*h*h/24. + h*h*h/6. - h*h/2. + h, a2 = h*h*h/12. - h*h/3. + 2.*h/3.;
	std::vector<double> y(N);
	y[0] = y0;
	double newy = y0;
	for(int n = 1; n < N; n++) {
		newy = (1-a1) * newy + (a1-a2) * z[2*n-2] + (a2-h/6) * z[2*n-1] + (h/6) * z[2*n] + a1;
		if(newy <= yth) {
			y[n] = yreset;
			newy = yreset;
			if(spiketimes != nullptr) spiketimes->push_back(n*h);
		}
		else y[n] = newy;
	}
	return y;
}

double int_lif::RK4LocError(std::vector<double>& y_ref, double h_ref, double h, int Npoints, std::vector<double>& z_ref, double params[2]) {
	double yth = params[0], yreset = params[1];
	double a1 = -h*h*h*h/24. + h*h*h/6. - h*h/2. + h, a2 = h*h*h/12. - h*h/3. + 2.*h/3.;
	int N_ref = y_ref.size();
	int methodSteps = 4;
	int init_ii = ceil(methodSteps * h/h_ref);
	std::vector<double> locerrorvec;
	for(int ii = init_ii; ii < N_ref; ii += int((N_ref-init_ii)/Npoints)) {
		double newy = (1-a1) * y_ref[int(ii-h/h_ref)] + (a1-a2) * z_ref[int(ii-h/h_ref)] + (a2-h/6) * z_ref[int(ii-h/h_ref/2)] + (h/6) * z_ref[int(ii)] + a1;
		if(newy <= yth) newy = yreset;
		locerrorvec.push_back(abs(newy - y_ref[ii]));
	}
	return utils::mean(locerrorvec);
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

std::vector<double> int_lif::currents::sine_wave(int n, int period, double amplitude, int phase, double offset) {
	std::vector<double> res(n);
	for(int ii = 0; ii < n; ii++) res[ii] = offset + amplitude * sin(2.*M_PI*(ii+phase)/period);
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