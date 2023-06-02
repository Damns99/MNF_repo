#include "integrate_LIF.h"

// single neuron Leaky Integrate and Fire model with external input
// {
//    Cm * dV/dt = -g * (V - Vrest) + I(t)
//    if V >= Vth => V -> V_reset and I(t:t+tref) = 0
// }

// empty overhead-timing method
// fill with constant
std::vector<double> int_lif::empty(double V0, double h, int N, std::vector<double>& I, double params[6]) {
	double Cm = params[0], g = params[1], Vrest = params[2], Vth = params[3], Vreset = params[4], tauref = params[5];
	double a = g/Cm*h;
	std::vector<double> V(N);
	int refr = 0;
	V[0] = V0;
	double new_V = V0;
	for(int n = 1; n < N; n++) {
		if(refr <= 0) new_V = V0;
		else {
			new_V = V0;
			refr--;
		}
		if(new_V >= Vth) {
			V[n] = Vreset;
			new_V = Vreset;
			refr = int(tauref / h);
		}
		else V[n] = new_V;
	}
	return V;
}

// 2 memory access method
// 
std::vector<double> int_lif::empty2(double V0, double h, int N, std::vector<double>& I, double params[6]) {
	double Cm = params[0], g = params[1], Vrest = params[2], Vth = params[3], Vreset = params[4], tauref = params[5];
	double a = g/Cm*h;
	std::vector<double> V(N);
	int refr = 0;
	V[0] = V0;
	double new_V = V0;
	for(int n = 1; n < N; n++) {
		if(refr <= 0) new_V = (1-a) * new_V + h/Cm * I[(3*n-1)%N];
		else {
			new_V = a * new_V;
			refr--;
		}
		if(new_V >= Vth) {
			V[n] = Vreset;
			new_V = Vreset;
			refr = int(tauref / h);
		}
		else V[n] = new_V;
	}
	return V;
}

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
			new_V = Vreset;
			refr = int(tauref / h);
		}
		else V[n] = new_V;
	}
	return V;
}

// for now no tauref must be 0
double int_lif::fwdEulerLocError(std::vector<double>& V_ref, double h_ref, double h, int Npoints, std::vector<double>& I_ref, double params[6]) {
	double Cm = params[0], g = params[1], Vrest = params[2], Vth = params[3], Vreset = params[4], tauref = params[5];
	double a = g/Cm*h;
	int N_ref = V_ref.size();
	int methodSteps = 1;
	int init_ii = ceil(methodSteps * h/h_ref);
	std::vector<double> locerrorvec;
	for(int ii = init_ii; ii < N_ref; ii += int((N_ref-init_ii)/Npoints)) {
		double new_V = (1-a) * V_ref[int(ii-h/h_ref)] + a * Vrest + h/Cm * I_ref[int(ii-h/h_ref)];
		if(new_V >= Vth) new_V = Vreset;
		locerrorvec.push_back(abs(new_V - V_ref[ii]));
	}
	return utils::mean(locerrorvec);
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
			new_V = Vreset;
			refr = int(tauref / h);
		}
		else V[n] = new_V;
	}
	return V;
}

// for now no tauref must be 0
double int_lif::bwdEulerLocError(std::vector<double>& V_ref, double h_ref, double h, int Npoints, std::vector<double>& I_ref, double params[6]) {
	double Cm = params[0], g = params[1], Vrest = params[2], Vth = params[3], Vreset = params[4], tauref = params[5];
	double a = g/Cm*h;
	int N_ref = V_ref.size();
	int methodSteps = 1;
	int init_ii = ceil(methodSteps * h/h_ref);
	std::vector<double> locerrorvec;
	for(int ii = init_ii; ii < N_ref; ii += int((N_ref-init_ii)/Npoints)) {
		double new_V = (V_ref[int(ii-h/h_ref)] + a * Vrest + h/Cm * I_ref[int(ii)]) / (1+a);
		if(new_V >= Vth) new_V = Vreset;
		locerrorvec.push_back(abs(new_V - V_ref[ii]));
	}
	return utils::mean(locerrorvec);
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
		/*double k1 = -g/Cm * (new_V - Vrest) + 1./Cm * I[n-1];
		double k2 = -g/Cm * (new_V + h * k1 - Vrest) + 1./Cm * I[n];
		new_V = new_V + h/2. * (k1 + k2);*/
		if(new_V >= Vth) {
			V[n] = Vreset;
			new_V = Vreset;
			refr = int(tauref / h);
		}
		else V[n] = new_V;
	}
	return V;
}

// for now no tauref must be 0
double int_lif::HeunLocError(std::vector<double>& V_ref, double h_ref, double h, int Npoints, std::vector<double>& I_ref, double params[6]) {
	double Cm = params[0], g = params[1], Vrest = params[2], Vth = params[3], Vreset = params[4], tauref = params[5];
	double a = g/Cm*h;
	int N_ref = V_ref.size();
	int methodSteps = 1;
	int init_ii = ceil(methodSteps * h/h_ref);
	std::vector<double> locerrorvec;
	for(int ii = init_ii; ii < N_ref; ii += int((N_ref-init_ii)/Npoints)) {
		double new_V = (a*a/2.-a+1) * V_ref[int(ii-h/h_ref)] - (a*a/2.-a) * Vrest + (1-a)/2. * h/Cm * I_ref[int(ii-h/h_ref)] + 1/2. * h/Cm * I_ref[int(ii)];
		/* double k1 = -g/Cm * (V_ref[int(ii-h/h_ref)] - Vrest) + 1./Cm * I_ref[int(ii-h/h_ref)];
		double k2 = -g/Cm * (V_ref[int(ii-h/h_ref)] + h * k1 - Vrest) + 1./Cm * I_ref[int(ii)];
		double new_V = V_ref[int(ii-h/h_ref)] + h/2. * (k1 + k2);*/
		if(new_V >= Vth) new_V = Vreset;
		locerrorvec.push_back(abs(new_V - V_ref[ii]));
	}
	return utils::mean(locerrorvec);
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
			new_V = Vreset;
			refr = int(tauref / h);
		}
		else V[n] = new_V;
	}
	return V;
}

// for now no tauref must be 0
double int_lif::RK4LocError(std::vector<double>& V_ref, double h_ref, double h, int Npoints, std::vector<double>& I_ref, double params[6]) {
	double Cm = params[0], g = params[1], Vrest = params[2], Vth = params[3], Vreset = params[4], tauref = params[5];
	double a = g/Cm*h;
	double c1 = (a*a*a*a/24.-a*a*a/6.+a*a/2.-a+1), c2 = 1-c1, c3 = (-a*a*a/24.+a*a/12.-a/6.+1/6.) * h/Cm, c4 = (a*a/12.-a/3.+2/3.) * h/Cm, c5 = 1/6. * h/Cm;
	int N_ref = V_ref.size();
	int methodSteps = 1;
	int init_ii = ceil(methodSteps * h/h_ref);
	std::vector<double> locerrorvec;
	for(int ii = init_ii; ii < N_ref; ii += int((N_ref-init_ii)/Npoints)) {
		double new_V = c1 * V_ref[int(ii-h/h_ref)] + c2 * Vrest + c3 * I_ref[int(ii-h/h_ref)] + c4 * I_ref[int(ii-h/h_ref/2)] + c5 * I_ref[int(ii)];
		if(new_V >= Vth) new_V = Vreset;
		locerrorvec.push_back(abs(new_V - V_ref[ii]));
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