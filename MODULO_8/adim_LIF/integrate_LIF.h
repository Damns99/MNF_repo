#ifndef INTEGRATE_LIF_H
#define INTEGRATE_LIF_H

#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <math.h>
#include "box_muller.h"

// single neuron Leaky Integrate and Fire model with external input
// {
//    Cm * dV/dt = -g * (V - Vrest) + I(t)
//    if V >= Vth => V -> V_reset and I(t:t+tref) = 0
// }

namespace int_lif {
	// Forward Euler method
	// dV(n)/dt = (V(n+1) - V(n)) / h
	std::vector<double> fwdEuler(double V0, double h, int N, std::vector<double>& I, double params[6], std::vector<double>* spiketimes = nullptr);
	double fwdEulerLocError(std::vector<double>& V_ref, double h_ref, double h, int Npoints, std::vector<double>& I_ref, double params[6]);

	// Backward Euler method
	// dV(n)/dt = (V(n) - V(n-1)) / h
	std::vector<double> bwdEuler(double V0, double h, int N, std::vector<double>& I, double params[6], std::vector<double>* spiketimes = nullptr);
	double bwdEulerLocError(std::vector<double>& V_ref, double h_ref, double h, int Npoints, std::vector<double>& I_ref, double params[6]);

	// Heun method
	// 
	std::vector<double> Heun(double V0, double h, int N, std::vector<double>& I, double params[6], std::vector<double>* spiketimes = nullptr);
	double HeunLocError(std::vector<double>& V_ref, double h_ref, double h, int Npoints, std::vector<double>& I_ref, double params[6]);

	// RungeKutta4 method
	// 
	std::vector<double> RK4(double V0, double h, int N, std::vector<double>& I, double params[6], std::vector<double>* spiketimes = nullptr);
	double RK4LocError(std::vector<double>& V_ref, double h_ref, double h, int Npoints, std::vector<double>& I_ref, double params[6]);

	// currents
	namespace currents {

		std::vector<double> square_wave(int n, int period, double amplitude, int phase, double offset);
		
		std::vector<double> sine_wave(int n, int period, double amplitude, int phase, double offset);
		std::vector<double> sine_wave(int n, double h, double period, double amplitude, double phase, double offset);

		std::vector<double> constant_pulse(int n, int start, int duration, double amplitude);

		std::vector<double> pulse_train(int n, int start, int duration, double amplitude, int interval, int numpulses);
		
		std::vector<double> white_noise(int n, double mean, double variance, long seed = -9542342);
		
		std::vector<double> ou_noise(int n, double mean, double sigma, double tau, double h, long seed = -9542342);

	}

	// utils
	namespace utils {

		std::vector<double> linspace(double x0, double x1, double n);

		std::vector<double> logspace(double e0, double e1, double n);

		double mean(const std::vector<double>& vec);

		double stdev(const std::vector<double>& vec);

		void statsCalcOnePass(const std::vector<double>& vec, double& mean, double& variance, double& skewness, double& kurtosis);

		void statsCalcTwoPass(const std::vector<double>& vec, double& mean, double& variance, double& skewness, double& kurtosis);

		double JarqueBera(const std::vector<double>& vec);

		void inline printPercent(int ii, int& percent, const int& N, std::string prefix = "") {
			if(100 * ii / N > percent) {
				percent = 100 * ii / N;
				std::cout << "\r" << prefix << percent << "%";
				std::flush(std::cout);
			}
		}
		
		double mse(const std::vector<double>& vec1, const std::vector<double>& vec2);
		
		double mae(const std::vector<double>& vec1, const std::vector<double>& vec2);

	}

}

#endif