#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <math.h>
#include <filesystem>
namespace fs = std::filesystem;

#include "cmdline_parser.h"
#include "integrate_LIF.h"
#include "text_io.h"

int main(int argc, char* argv[]) {
	// Voltages in [Vrest = -70e-3 V], Currents in [g*Vrest = -7e-10 A], Times in [Cm/g = 1e-2 s]
	double yth, yreset;
	std::string method, current;
	double simtime;
	int N;
	double x0, y0;
	double amplitude, period, offset, phase, start, duration, currmean, currvariance;
	int numpulses;
	std::string filename;
	
	cmdlineParser::CmdlineParser parser;
    parser.addPosParameter<double>("yth", &yth, 0.7, "y threshold for spiking [double]");
    parser.addPosParameter<double>("yreset", &yreset, 1., "y reset value after spike [double]");
    parser.addPosParameter<std::string>("method", &method, "fwdEuler", "select method: fwdEuler, bwdEuler, Heun, RK4 [string]");
    parser.addPosParameter<double>("simtime", &simtime, 1.e1, "time to simulate beginning from x0 [double]");
    parser.addPosParameter<int>("N", &N, 1e3, "number of steps, (x0,y0) included [int]");
    parser.addPosParameter<double>("x0", &x0, 0., "x0 [double]");
    parser.addPosParameter<double>("y0", &y0, 0., "y0 [double]");
    parser.addPosParameter<std::string>("current", &current, "sine_wave", "select cuurent: sine_wave, pulse_train, white_noise]");
	parser.addOptParameter<double>("amplitude", &amplitude, 0., "current amplitude (used iff current = sine_wave, pulse_train) [double]");
	parser.addOptParameter<double>("period", &period, 0., "current period (used iff current = sine_wave, pulse_train) [double]");
	parser.addOptParameter<double>("offset", &offset, 0., "current offset (used iff current = sine_wave) [double]");
	parser.addOptParameter<double>("phase", &phase, 0., "current phase (used iff current = sine_wave) [double]");
	parser.addOptParameter<double>("start", &start, 0., "current pulses' start (used iff current = pulse_train) [double]");
	parser.addOptParameter<double>("duration", &duration, 0., "current pulses' duration (used iff current = pulse_train) [double]");
	parser.addOptParameter<int>("numpulses", &numpulses, 0., "current number of pulses (used iff current = pulse_train) [double]");
	parser.addOptParameter<double>("currmean", &currmean, 0., "current mean (used iff current = white_noise) [double]");
	parser.addOptParameter<double>("currvariance", &currvariance, 0., "current variance (used iff current = white_noise) [double]");
    parser.addPosParameter<std::string>("filename", &filename, "out", "output file name without extension [string]");
	
	if(parser.parseAll(argc, argv) == HELP_RETURN) return 0;
    parser.kickOff(argv[0]);
	
	double params[2]; params[0] = yth; params[1] = yreset;
	double h = simtime/(N-1);
	double h1 = simtime/((N+1)/2-1);
	
	std::vector<double> x, y, z, spkt;
	x.push_back(x0);
	y.push_back(y0);
	
	std::string currparstring;
	if(current == "pulse_train") {
		currparstring = "start: "+std::to_string(start)+" duration: "+std::to_string(duration)+" amplitude: "+std::to_string(amplitude)+" period: "+std::to_string(period)+" numpulses: "+std::to_string(numpulses);
	}
	else if(current == "white_noise") {
		currparstring = "currmean: "+std::to_string(currmean)+" currvariance: "+std::to_string(currvariance);
	}
	else if(current == "sine_wave") {
		currparstring = "period: "+std::to_string(period)+" amplitude: "+std::to_string(amplitude)+" phase: "+std::to_string(phase)+" offset: "+std::to_string(offset);
	}
	else {
		std::cout << "Error: invalid current" << std::endl;
		return 1;
	}
	
	fs::current_path(fs::current_path() / "measures");
	textIo::textOut(filename+"_x_y.txt", '\t', '#', "N: "+std::to_string(N)+"    h: "+std::to_string(h)+"    method: "+method+"    current: "+current+" "+currparstring, 0, false, x);	
	textIo::textOut(filename+"_spkt.txt", '\t', '#', "N: "+std::to_string(N)+"    h: "+std::to_string(h)+"    method: "+method+"    current: "+current+" "+currparstring, 0, false, spkt);
	textIo::textOut(filename+"_x_y.txt", '\t', '#', "x\ty", 0, true, x);
	
	while(N > 0) {
		int NN = N < 1e3-1 ? N : 1e3-1;
		x0 = x[x.size()-1];
		y0 = y[y.size()-1];
		
		std::vector<double> z;
		if(current == "pulse_train") {
			z = int_lif::currents::pulse_train(NN, int(start/h), int(duration/h), amplitude, int((period-duration)/h), numpulses);
		}
		else if(current == "white_noise") {
			z = int_lif::currents::white_noise(NN, currmean, currvariance);
		}
		else if(current == "sine_wave") {
			z = int_lif::currents::sine_wave(NN, h, period, amplitude, phase+x0, offset);
		}
		
		if(method == "fwdEuler") {
			x = int_lif::utils::linspace(x0, x0 + (NN-1)*h, NN);
			y = int_lif::fwdEuler(y0, h, NN, z, params, &spkt, x0);
		}
		else if(method == "bwdEuler") {
			x = int_lif::utils::linspace(x0, x0 + (NN-1)*h, NN);
			y = int_lif::bwdEuler(y0, h, NN, z, params, &spkt, x0);
		}
		else if(method == "Heun") {
			x = int_lif::utils::linspace(x0, x0 + (NN-1)*h, NN);
			y = int_lif::Heun(y0, h, NN, z, params, &spkt, x0);
		}
		else if(method == "RK4") {
			if(N%2 == 0) {
				std::cout << "Error: N must be odd for RK4 method to work" << std::endl;
				return 1;
			}
			int N1 = (NN+1)/2;
			x = int_lif::utils::linspace(x0, x0 + (NN-1)*h, N1);
			y = int_lif::RK4(y0, h1, N1, z, params, &spkt, x0);
		}
		else if(method == "intTrapezioidal") {
			x = int_lif::utils::linspace(x0, x0 + (NN-1)*h, NN);
			y = int_lif::intTrapezioidal(y0, h, NN, z, params, &spkt, x0);
		}
		else if(method == "intSimpson") {
			if(N%2 == 0) {
				std::cout << "Error: N must be odd for intSimpson method to work" << N << std::endl;
				return 1;
			}
			int N1 = (NN+1)/2;
			x = int_lif::utils::linspace(x0, x0 + (NN-1)*h, N1);
			y = int_lif::intSimpson(y0, h1, N1, z, params, &spkt, x0);
		}
		else {
			std::cout << "Error: invalid method" << std::endl;
			return 1;
		}
	
		textIo::textOut(filename+"_x_y.txt", '\t', '#', "", x.size()-1, true, x, y);
		
		if(NN > 1) N -= NN-1;
		else N = 0;
	}
	
	std::vector<double> x_tmp(1,x[x.size()-1]), y_tmp(1,y[y.size()-1]);
	textIo::textOut(filename+"_x_y.txt", '\t', '#', "", 1, true, x_tmp, y_tmp);
	textIo::textOut(filename+"_spkt.txt", '\t', '#', "spkt", spkt.size(), true, spkt);
	
}