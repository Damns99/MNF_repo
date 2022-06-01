#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <filesystem>
namespace fs = std::filesystem;

#include "lattice.h"
#include "bootstrap.h"
#include "text_io.h"
#include "cmdline_parser.h"

double mean(const std::vector<double>& vec) {
	double s = 0.;
	int length = vec.size();
	for(auto& ii: vec) s += ii / length;
	return s;
}

double twoPointConnectedCorrelation(int tau, double* yy, int len) {
	double res = 0., sum = 0.;
	for(int i = 0; i < len; i++) {
		res += yy[i] * yy[index1d(i + tau, len)] / len;
		sum += yy[i] / len;
	}
	res -= sum * sum;
	return res;
}

double twoPointConnectedCorrelation2(int tau, double* yy, int len) {
	double res = 0., sum = 0.;
	for(int i = 0; i < len; i++) {
		res += yy[i] * yy[i] * yy[index1d(i + tau, len)] * yy[index1d(i + tau, len)] / len;
		sum += yy[i] * yy[i] / len;
	}
	res -= sum * sum;
	return res;
}

int main(int argc, char* argv[]) {
	std::string folder, outfilename, conffilename;
	int nmeas, append;
	long seed = -434567;
	
	cmdlineParser::CmdlineParser parser;
    parser.addPosParameter<std::string>("folder", &folder, "0220520212823", "[std::string] Working folder name.");
    parser.addPosParameter<std::string>("conffile", &conffilename, "double_hole_conf", "[std::string] Input file name for the configurations.");
    parser.addPosParameter<int>("nmeas", &nmeas, 0, "[int] Number of configurations to analyze.");
    parser.addOptParameter<std::string>("outfile", &outfilename, "two_point_results", "[std::string] Output file name for the results.");
    parser.addOptParameter<int>("append", &append, 0, "[int] If != 0 append means to measfile instead of overwrtiting them.");
	
	if (parser.parseAll(argc, argv) == HELP_RETURN) return 0;
    parser.kickOff(argv[0]);
	
	fs::current_path(fs::current_path() / "measures" / folder);
	
	std::ofstream outfile1;
	outfile1.open(outfilename + "_tmp.txt", std::fstream::out);
	outfile1 << "#part \ttau \ttpccorr \ttpccorr2" << std::endl;
	for(int ii = 0; ii < nmeas; ii++) {
		load(conffilename + "-" + std::to_string(ii) + ".txt");
		for(int pp = 0; pp < nparticles; pp++) for(int jj = 0; jj < p_length / 2; jj++) {
			double conncorr = twoPointConnectedCorrelation(jj, y + pp * p_length, p_length);
			double conncorr2 = twoPointConnectedCorrelation2(jj, y + pp * p_length, p_length);
			outfile1 << pp << '\t' << jj << '\t' << conncorr << '\t' << conncorr2 << std::endl;
		}
	}
	outfile1.close();
	
	std::vector<double> particle, tau, tpccorr, tpccorr2;
	int ntemp = textIo::textIn(outfilename + "_tmp.txt", '\t', '#', &particle, &tau, &tpccorr, &tpccorr2);
	std::vector<double> ptcpcorr[nparticles][p_length / 2], ptcpcorr2[nparticles][p_length / 2];
	for(int ii = 0; ii < ntemp; ii++) {
		ptcpcorr[int(particle[ii])][int(tau[ii])].push_back(tpccorr[ii]);
		ptcpcorr2[int(particle[ii])][int(tau[ii])].push_back(tpccorr2[ii]);
	}
	
	std::ofstream outfile2;
	if (append) outfile2.open(outfilename + ".txt", std::fstream::out | std::fstream::app);
	else {
		outfile2.open(outfilename + ".txt", std::fstream::out);
		outfile2 << "#part \ttau \ttpccorr \tdtpccorr \ttpccorr2 \tdtpccorr2" << std::endl;
	}
	for(int pp = 0; pp < nparticles; pp++) for(int jj = 0; jj < p_length / 2; jj++) {
		double c, dc, c2, dc2;
		c = mean(ptcpcorr[pp][jj]);
		dc = bootstrap::bootstrapError(mean, ptcpcorr[pp][jj], 256, 64, seed *nparticles + jj);
		c2 = mean(ptcpcorr2[pp][jj]);
		dc2 = bootstrap::bootstrapError(mean, ptcpcorr2[pp][jj], 256, 64, seed *nparticles + jj);
		outfile2 << pp << '\t' << jj << '\t' << c << '\t' << dc << '\t' << c2 << '\t' << dc2 << std::endl;
	}
	outfile2.close();
}