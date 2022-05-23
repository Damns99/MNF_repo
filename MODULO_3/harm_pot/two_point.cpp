#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <filesystem>
namespace fs = std::filesystem;

#include "cmdline_parser.h"
#include "lattice.h"

double twoPointConnectedCorrelation(int tau, double* yy, int len) {
	double res = 0., sum = 0.;
	for(int i = 0; i < len; i++) {
		res += yy[i] * yy[index1d(i + tau, len)] / len;
		sum += yy[i] / len;
	}
	res -= sum * sum;
	return res;
}

int main(int argc, char* argv[]) {
	std::string folder, outfilename, conffilename;
	int nmeas, append;
	
	cmdlineParser::CmdlineParser parser;
    parser.addPosParameter<std::string>("folder", &folder, "0220520212823", "[std::string] Working folder name.");
    parser.addPosParameter<std::string>("conffile", &conffilename, "harm_pot_conf", "[std::string] Input file name for the configurations.");
    parser.addPosParameter<int>("nmeas", &nmeas, 0, "[int] Number of configurations to analyze.");
    parser.addOptParameter<std::string>("outfile", &outfilename, "two_point_results.txt", "[std::string] Output file for the results.");
    parser.addOptParameter<int>("append", &append, 0, "[int] If != 0 append means to measfile instead of overwrtiting them.");
	
	if (parser.parseAll(argc, argv) == HELP_RETURN) return 0;
    parser.kickOff(argv[0]);
	
	fs::current_path(fs::current_path() / "measures" / folder);
	
	std::ofstream outfile;
	if (append) outfile.open(outfilename, std::fstream::out | std::fstream::app);
	else {
		outfile.open(outfilename, std::fstream::out);
		outfile << "#part \ttau \ttpccorr" << std::endl;
	}
	for(int ii = 0; ii < nmeas; ii++) {
		load(conffilename + "-" + std::to_string(ii) + ".txt");
		for(int pp = 0; pp < nparticles; pp++) for(int jj = 0; jj < p_length / 2; jj++) {
			double conncorr = twoPointConnectedCorrelation(jj, y + pp * p_length, p_length);
			outfile << pp << '\t' << jj << '\t' << conncorr << std::endl;
		}
	}
	outfile.close();
}