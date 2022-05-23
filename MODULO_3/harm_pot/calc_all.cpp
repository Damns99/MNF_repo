#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <filesystem>
namespace fs = std::filesystem;

#include "bootstrap.h"
#include "text_io.h"
#include "cmdline_parser.h"

double mean(const std::vector<double>& vec) {
	double s = 0.;
	int length = vec.size();
	for(auto& ii: vec) s += ii / length;
	return s;
}

double derivative(const std::vector<double>& vec) {
	double s1 = 0., s2 = 0.;
	int length = vec.size();
	for(auto& ii: vec) {
		s1 += ii / length;
		s2 += (ii * ii) / length;
	}
	return s2 - s1 * s1;
}

int main(int argc, char* argv[]) {
	long seed = -7745;
    std::string folder, outfilename, measfilename;
    std::vector<double> particle, acc;
    std::vector<double> obs1, obs2;
	int npart;
	double x;
	int append;
	
	cmdlineParser::CmdlineParser parser;
    parser.addPosParameter<std::string>("folder", &folder, "0220520212823", "[std::string] Working folder name.");
    parser.addPosParameter<std::string>("measfile", &measfilename, "harm_pot_meas.txt", "[std::string] Input file for the measures.");
    parser.addPosParameter<double>("x", &x, 0, "[double] x coordinate for whatever parameter you are varying.");
    parser.addOptParameter<std::string>("outfile", &outfilename, "results.txt", "[std::string] Output file for the results.");
    parser.addOptParameter<int>("append", &append, 0, "[int] If != 0 append means to measfile instead of overwrtiting them.");
	if (parser.parseAll(argc, argv) == HELP_RETURN) return 0;
    parser.kickOff(argv[0]);
	
    fs::current_path(fs::current_path() / "measures" / folder);
    int nmeas = textIo::textIn(measfilename, '\t', '#', &particle, &obs1, &obs2, &acc);
	
	std::string tmp;
	std::ifstream tmpif;
	tmpif.open(measfilename, std::fstream::in);
	tmpif >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> npart;
	tmpif.close();
	int pnmeas = nmeas / npart;
	
	std::vector<double> pobs1[npart], pobs2[npart];
	for(int ii = 0; ii < npart; ii++) {
		pobs1[ii].reserve(pnmeas);
		pobs2[ii].reserve(pnmeas);
	}
	for(int ii = 0; ii < nmeas; ii++) {
		pobs1[int(particle[ii])].push_back(obs1[ii]);
		pobs2[int(particle[ii])].push_back(obs2[ii]);
	}
	
	std::ofstream outfile;
	if (append) outfile.open(outfilename, std::fstream::out | std::fstream::app);
	else {
		outfile.open(outfilename, std::fstream::out);
		outfile << "#x \tpart \t<obs1> \terror \t<obs2> \terror" << std::endl;
	}
	for(int ii = 0; ii < npart; ii++) {
		double o1m, do1m, o2m, do2m;
		o1m = mean(pobs1[ii]);
		do1m = bootstrap::bootstrapError(mean, pobs1[ii], 256, 64, seed * 1 + ii);
		o2m = mean(pobs2[ii]);
		do2m = bootstrap::bootstrapError(mean, pobs2[ii], 256, 64, seed * 2 + ii);
		
		outfile << x << '\t' << ii << '\t';
		outfile << o1m << " \t" << do1m << " \t" << o2m << " \t" << do2m << '\t';
	}
	outfile.close();
}