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

#include <TError.h>

std::string makeFixedLength(const int i, const int l)
{
    std::stringstream s;
    s << std::setfill('0') << std::setw(l) << i;
    return s.str();
}

std::string currentTimeDate() {
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);

    std::stringstream ss;
    ss << std::put_time(std::localtime(&in_time_t), "%y%m%d%H%M%S");
    return ss.str();
}

void inline printPercent(int ii, int& percent, const int& nmeas) {
	if(100 * ii / nmeas > percent) {
		percent = 100 * ii / nmeas;
		std::cout << "\r" << percent << "%";
		std::flush(std::cout);
	}
}

int main(int argc, char* argv[]) {
	gErrorIgnoreLevel = kWarning;
	
    std::string measfilename, outfilename, infilename;
    long seed;
    int nmeas, ncycles, length, geom, init_mode, append, do_snapshots;
    double extrafield, beta;
	
	cmdlineParser::CmdlineParser parser;
    parser.addPosParameter<int>("init_mode", &init_mode, HOT, "[int] Init mode for the first lattice state. HOT = 1, COLD = 2, FROM_FILE = 3");
    parser.addPosParameter<std::string>("measfile", &measfilename, "metro_ising_meas.txt", "[std::string] Output file for the measures.");
    parser.addPosParameter<std::string>("outfile", &outfilename, "metro_ising_conf.txt", "[std::string] Output file for the last iteration lattice state.");
	parser.addOptParameter<std::string>("infile", &infilename, "metro_ising_conf.txt", "[std::string] Intput file for the first iteration lattice state.");
    parser.addOptParameter<long>("seed", &seed, -42, "[long] Random number generator (Ran2) seed.");
    parser.addOptParameter<int>("nm", &nmeas, 1024, "[int] Number of measures to take.");
    parser.addOptParameter<int>("nc", &ncycles, 1, "[int] Number of cycles of iteration between each measure.");
    parser.addOptParameter<int>("len", &length, 10, "[int] Side length of the square ");
    parser.addOptParameter<double>("extf", &extrafield, 0., "[double] Extern magnetic field (in units of JACC).");
    parser.addOptParameter<double>("beta", &beta, 0.3, "[double] One over temperature of the system (in units of JACC).");
    parser.addOptParameter<int>("append", &append, 0, "[int] If != 0 append measures to measfile instead of overwrtiting them.");
    parser.addOptParameter<int>("geom", &geom, SQUARE_GEOM, "[int] Select lattice geomerty: SQUARE = 0, EXAGONAL = 1, TRIANGLUAR = 2");
    parser.addOptParameter<int>("snapshots", &do_snapshots, 0, "[int] If != 0 takes snapshots of the lattice at each measure.");
	
    if (parser.parseAll(argc, argv) == HELP_RETURN) return 0;
	
    parser.kickOff(argv[0]);
	
	createLattice(length, geom, beta, extrafield, seed, init_mode, infilename);
	
	std::string folder = "measures" + currentTimeDate();
	fs::create_directory(folder);
	fs::current_path(fs::current_path() / folder);
	
	int lognmeas = log10(nmeas) + 1;
	if(do_snapshots) {
		std::string snapname = "snapshot" + makeFixedLength(0, lognmeas) + ".png";
		snapshot(snapname);
	}
	
	std::ofstream measfile;
	if (append) measfile.open(measfilename, std::fstream::out | std::fstream::app);
	else {
		measfile.open(measfilename, std::fstream::out);
		measfile << "#energy \tmagnetization \tacceptance" << std::endl;
	}
	int percent = 0;
	
	for(int ii = 0; ii < nmeas; ii++) {
		printPercent(ii, percent, nmeas);
		int acc = 0;
		for(int jj = 0; jj < ncycles; jj++) for(int kk = 0; kk < length * length; kk ++) acc += updateMetropolis();
		measfile << energy << '\t' << magnetization << '\t' << 1. * acc / ncycles << '\n';
		if(do_snapshots) {
			std::string snapname = "snapshot" + makeFixedLength(ii + 1, lognmeas) + ".png";
			snapshot(snapname);
		}
	}
	
	measfile.close();
	save(outfilename);
	
	std::cout << std::endl;
	return 0;
}