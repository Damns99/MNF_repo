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

// Utils

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

// Rules

int updateFreeParticle(int part) {
	double delta = 1.;
	int x = gen.randL(0, plength) + part * plength;
	int yp = delta * gen.randF() + y[x];
	double tmp0 = (yp - y[x]), tmp1 = (yp + y[x]), tmp2 = (y[links[2 * x + 0]] + y[links[2 * x + 1]]);
	double ds = tmp0 * (tmp1 - tmp2);
	double r = exp(-ds);
	double rr = gen.randF();
	if (rr < r) {
		y[x] = yp;
		obs1[part] += (tmp0 * tmp1) / (1. * p_length);
		obs2[part] += (ds) / (1. * p_length);
		return 1;
	}
	return 0;
}

template<int T>
int updateFreeParticles() {
	return updateFreeParticle(T);
}

// Main

int main(int argc, char* argv[]) {
	gErrorIgnoreLevel = kWarning;
	
    std::string measfilename, outfilename, infilename, outfoldername;
    long seed;
    int nmeas, ncycles, nskip, length, geom, init_mode, append, do_snapshots, do_timing;
    double extrafield, beta;
	
	cmdlineParser::CmdlineParser parser;
    parser.addPosParameter<int>("init_mode", &init_mode, HOT, "[int] Init mode for the first lattice state. HOT = 1, COLD = 2, FROM_FILE = 3");
    parser.addPosParameter<std::string>("measfile", &measfilename, "metro_ising_meas.txt", "[std::string] Output file for the measures.");
    parser.addPosParameter<std::string>("outfile", &outfilename, "metro_ising_conf", "[std::string] Output file for the last iteration lattice state.");
	parser.addOptParameter<std::string>("infile", &infilename, "metro_ising_conf.txt", "[std::string] Intput file for the first iteration lattice state.");
	parser.addOptParameter<std::string>("folder", &outfoldername, "", "[std::string] Output folder name.");
    parser.addOptParameter<long>("seed", &seed, -std::time(NULL), "[long] Random number generator (Ran2) seed.");
    parser.addOptParameter<int>("nm", &nmeas, 1024, "[int] Number of measures to take.");
    parser.addOptParameter<int>("nc", &ncycles, 1, "[int] Number of cycles of iteration between each measure.");
    parser.addOptParameter<int>("ns", &nskip, 0, "[int] Number of measures to skip at start.");
    parser.addOptParameter<int>("len", &length, 10, "[int] Side length of the square ");
    parser.addOptParameter<double>("extf", &extrafield, 0., "[double] Extern magnetic field (in units of JACC).");
    parser.addOptParameter<double>("beta", &beta, 0.3, "[double] One over temperature of the system (in units of JACC).");
    parser.addOptParameter<int>("append", &append, 0, "[int] If != 0 append measures to measfile instead of overwrtiting them.");
    parser.addOptParameter<int>("geom", &geom, SQUARE_GEOM, "[int] Select lattice geomerty: SQUARE = 0, EXAGONAL = 1, TRIANGLUAR = 2");
    parser.addOptParameter<int>("snapshots", &do_snapshots, 0, "[int] If != 0 takes snapshots of the lattice at each measure.");
    parser.addOptParameter<int>("timing", &do_timing, 0, "[int] If != 0 prints time to take each update cycle and measure.");
	
    if (parser.parseAll(argc, argv) == HELP_RETURN) return 0;
	
    parser.kickOff(argv[0]);
	
	createLattice(length, geom, beta, extrafield, seed, init_mode, infilename);
	
	fs::current_path(fs::current_path() / "measures");
	std::string folder;
	if(outfoldername != "") folder = outfoldername;
	else folder = currentTimeDate();
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
		measfile << "# nm = " << nmeas << " nc = " << ncycles << " len = " << length << " extf = " << extrafield << " beta = " << beta << std::endl;
		measfile << "#energy \tmagnetization \tacceptance" << std::endl;
	}
	int percent = 0;
    std::chrono::time_point<std::chrono::steady_clock> start, end;
	
	for(int ii = 0; ii < nskip; ii++) {
		for(int jj = 0; jj < ncycles; jj++) for(int kk = 0; kk < length * length; kk ++) updateMetropolis();
	}
	for(int ii = 0; ii < nmeas; ii++) {
		printPercent(ii, percent, nmeas);
		int acc = 0;
        if(do_timing) start = std::chrono::steady_clock::now();
		for(int jj = 0; jj < ncycles; jj++) for(int kk = 0; kk < length * length; kk ++) acc += updateMetropolis();
        if(do_timing) end = std::chrono::steady_clock::now();
		measfile << energy << '\t' << magnetization << '\t' << 1. * acc / ncycles / length / length << '\n';
		if(do_snapshots) {
			std::string snapname = "snapshot" + makeFixedLength(ii + 1, lognmeas) + ".png";
			snapshot(snapname);
			save(outfilename + "-" + std::to_string(ii) + ".txt");
		}
        if(do_timing) std::cout << "updated + measured in: " << (end - start).count() / 1000000. << std::endl;
	}
	
	measfile.close();
	save(outfilename + ".txt");
	
	std::cout << std::endl;
	return 0;
}