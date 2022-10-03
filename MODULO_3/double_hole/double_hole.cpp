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

int currentIndex = -1;
int nextIndex() {
	currentIndex++;
	if(currentIndex >= p_length) currentIndex = 0;
	return currentIndex;
}

double lambda = 1.;
int antiperiodic = 0;
int updateDoubleHoles(int part) {
	double eta = beta / p_length, alpha = lambda;
	double delta = 2. * sqrt(eta);
	int x = gen.randL(0, p_length) + part * p_length;
	// int x = nextIndex() + part * p_length;
	double yp = delta * (2. * gen.randF() - 1.) + y[x];
	double tmp0 = (yp - y[x]), tmp1 = (yp + y[x]), tmp2 = (y[links[2 * x + 0]] + y[links[2 * x + 1]]);
	if(antiperiodic && x == 0) tmp2 = (- y[links[2 * x + 0]] + y[links[2 * x + 1]]);
	else if(antiperiodic && x == p_length - 1) tmp2 = (y[links[2 * x + 0]] - y[links[2 * x + 1]]);
	double tmp3 = (yp * yp + y[x] * y[x] - 2.);
	double ds = tmp0 * (tmp1 - tmp2) / eta + tmp0 * tmp1 * tmp3 * eta * alpha;
	double r = exp(-ds);
	double rr = gen.randF();
	if (rr < r) {
		y[x] = yp;
		obs1[part] += (tmp0 * tmp1) / p_length;
		obs2[part] += 2. * (tmp0 * (tmp1 - tmp2)) / p_length;
		return 1;
	}
	return 0;
}

template<int T>
int updateDoubleHoles() {
	return updateDoubleHoles(T);
}

// Main

int main(int argc, char* argv[]) {
	gErrorIgnoreLevel = kWarning;
	
    std::string measfilename, outfilename, infilename, outfoldername;
    long seed;
    int nmeas, ncycles, nskip, length, nparticles, init_mode, append, do_snapshots, do_timing;
    double beta;
	
	cmdlineParser::CmdlineParser parser;
    parser.addPosParameter<int>("init_mode", &init_mode, HOT, "[int] Init mode for the first lattice state. HOT = 1, COLD = 2, FROM_FILE = 3");
    parser.addPosParameter<std::string>("measfile", &measfilename, "double_hole_meas.txt", "[std::string] Output file for the measures.");
    parser.addPosParameter<std::string>("outfile", &outfilename, "double_hole_conf", "[std::string] Output file for the last iteration lattice state.");
	parser.addOptParameter<std::string>("infile", &infilename, "double_hole_conf.txt", "[std::string] Intput file for the first iteration lattice state.");
	parser.addOptParameter<std::string>("folder", &outfoldername, "", "[std::string] Output folder name.");
    parser.addOptParameter<long>("seed", &seed, -std::time(NULL), "[long] Random number generator (Ran2) seed.");
    parser.addOptParameter<int>("nm", &nmeas, 1024, "[int] Number of measures to take.");
    parser.addOptParameter<int>("nc", &ncycles, 1, "[int] Number of cycles of iteration between each measure.");
    parser.addOptParameter<int>("ns", &nskip, 0, "[int] Number of measures to skip at start.");
    parser.addOptParameter<int>("len", &length, 10, "[int] Total length of the lattice.");
    parser.addOptParameter<int>("npart", &nparticles, 1, "[int] Number of particles.");
    parser.addOptParameter<double>("beta", &beta, 3., "[double] One over temperature of the system.");
    parser.addOptParameter<double>("lambda", &lambda, 1., "[double] Height of the wall.");
    parser.addOptParameter<int>("append", &append, 0, "[int] If != 0 append measures to measfile instead of overwrtiting them.");
    parser.addOptParameter<int>("snapshots", &do_snapshots, 0, "[int] If != 0 takes snapshots of the lattice at each measure.");
    parser.addOptParameter<int>("timing", &do_timing, 0, "[int] If != 0 prints time to take each update cycle and measure.");
    parser.addOptParameter<int>("antiperiodic", &antiperiodic, 0, "[int] If != 0 uses antiperiodic boundary conditions, periodic otherwise.");
	
    if (parser.parseAll(argc, argv) == HELP_RETURN) return 0;
	
    parser.kickOff(argv[0]);
	
	createLattice(length, nparticles, beta, seed, init_mode, infilename);
	
	fs::current_path(fs::current_path() / "measures");
	std::string folder;
	if(outfoldername != "") folder = outfoldername;
	else folder = currentTimeDate();
	fs::create_directory(folder);
	std::cout << "Created folder " << folder << std::endl;
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
		measfile << "# nm = " << nmeas << " nc = " << ncycles << " len = " << length << " npart = " << nparticles << " beta = " << beta << std::endl;
		measfile << "#part \tobs1 \tobs2 \tacc" << std::endl;
	}
	int percent = 0;
    std::chrono::time_point<std::chrono::steady_clock> start, end;
	
	if(nparticles > 0) addRule(updateDoubleHoles<0>, length / nparticles);
	if(nparticles > 1) addRule(updateDoubleHoles<1>, length / nparticles);
	if(nparticles > 2) addRule(updateDoubleHoles<2>, length / nparticles);
	if(nparticles > 3) addRule(updateDoubleHoles<3>, length / nparticles);
	if(nparticles > 4) addRule(updateDoubleHoles<4>, length / nparticles);
	if(nparticles > 5) addRule(updateDoubleHoles<5>, length / nparticles);
	if(nparticles > 6) addRule(updateDoubleHoles<6>, length / nparticles);
	if(nparticles > 7) addRule(updateDoubleHoles<7>, length / nparticles);
	
	for(int ii = 0; ii < nskip; ii++) {
		for(int jj = 0; jj < ncycles; jj++) update();
	}
	for(int ii = 0; ii < nmeas; ii++) {
		printPercent(ii, percent, nmeas);
        if(do_timing) start = std::chrono::steady_clock::now();
		double acc = 0.;
		for(int jj = 0; jj < ncycles; jj++) acc += update() / ncycles;
        if(do_timing) end = std::chrono::steady_clock::now();
		for(int pp = 0; pp < nparticles; pp++) measfile << pp << '\t' << obs1[pp] << '\t' << obs2[pp] << '\t' << acc << '\n';
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