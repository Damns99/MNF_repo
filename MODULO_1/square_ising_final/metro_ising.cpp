#include <iostream>
#include <vector>
#include <string>
#include <fstream>
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
	
    std::string measfilename, outfilename, infilename, outfoldername;
    long seed;
    int nmeas, ncycles, nskip, length, init_mode, append, do_snapshots;
    double extrafield, beta;
	
	cmdlineParser::CmdlineParser parser;
    parser.addPosParameter<std::string>("measfile", &measfilename, "metro_ising_meas", "[std::string] Output file for the measures (no file extension).");
    parser.addPosParameter<std::string>("outfile", &outfilename, "metro_ising_conf", "[std::string] Output file for the last iteration lattice state.");
	parser.addOptParameter<std::string>("infile", &infilename, "metro_ising_conf", "[std::string] Intput file for the first iteration lattice state.");
	parser.addOptParameter<std::string>("folder", &outfoldername, "", "[std::string] Output folder name.");
    parser.addOptParameter<int>("init", &init_mode, HOT, "[int] Init mode for the first lattice state. HOT = 1, COLD = 2, FROM_FILE = 3");
    parser.addOptParameter<long>("seed", &seed, -std::time(NULL), "[long] Random number generator (Ran2) seed.");
    parser.addOptParameter<int>("nm", &nmeas, 1024, "[int] Number of measures to take.");
    parser.addOptParameter<int>("nc", &ncycles, 1, "[int] Number of cycles of iteration between each measure.");
    parser.addOptParameter<int>("ns", &nskip, 0, "[int] Number of measures to skip at start.");
    parser.addOptParameter<int>("len", &length, 128, "[int] Side length of the square ");
    parser.addOptParameter<double>("extf", &extrafield, 0., "[double] Extern magnetic field (in units of JACC).");
    parser.addOptParameter<double>("beta", &beta, 0.3, "[double] One over temperature of the system (in units of JACC).");
    parser.addOptParameter<int>("append", &append, 0, "[int] If != 0 append measures to measfile instead of overwrtiting them.");
    parser.addOptParameter<int>("snapshots", &do_snapshots, 0, "[int] If != 0 takes snapshots of the lattice at each measure.");
	
    if (parser.parseAll(argc, argv) == HELP_RETURN) return 0;
	
    parser.kickOff(argv[0]);
	
	createLattice(length, beta, extrafield, seed, init_mode, infilename + ".txt");
	
	fs::current_path(fs::current_path() / "measures");
	std::string folder;
	if(outfoldername != "") folder = outfoldername;
	else folder = currentTimeDate();
	fs::create_directory(folder);
	fs::current_path(fs::current_path() / folder);
	fs::create_directory("configurations");
	fs::create_directory("snapshots");
	
	int lognmeas = log10(nmeas) + 1;
	if(do_snapshots) {
		auto curr_path = fs::current_path();
		fs::current_path(fs::current_path() / "snapshots");
		std::string snapname = "snapshot" + makeFixedLength(0, lognmeas) + ".png";
		snapshot(snapname);
		fs::current_path(curr_path);
	}
	
	std::ofstream measfile;
	if (append) measfile.open(measfilename + ".txt", std::fstream::out | std::fstream::app);
	else {
		measfile.open(measfilename + ".txt", std::fstream::out);
		measfile << "# nm = " << nmeas << " nc = " << ncycles << " len = " << length << " extf = " << extrafield << " beta = " << beta << std::endl;
		measfile << "#energy \tmagnetization \tacceptance" << std::endl;
	}
	int percent = 0;
	
	for(int ii = 0; ii < nskip; ii++) {
		for(int jj = 0; jj < ncycles; jj++) for(int kk = 0; kk < length * length; kk ++) updateMetropolis();
	}
	for(int ii = 0; ii < nmeas; ii++) {
		printPercent(ii, percent, nmeas);
		int acc = 0;
		
		for(int jj = 0; jj < ncycles; jj++) for(int kk = 0; kk < length * length; kk ++) acc += updateMetropolis();
		calculateEnergyMagnetization();
		measfile << energy << '\t' << magnetization << '\t' << 1. * acc / ncycles / (length * length) << '\n';
		
		if(do_snapshots && ii % do_snapshots == 0) {
			auto curr_path = fs::current_path();
			fs::current_path(fs::current_path() / "snapshots");
			std::string snapname_ = "snapshot" + makeFixedLength(ii + 1, lognmeas) + ".png";
			snapshot(snapname_);
			fs::current_path(curr_path / "configurations");
			std::string outname_ = outfilename + "-" + makeFixedLength(ii + 1, lognmeas) + ".txt";
			save(outname_);
			fs::current_path(curr_path);
		}
	}
	
	measfile.close();
	save(outfilename + ".txt");
	
	std::cout << std::endl;
	std::cout << "Measure saved in folder " << folder << std::endl;
	return 0;
}