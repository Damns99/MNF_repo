#include <iostream>
#include <vector>
#include <math.h>
#include <string>
#include <cassert>
#include <fstream>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <filesystem>
namespace fs = std::filesystem;

#include "cmdline_parser.h"
#include "ran2.h"
#include "text_io.h"

#include <TCanvas.h>
#include <TAxis.h>
#include <TFrame.h>
#include <TLegend.h>
#include <TH2I.h>
#include <TH1.h>
#include <TError.h>

constexpr int MAX_LENGTH = 65536; // = 2^16

constexpr int KB = 1;
constexpr int JACC = 1;

constexpr int HOT = 1;
constexpr int COLD = 2;
constexpr int FROM_FILE = 3;

constexpr int SQUARE_GEOM = 0;
constexpr int EXAGON_GEOM = 1;
constexpr int TRIANG_GEOM = 2;

#define pacman(I, L) (((I) >= 0) ? ((I) % (L)) : ((L) - (((L) - (I)) % (L))))

#define index1d(I, L) (pacman((I), (L)))
#define index2d(I, J, L) ((L) * index1d((I), (L)) + index1d((J), (L)))
#define index3d(I, J, K, L) ((L) * index2d((I), (J), (L)) + index1d((K), (L)))

#define linked_spin(V, I, D) (V).spin[(V).links[(V).links_per_spin * (I) + (D)]]

class Lattice2D {
	public:
		int spin[MAX_LENGTH];
		int links[MAX_LENGTH];
		int length;
		int links_per_spin;
		double beta;
		double extrafield;
		double magnetization;
		double energy;
		ran2::RandomGenerator gen;
		
		Lattice2D() {;}
		
		Lattice2D(std::string input_file) {
			load(input_file);
			calculateEnergyMagnetizaton();
		}
		
		Lattice2D(int len, int geom, double inv_temp, double extf, long seed, int init_mode, const std::string init_file = "") {
			length = len;
			assert(length > 0 && length * length <= MAX_LENGTH);
			
			gen = ran2::RandomGenerator(seed);
			
			geomInit(geom);
			
			beta = inv_temp;
			assert(beta >= 0.);
			extrafield = extf;
			
			spinInit(init_mode, init_file);
			
			calculateEnergyMagnetizaton();
		}
		
		void geomInit(int geom) {
			switch(geom) {
				case SQUARE_GEOM:
					links_per_spin = 4;
					assert(length * length * links_per_spin <= MAX_LENGTH);
					for(int i = 0; i < length; i++) {
						for(int j = 0; j < length; j++) {
							int tmp = links_per_spin * index2d(i, j, length);
							links[tmp + 0] = index2d(i - 1, j, length);
							links[tmp + 1] = index2d(i + 1, j, length);
							links[tmp + 2] = index2d(i, j - 1, length);
							links[tmp + 3] = index2d(i, j + 1, length);
						}
					}
					break;
				case EXAGON_GEOM:
					links_per_spin = 3;
					assert(length * length * links_per_spin <= MAX_LENGTH);
					for(int i = 0; i < length; i++) {
						for(int j = 0; j < length; j++) {
							int tmp = links_per_spin * index2d(i, j, length);
							links[tmp + 0] = index2d(i - 1, j, length);
							links[tmp + 1] = index2d(i + 1, j, length);
							links[tmp + 2] = index2d(i, j - 1 + 2 * ((i + j) % 2), length);
						}
					}
					break;
				case TRIANG_GEOM:
					links_per_spin = 6;
					assert(length * length * links_per_spin <= MAX_LENGTH);
					for(int i = 0; i < length; i++) {
						for(int j = 0; j < length; j++) {
							int tmp = links_per_spin * index2d(i, j, length);
							links[tmp + 0] = index2d(i - 1, j, length);
							links[tmp + 1] = index2d(i + 1, j, length);
							links[tmp + 2] = index2d(i, j - 1, length);
							links[tmp + 3] = index2d(i, j + 1, length);
							links[tmp + 4] = index2d(i - 1, j - 1 + 2 * (i % 2), length);
							links[tmp + 5] = index2d(i + 1, j - 1 + 2 * (i % 2), length);
						}
					}
					break;
				default:
					links_per_spin = 0;
					break;
			}
			
		}
		
		void spinInit(int init_mode, const std::string init_file = "") {
			switch(init_mode) {
				case HOT:
					for(int i = 0; i < length * length; i++) spin[i] = 2 * gen.randL(0, 2) - 1;
					break;
				case COLD:
					for(int i = 0; i < length * length; i++) spin[i] = 1;
					break;
				case FROM_FILE:
					spinFromFile(init_file);
					break;
				default:
					break;
			}
		}
		
		void calculateEnergyMagnetizaton() {
			energy = 0.;
			magnetization = 0.;
			for (int i = 0; i < length * length; i++) {
				for (int j = 0; j < links_per_spin; j++) {
					energy += - JACC / 2. * spin[i] * spin[links[links_per_spin * i + j]] / (1. * length * length);
				}
				energy += - spin[i] * extrafield / (1. * length * length);
				magnetization += spin[i] / (1. * length * length);
			}
		}
		
		void spinToFile(const std::string output_file) {
			std::ofstream outfile;
			outfile.open(output_file, std::fstream::out);
			for(int i = 0; i < length * length; i++) outfile << spin[i] << " ";
			outfile << std::endl;
			outfile.close();
		}
		
		void spinFromFile(const std::string init_file) {
			std::ifstream infile;
			infile.open(init_file, std::fstream::in);
			for(int i = 0; i < length * length && infile.peek() != EOF; i++) infile >> spin[i];
			infile.close();
		}
		
		void save(const std::string output_file) {
			std::ofstream outfile;
			outfile.open(output_file, std::fstream::out);
			outfile << length << std::endl;
			for(int i = 0; i < length * length; i++) outfile << spin[i] << " ";
			outfile << std::endl << links_per_spin << std::endl;
			for(int i = 0; i < length * length * links_per_spin; i++) outfile << links[i] << " ";
			outfile << std::endl << beta << " " << extrafield << std::endl;
			outfile.close();
			gen.toFile("rand_" + output_file);
		}
		
		void load(const std::string input_file) {
			std::ifstream infile;
			infile.open(input_file, std::fstream::in);
			infile >> length;
			for(int i = 0; i < length * length; i++) infile >> spin[i];
			infile >> links_per_spin;
			for(int i = 0; i < length * length * links_per_spin; i++) infile >> links[i];
			infile >> beta >> extrafield;
			infile.close();
			gen.fromFile("rand_" + input_file);
			calculateEnergyMagnetizaton();
		}
		
		void snapshot(const std::string out_image = "snapshot.png") {
			std::stringstream title;
			title << "beta = " << beta << "    ener = " << energy << "    magn = " << magnetization;
			auto c1 = new TCanvas("c1", out_image.c_str(), 1000, 1000);
			auto h1 = new TH2I("h1", title.str().c_str(), length, 0, length, length, 0, length);
			for(int i = 0; i < length; i++) for(int j = 0; j < length; j++) h1->Fill(i, j, spin[index2d(i, j, length)]);
			h1->SetBit(TH1::kNoStats);
			h1->SetMinimum(-1);
			h1->SetMaximum(+1);
			h1->Draw("COL");
			c1->SaveAs(out_image.c_str());
			delete c1;
			delete h1;
		}
		
		int updateMetropolis() {
			int x = gen.randL(0, length * length);
			int newspin = - spin[x];
			double denergy = -2. * newspin * extrafield;
			for (int j = 0; j < links_per_spin; j++) {
				denergy += -2. * JACC * newspin * spin[links[links_per_spin * x + j]];
			}
			double r = exp(-beta * denergy);
			double rr = gen.randF();
			if (rr < r) {
				spin[x] = newspin;
				energy += denergy / (1. * length * length);
				magnetization += 2. * newspin / (1. * length * length);
				return 1;
			}
			return 0;
		}
		
};

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

int main(int argc, char* argv[]) {
	gErrorIgnoreLevel = kWarning;
	
    std::string measfilename, outfilename, infilename;
    long seed;
    int nmeas, ncycles, length, nlinks, geom, init_mode, append, do_snaphots;
    double extrafield, beta;
	
	cmdlineParser::CmdlineParser parser;
    parser.addPosParameter<int>("init_mode", &init_mode, HOT, "[int] Init mode for the first lattice state. HOT = 1, COLD = 2, FROM_FILE = 3");
    parser.addPosParameter<std::string>("measfile", &measfilename, "metro_ising_meas.txt", "[std::string] Output file for the measures.");
    parser.addPosParameter<std::string>("outfile", &outfilename, "metro_ising_conf.txt", "[std::string] Output file for the last iteration lattice state.");
	parser.addOptParameter<std::string>("infile", &infilename, "metro_ising_conf.txt", "[std::string] Intput file for the first iteration lattice state.");
    parser.addOptParameter<long>("seed", &seed, -42, "[long] Random number generator (Ran2) seed.");
    parser.addOptParameter<int>("nm", &nmeas, 1024, "[int] Number of measures to take.");
    parser.addOptParameter<int>("nc", &ncycles, 128, "[int] Number of cycles of iteration between each measure.");
    parser.addOptParameter<int>("len", &length, 10, "[int] Side length of the square lattice.");
    parser.addOptParameter<double>("extf", &extrafield, 0., "[double] Extern magnetic field (in units of JACC).");
    parser.addOptParameter<double>("beta", &beta, 0.3, "[double] One over temperature of the system (in units of JACC).");
    parser.addOptParameter<int>("append", &append, 0, "[int] If != 0 append measures to measfile instead of overwrtiting them.");
    parser.addOptParameter<int>("geom", &geom, SQUARE_GEOM, "[int] Select lattice geomerty: SQUARE = 0, EXAGONAL = 1, TRIANGLUAR = 2");
    parser.addOptParameter<int>("snapshots", &do_snaphots, 0, "[int] If != 0 takes snapshots of the lattice at each measure.");
	
    if (parser.parseAll(argc, argv) == HELP_RETURN) return 0;
	
    parser.kickOff(argv[0]);
	
	auto lattice = Lattice2D(length, geom, beta, extrafield, seed, init_mode, infilename);
	
	std::string folder = "measures" + currentTimeDate();
	fs::create_directory(folder);
	fs::current_path(fs::current_path() / folder);
	
	int lognmeas = log10(nmeas) + 1;
	if(do_snaphots) {
		std::string snapname = "snapshot" + makeFixedLength(0, lognmeas) + ".png";
		lattice.snapshot(snapname);
	}
	
	std::ofstream measfile;
	if (append) measfile.open(measfilename, std::fstream::out | std::fstream::app);
	else {
		measfile.open(measfilename, std::fstream::out);
		measfile << "#energy \tmagnetization \tacceptance" << std::endl;
	}
	int percent = 0;
	for(int ii = 0; ii < nmeas; ii++) {
		if(100 * ii / nmeas > percent) {
			percent = 100 * ii / nmeas;
			std::cout << "\r" << percent << "%";
			std::flush(std::cout);
		}
		int acc = 0;
		for(int jj = 0; jj < ncycles; jj++) acc += lattice.updateMetropolis();
		measfile << lattice.energy << '\t' << lattice.magnetization << '\t' << 1. * acc / ncycles << '\n';
		if(do_snaphots) {
			std::string snapname = "snapshot" + makeFixedLength(ii + 1, lognmeas) + ".png";
			lattice.snapshot(snapname);
		}
	}
	measfile.close();
	lattice.save(outfilename);
	
	std::cout << std::endl;
	return 0;
}