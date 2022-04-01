#include <iostream>
#include <vector>
#include <math.h>
#include <string>
#include <cassert>
#include <fstream>

#include "cmdline_parser.h"
#include "ran2.h"
#include "text_io.h"

#include <TCanvas.h>
#include <TAxis.h>
#include <TFrame.h>
#include <TLegend.h>
#include <TH2I.h>
#include <TH1.h>

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

#define index1d(I, L) pacman((I), (L))
#define index2d(I, J, L) (L) * index1d((I), (L)) + index1d((J), (L))
#define index3d(I, J, K, L) (L) * index2d((I), (J), (L)) + index1d((K), (L))

#define linked_spin(V, I, D) (V).spin[(V).links[(V).links_per_spin * (I) + (D)]]

class Lattice2D {
	public:
		double spin[MAX_LENGTH];
		int links[MAX_LENGTH];
		int length;
		int links_per_spin;
		double temperature;
		double extrafield;
		double magnetization;
		double energy;
		ran2::RandomGenerator gen;
		
		Lattice2D() {;}
		
		Lattice2D(std::string input_file) {
			load(input_file);
			
			// calculate energy magnetization
		}
		
		Lattice2D(int len, int geom, double temp, double extf, long seed, int init_mode, const std::string init_file = "") {
			length = len;
			assert(length > 0 && length * length <= MAX_LENGTH);
			
			gen = ran2::RandomGenerator(seed);
			
			geomInit(geom);
			
			temperature = temp;
			assert(temperature >= 0.);
			extrafield = extf;
			
			spinInit(init_mode, init_file);
			
			// calculate energy magnetization
		}
		
		void geomInit(int geom) {
			switch(geom) {
				case SQUARE_GEOM:
					links_per_spin = 4;
					assert(length * length * links_per_spin <= MAX_LENGTH);
					for(int i = 0; i < length; i++) {
						for(int j = 0; j < length; j++) {
							links[links_per_spin * index2d(i, j, length) + 0] = index2d(i - 1, j, length);
							links[links_per_spin * index2d(i, j, length) + 1] = index2d(i + 1, j, length);
							links[links_per_spin * index2d(i, j, length) + 2] = index2d(i, j - 1, length);
							links[links_per_spin * index2d(i, j, length) + 3] = index2d(i, j + 1, length);
						}
					}
					break;
				case EXAGON_GEOM:
					links_per_spin = 3;
					assert(length * length * links_per_spin <= MAX_LENGTH);
					for(int i = 0; i < length; i++) {
						for(int j = 0; j < length; j++) {
							links[links_per_spin * index2d(i, j, length) + 0] = index2d(i - 1, j, length);
							links[links_per_spin * index2d(i, j, length) + 1] = index2d(i + 1, j, length);
							links[links_per_spin * index2d(i, j, length) + 2] = index2d(i, j - 1 + 2 * ((i + j) % 2), length);
						}
					}
					break;
				case TRIANG_GEOM:
					links_per_spin = 6;
					assert(length * length * links_per_spin <= MAX_LENGTH);
					for(int i = 0; i < length; i++) {
						for(int j = 0; j < length; j++) {
							links[links_per_spin * index2d(i, j, length) + 0] = index2d(i - 1, j, length);
							links[links_per_spin * index2d(i, j, length) + 1] = index2d(i + 1, j, length);
							links[links_per_spin * index2d(i, j, length) + 2] = index2d(i, j - 1, length);
							links[links_per_spin * index2d(i, j, length) + 3] = index2d(i, j + 1, length);
							links[links_per_spin * index2d(i, j, length) + 4] = index2d(i - 1, j - 1 + 2 * (i % 2), length);
							links[links_per_spin * index2d(i, j, length) + 5] = index2d(i + 1, j - 1 + 2 * (i % 2), length);
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
			outfile << std::endl << temperature << " " << extrafield << std::endl;
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
			infile >> temperature >> extrafield;
			infile.close();
			gen.fromFile("rand_" + input_file);
			
			// calculate energy magnetization
		}
		
		void snapshot(const std::string out_image = "snapshot.png") {
			auto c1 = new TCanvas("c1","c1",600,400);
			auto hcol1 = new TH2I("hcol1","",length,0,length,length,0,length);
			for(int i = 0; i < length; i++) for(int j = 0; j < length; j++) hcol1->Fill(i, j, spin[index2d(i, j, length)]);
			hcol1->SetBit(TH1::kNoStats);
			hcol1->SetMinimum(-1);
			hcol1->SetMaximum(+1);
			hcol1->Draw("COL");
			c1->SaveAs(out_image.c_str());
		}
		
		// update
		
};

int main(int argc, char* argv[]) {
	
    std::string measfilename, outfilename, infilename;
    long seed;
    int nmeas, ncycles, length, nlinks, geom, init_mode, append;
    double extrafield, temperature;
	
	cmdlineParser::CmdlineParser parser;
    parser.addPosParameter<std::string>("measfile", &measfilename, "metro_ising_meas.txt", "[std::string] Output file for the measures.");
    parser.addPosParameter<std::string>("outfile", &outfilename, "metro_ising_conf.txt", "[std::string] Output file for the last iteration lattice state.");
    parser.addPosParameter<int>("init_mode", &init_mode, COLD, "[int] Init mode for the first lattice state. HOT = 1, COLD = 2, FROM_FILE = 3");
	parser.addOptParameter<std::string>("infile", &infilename, "metro_ising_conf.txt", "[std::string] Intput file for the first iteration lattice state.");
    parser.addOptParameter<long>("seed", &seed, -42, "[long] Random number generator (Ran2) seed.");
    parser.addOptParameter<int>("nm", &nmeas, 1024, "[int] Number of measures to take.");
    parser.addOptParameter<int>("nc", &ncycles, 128, "[int] Number of cycles of iteration between each measure.");
    parser.addOptParameter<int>("len", &length, 10, "[int] Side length of the square lattice.");
    parser.addOptParameter<double>("extf", &extrafield, 0., "[double] Extern magnetic field.");
    parser.addOptParameter<double>("temp", &temperature, 0., "[double] Temperature of the system (in units of extrafield/k_B)");
    parser.addOptParameter<int>("append", &append, 0, "[int] If != 0 append measures to measfile instead of overwrtiting them.");
    parser.addOptParameter<int>("geom", &geom, SQUARE_GEOM, "[int] Select lattice geomerty: SQUARE = 0, EXAGONAL = 1, TRIANGLUAR = 2");
	
    if (parser.parseAll(argc, argv) == HELP_RETURN) return 0;
	
    parser.kickOff(argv[0]);
	
	auto lattice = Lattice2D(length, geom, temperature, extrafield, seed, init_mode, infilename);
	lattice.snapshot();
	
	return 0;
}