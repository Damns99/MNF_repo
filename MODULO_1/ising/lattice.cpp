#include "lattice.h"

Lattice2D::Lattice2D() {;}
		
Lattice2D::Lattice2D(std::string input_file) {
	load(input_file);
	calculateEnergyMagnetizaton();
}
		
Lattice2D::Lattice2D(int len, int geom, double inv_temp, double extf, long seed, int init_mode, const std::string init_file = "") {
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
		
void Lattice2D::geomInit(int geom) {
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
		
void Lattice2D::spinInit(int init_mode, const std::string init_file = "") {
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
		
void Lattice2D::calculateEnergyMagnetizaton() {
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
		
void Lattice2D::spinToFile(const std::string output_file) {
	std::ofstream outfile;
	outfile.open(output_file, std::fstream::out);
	for(int i = 0; i < length * length; i++) outfile << spin[i] << " ";
	outfile << std::endl;
	outfile.close();
}
		
void Lattice2D::spinFromFile(const std::string init_file) {
	std::ifstream infile;
	infile.open(init_file, std::fstream::in);
	for(int i = 0; i < length * length && infile.peek() != EOF; i++) infile >> spin[i];
	infile.close();
}
		
void Lattice2D::save(const std::string output_file) {
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
		
void Lattice2D::load(const std::string input_file) {
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
		
int Lattice2D::updateMetropolis() {
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