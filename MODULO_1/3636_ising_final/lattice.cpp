#include "lattice.h"

int lnkd_spin(int i, int j, int L, int link) {
	int type = node_type(i, j);
	switch(type) {
		case 0:
			return 0;
		case 1:
			switch(link) {
				case 0:
					return spin[index2d(i+1, j, L)];
				case 1:
					return spin[index2d(i, j+1, L)];
				case 2:
					return spin[index2d(i-1, j, L)];
				case 3:
					return spin[index2d(i, j-1, L)];
				default:
					return 0;
			}
		case 2:
			switch(link) {
				case 0:
					return spin[index2d(i+1, j, L)];
				case 1:
					return spin[index2d(i+1, j+1, L)];
				case 2:
					return spin[index2d(i-1, j, L)];
				case 3:
					return spin[index2d(i-1, j-1, L)];
				default:
					return 0;
			}
		case 3:
			switch(link) {
				case 0:
					return spin[index2d(i+1, j+1, L)];
				case 1:
					return spin[index2d(i, j+1, L)];
				case 2:
					return spin[index2d(i-1, j-1, L)];
				case 3:
					return spin[index2d(i, j-1, L)];
				default:
					return 0;
			}
		default:
			return 0;
	}
}

void latticeFromFile(std::string input_file) {
	load(input_file);
	calculateEnergyMagnetization();
}

void createLattice(int len, double inv_temp, double extf, long seed, int init_mode, const std::string init_file = "") {
	length = len;
	assert(length > 0 && length * length <= MAX_LENGTH && (length % 2) == 0);
	
	gen = ran2::RandomGenerator(seed);
	
	beta = inv_temp;
	assert(beta >= 0.);
	extrafield = extf;
	
	spinInit(init_mode, init_file);
	
	calculateEnergyMagnetization();
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
	outfile << std::endl << beta << " " << extrafield << std::endl;
	outfile.close();
	gen.toFile("rand_" + output_file);
}

void load(const std::string input_file) {
	std::ifstream infile;
	infile.open(input_file, std::fstream::in);
	infile >> length;
	for(int i = 0; i < length * length; i++) infile >> spin[i];
	infile >> beta >> extrafield;
	infile.close();
	gen.fromFile("rand_" + input_file);
	calculateEnergyMagnetization();
}

void calculateEnergyMagnetization() {
	energy = 0.;
	magnetization = 0.;
	
	for (int i = 0; i < length; i++) for (int j = 0; j < length; j++) {
		int curr_spin = spin[index2d(i,j,length)];
		int neighbours = 0;
		for (int k = 0; k < 4; k++) neighbours += lnkd_spin(i, j, length, k);
		energy += - curr_spin * neighbours / 2.;
		magnetization += curr_spin;
	}
	
	double normalization = 1. / (length * length) * (16. / 9.);
	magnetization *= normalization;
	energy *= JACC * normalization;
	energy += - extrafield * magnetization;
}

int updateMetropolis() {
	int x = gen.randL(0, length * length);
	int newspin = - spin[x];
	int i = x / length, j = x % length;
	int neighbours = 0;
	for (int k = 0; k < 4; k++) neighbours += lnkd_spin(i, j, length, k);
	
	double denergy = -2. * newspin * (JACC * neighbours + extrafield);
	
	double r = exp(-beta * denergy);
	double rr = gen.randF();
	if (rr < r) {	
		spin[x] = newspin;
		return 1;
	}
	return 0;
}