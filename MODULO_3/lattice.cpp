#include "lattice.h"
		
void latticeFromFile(std::string input_file) {
	load(input_file);
	calculateObs();
}
		
void createLattice(int len, int n_part, double inv_temp, long seed, int init_mode, const std::string init_file = "") {
	length = len;
	assert(length > 0 && 2 * length <= MAX_LENGTH);
	nparticles = n_part;
	assert(nparticles > 0 && nparticles <= MAX_PARTICLES && nparticles % length == 0);
	p_length = length / nparticles;
	
	gen = ran2::RandomGenerator(seed);
	
	beta = inv_temp;
	assert(beta >= 0.);
	
	yInit(init_mode, init_file);
	
	calculateObs();
}
		
void yInit(int init_mode, const std::string init_file = "") {
	for(int i = 0; i < length; i++) {
		links[2 * i + 0] = index1d((i % p_length) - 1, p_length) + ((i / p_length) * p_length);
		links[2 * i + 1] = index1d((i % p_length) + 1, p_length) + ((i / p_length) * p_length);
	}
	switch(init_mode) {
		case HOT:
			for(int i = 0; i < length; i++) y[i] = 2. * gen.randF() - 1;
			break;
		case COLD:
			for(int i = 0; i < length; i++) y[i] = 0.;
			break;
		case FROM_FILE:
			yFromFile(init_file);
			break;
		default:
			break;
	}
}
		
void calculateObs() {
	for (int i = 0; i < nparticles; i++) {
		obs1[i] = 0.;
		obs2[i] = 0.;
		for (int j = 0; j < p_length; j++) {
			int x = p_length * i + j;
			obs1[i] += y[x] * y[x];
			obs2[i] += (y[x] - y[links[2 * x + 0]]) * (y[x] - y[links[2 * x + 0]]);
		}
	}
}
		
void yToFile(const std::string output_file) {
	std::ofstream outfile;
	outfile.open(output_file, std::fstream::out);
	for(int i = 0; i < length; i++) outfile << y[i] << " ";
	outfile << std::endl;
	outfile.close();
}
		
void yFromFile(const std::string init_file) {
	std::ifstream infile;
	infile.open(init_file, std::fstream::in);
	for(int i = 0; i < length && infile.peek() != EOF; i++) infile >> y[i];
	infile.close();
}
		
void save(const std::string output_file) {
	std::ofstream outfile;
	outfile.open(output_file, std::fstream::out);
	outfile << length << std::endl;
	for(int i = 0; i < length; i++) outfile << y[i] << " ";
	outfile << std::endl << nparticles << std::endl;
	for(int i = 0; i < length * 2; i++) outfile << links[i] << " ";
	outfile << std::endl << beta << std::endl;
	outfile.close();
	gen.toFile("rand_" + output_file);
}
		
void load(const std::string input_file) {
	std::ifstream infile;
	infile.open(input_file, std::fstream::in);
	infile >> length;
	for(int i = 0; i < length; i++) infile >> y[i];
	infile >> nparticles;
	p_length = length / nparticles;
	for(int i = 0; i < length * 2; i++) infile >> links[i];
	infile >> beta;
	infile.close();
	gen.fromFile("rand_" + input_file);
	calculateObs();
}

void addRule(Rule r, int rep) {
	assert(nrules + 1 <= MAX_RULES);
	rules[nrules] = r;
	repetitions[nrules] = rep;
	nrules++;
}

void update() {
	for(int ii = 0; ii < nrules; ii++) for(int jj = 0; jj < repetitions[ii]; jj++) rules[ii]();
}
