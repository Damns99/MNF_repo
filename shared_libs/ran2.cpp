// RAN2 from Numerical Recipes
#include "ran2.h"

ran2::RandomGenerator::RandomGenerator() {
	RandomGenerator(42);
}

ran2::RandomGenerator::RandomGenerator(long seed) {
	idum = seed < 0 ? seed : -seed;
}

void ran2::RandomGenerator::reset(long seed) {
	idum = seed;
	idum2 = 123456789;
	iy = 0;
	for(auto& ii: iv) ii = 0;
}

long ran2::RandomGenerator::ran2Base() {
	int j;
	long k;

	if (idum <= 0) {
		if (-(idum) < 1) idum=1;
		else idum = -(idum);
		idum2=(idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(idum)/IQ1;
			idum=IA1*(idum-k*IQ1)-k*IR1;
			if (idum < 0) idum += IM1;
			if (j < NTAB) iv[j] = idum;
		}
	iy=iv[0];
	}
	
	k=(idum)/IQ1;
	idum=IA1*(idum-k*IQ1)-k*IR1;
	if (idum < 0) idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = idum;
	if (iy < 1) iy += IMM1;
	return iy;
}

long ran2::RandomGenerator::randL() {
	return ran2Base();
}

long ran2::RandomGenerator::randL(const long a, const long b) {
	return (ran2Base() % (b - a)) + a;
}

float ran2::RandomGenerator::randF() {
	float temp;
	if((temp = AM * ran2Base()) > RNMX) return RNMX;
	else return temp;
}

float ran2::RandomGenerator::randF(const float a, const float b) {
	float temp;
	if((temp = (b - a) / IM1 * ran2Base() + a) > b) return b;
	else return temp;
}

void ran2::RandomGenerator::fromFile(const std::string filepath) {
	std::ifstream infile;
	infile.open(filepath, std::fstream::in);
	infile >> idum >> idum2 >> iy;
	for(auto& ii: iv) infile >> ii;
}

void ran2::RandomGenerator::toFile(const std::string filepath) {
	std::ofstream outfile;
	outfile.open(filepath, std::fstream::out);
	outfile << idum << std::endl << idum2 << std::endl  << iy << std::endl;
	for(auto& ii: iv) outfile << ii << std::endl;
}
