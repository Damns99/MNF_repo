#include "cmdline_parser.h"
#include <iostream>
#include <vector>
#include "ran2.h"
#include <math.h>
#include <string>

#define MAX_LENGTH 65536 // = 2^16

#define KB 1
#define JACC 1

#define HOT 1
#define COLD 2

#define SQUARE_GEOM 0
#define EXAGON_GEOM 1
#define TRIANG_GEOM 2

#define pacman(I, L) (((I) >= 0) ? ((I) % (L)) : ((L) - (((L) - (I)) % (L))))

#define index1d(I, L) pacman((I), (L))
#define index2d(I, J, L) (L) * index1d((I), (L)) + index1d((J), (L))
#define index3d(I, J, K, L) (L) * index2d((I), (J), (L)) + index1d((K), (L))

class Lattice {
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
		
		Lattice(int len, int geom, double temp, double extf, char init[], long* seed) {
			;
		}
};