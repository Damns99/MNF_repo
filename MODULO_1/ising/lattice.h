#ifndef LATTICE_H
#define LATTICE_H

#include <iostream>
#include <math.h>
#include <string>
#include <cassert>
#include <fstream>
#include <sstream>

#include "ran2.h"

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
		
		Lattice2D();
		
		Lattice2D(std::string input_file);
		
		Lattice2D(int len, int geom, double inv_temp, double extf, long seed, int init_mode, const std::string init_file);
		
		void geomInit(int geom);
		
		void spinInit(int init_mode, const std::string init_file);
		
		void calculateEnergyMagnetizaton();
		
		void spinToFile(const std::string output_file);
		
		void spinFromFile(const std::string init_file);
		
		void save(const std::string output_file);
		
		void load(const std::string input_file);
		
		void snapshot(const std::string out_image);
		
		int updateMetropolis();
		
};

void cudaSquareUpdateMetropolis(Lattice2D* lattice);

#endif