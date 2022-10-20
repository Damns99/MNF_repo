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

constexpr int JACC = 1;

constexpr int HOT = 1;
constexpr int COLD = 2;
constexpr int FROM_FILE = 3;

inline int index1d(int i, int L) {
    return ((i %= L) < 0) ? i + L : i;
}
inline int index2d(int i, int j, int L) {
	return L * index1d(i, L) + index1d(j, L);
}

inline int node_type(int i, int j) {
	int k = j%7;
	int jj1 = (6*k)%14, jj2 = (6*k+3)%14, ii = i%14;
	if((jj1 == ii) || (jj2 == ii)) return 0;
	return 1;
}

inline int spin[MAX_LENGTH];
inline int length;
inline double beta;
inline double extrafield;
inline double magnetization;
inline double energy;
inline ran2::RandomGenerator gen;



int lnkd_spin(int i, int j, int L, int link);
	
void latticeFromFile(std::string input_file);

void createLattice(int len, double inv_temp, double extf, long seed, int init_mode, const std::string init_file);

void spinInit(int init_mode, const std::string init_file);

void calculateEnergyMagnetization();

void spinToFile(const std::string output_file);

void spinFromFile(const std::string init_file);

void save(const std::string output_file);

void load(const std::string input_file);

void snapshot(const std::string out_image);

int updateMetropolis();

#endif