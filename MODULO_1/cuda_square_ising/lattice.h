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

#define pacman(I, L) (((I) >= 0) ? ((I) % (L)) : ((L) - (((L) - (I)) % (L))))

#define index1d(I, L) (pacman((I), (L)))
#define index2d(I, J, L) ((L) * index1d((I), (L)) + index1d((J), (L)))

inline int spin[MAX_LENGTH];
inline int length;
inline double beta;
inline double extrafield;
inline double magnetization;
inline double energy;
inline ran2::RandomGenerator gen;
		
void latticeFromFile(std::string input_file);

void createLattice(int len, double inv_temp, double extf, long seed, int init_mode, const std::string init_file);

void spinInit(int init_mode, const std::string init_file);

void calculateEnergyMagnetizaton();

void spinToFile(const std::string output_file);

void spinFromFile(const std::string init_file);

void save(const std::string output_file);

void load(const std::string input_file);

void snapshot(const std::string out_image);

int updateMetropolis();
int updateMetropolis2();

#endif