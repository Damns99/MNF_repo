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
constexpr int MAX_RULES = 8;
constexpr int MAX_PARTICLES = 3;

constexpr int HOT = 1;
constexpr int COLD = 2;
constexpr int FROM_FILE = 3;

inline int index1d(int i, int L) {
    return ((i %= L) < 0) ? i + L : i;
}

typedef int (*Rule)();

inline double y[MAX_LENGTH];
inline int links[MAX_LENGTH];
inline int length;
inline int nparticles;
inline int p_length;
inline double beta;
inline ran2::RandomGenerator gen;
inline Rule rules[MAX_RULES];
inline int repetitions[MAX_RULES];
inline int nrules;
inline double obs1[MAX_PARTICLES];
inline double obs2[MAX_PARTICLES];

void latticeFromFile(std::string input_file);

void createLattice(int len, int n_part, double inv_temp, long seed, int init_mode, const std::string init_file);

void yInit(int init_mode, const std::string init_file);

void yToFile(const std::string output_file);

void yFromFile(const std::string init_file);

void save(const std::string output_file);

void load(const std::string input_file);

void snapshot(const std::string out_image);

void calculateObs();

void addRule(Rule r, int rep);

double update();

#endif