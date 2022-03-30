#include <iostream>
#include <vector>
#include <chrono>
#include "ran2.h"
#include <math.h>

#define SEED -42
#define LENGTH 128
#define REPETITIONS 10000000
#define LOOPS 10

#define pacman(I, L) (((I) >= 0) ? ((I) % (L)) : ((L) - (((L) - (I)) % (L))))

#define index1d(I, L) pacman((I), (L))
#define index2d(I, J, L) (L) * index1d((I), (L)) + index1d((J), (L))
#define index3d(I, J, K, L) (L) * index2d((I), (J), (L)) + index1d((K), (L))

#define tolink(I, J, L) (I) * (L) + (J)

class my_vector {
	public:
		std::vector<double> data;
		my_vector(int n, double el) {
			data.reserve(n);
			for(int i = 0; i < n; i++) data.push_back(el);
		}
};

class lattice {
	public:
		std::vector<double> data;
		std::vector<int> links;
		int links_per_element;
		lattice(int n, double el) {
			data.reserve(n);
			links_per_element = 2;
			links.reserve(links_per_element * n);
			for(int i = 0; i < n; i++) {
				data.push_back(el);
				links.push_back(index1d(i - 1, n));
				links.push_back(index1d(i + 1, n));
			}
		}
		double linked_data(int i, int direction) {return data[links[links_per_element * i + direction]];}
};

int main() {
	{
		double m = 0., m2 = 0.;
		for (int iii = 0; iii < LOOPS; iii++) {
			double vec[LENGTH];
			for(int i = 0; i < LENGTH; i++) vec[i] = 1.;
			ran2::RandomGenerator gen(SEED);
			const auto start = std::chrono::steady_clock::now();
			for(int i = 0; i < REPETITIONS; i++) {
				int r = gen.randL(0, LENGTH);
				vec[r] = vec[r] + (vec[index1d(r - 1, LENGTH)] + vec[index1d(r + 1, LENGTH)]) / 100.;
			}
			const auto end = std::chrono::steady_clock::now();
			double elapsed = (end - start).count() / 1000000.;
			m += elapsed / LOOPS;
			m2 += elapsed * elapsed / LOOPS;
		}
		std::cout << "double[], l=" << LENGTH << ", rep=" << REPETITIONS << " : \t" << m << " +- " << sqrt(m2 - m * m) / sqrt(LOOPS - 1) << " ms" << std::endl;
	}
	{
		double m = 0., m2 = 0.;
		for (int iii = 0; iii < LOOPS; iii++) {
			std::vector<double> vec(LENGTH, 1.);
			ran2::RandomGenerator gen(SEED);
			const auto start = std::chrono::steady_clock::now();
			for(int i = 0; i < REPETITIONS; i++) {
				int r = gen.randL(0, LENGTH);
				vec[r] = vec[r] + (vec[index1d(r - 1, LENGTH)] + vec[index1d(r + 1, LENGTH)]) / 100.;
			}
			const auto end = std::chrono::steady_clock::now();
			double elapsed = (end - start).count() / 1000000.;
			m += elapsed / LOOPS;
			m2 += elapsed * elapsed / LOOPS;
		}
		std::cout << "vector, l=" << LENGTH << ", rep=" << REPETITIONS << " : \t\t" << m << " +- " << sqrt(m2 - m * m) / sqrt(LOOPS - 1) << " ms" << std::endl;
	}
	{
		double m = 0., m2 = 0.;
		for (int iii = 0; iii < LOOPS; iii++) {
			my_vector vec(LENGTH, 1.);
			ran2::RandomGenerator gen(SEED);
			const auto start = std::chrono::steady_clock::now();
			for(int i = 0; i < REPETITIONS; i++) {
				int r = gen.randL(0, LENGTH);
				vec.data[r] = vec.data[r] + (vec.data[index1d(r - 1, LENGTH)] + vec.data[index1d(r + 1, LENGTH)]) / 100.;
			}
			const auto end = std::chrono::steady_clock::now();
			double elapsed = (end - start).count() / 1000000.;
			m += elapsed / LOOPS;
			m2 += elapsed * elapsed / LOOPS;
		}
		std::cout << "my_vector, l=" << LENGTH << ", rep=" << REPETITIONS << " : \t" << m << " +- " << sqrt(m2 - m * m) / sqrt(LOOPS - 1) << " ms" << std::endl;
	}
	{
		double m = 0., m2 = 0.;
		for (int iii = 0; iii < LOOPS; iii++) {
			lattice vec(LENGTH, 1.);
			ran2::RandomGenerator gen(SEED);
			const auto start = std::chrono::steady_clock::now();
			for(int i = 0; i < REPETITIONS; i++) {
				int r = gen.randL(0, LENGTH);
				vec.data[r] = vec.data[r] + (vec.linked_data(r, 0) + vec.linked_data(r, 1)) / 100.;
			}
			const auto end = std::chrono::steady_clock::now();
			double elapsed = (end - start).count() / 1000000.;
			m += elapsed / LOOPS;
			m2 += elapsed * elapsed / LOOPS;
		}
		std::cout << "lattice, l=" << LENGTH << ", rep=" << REPETITIONS << " : \t\t" << m << " +- " << sqrt(m2 - m * m) / sqrt(LOOPS - 1) << " ms" << std::endl;
	}
	{
		double m = 0., m2 = 0.;
		for (int iii = 0; iii < LOOPS; iii++) {
			lattice vec(LENGTH, 1.);
			ran2::RandomGenerator gen(SEED);
			const auto start = std::chrono::steady_clock::now();
			for(int i = 0; i < REPETITIONS; i++) {
				int r = gen.randL(0, LENGTH);
				vec.data[r] = vec.data[r] + (vec.data[vec.links[vec.links_per_element * r + 0]] + vec.data[vec.links[vec.links_per_element * r + 1]]) / 100.;
			}
			const auto end = std::chrono::steady_clock::now();
			double elapsed = (end - start).count() / 1000000.;
			m += elapsed / LOOPS;
			m2 += elapsed * elapsed / LOOPS;
		}
		std::cout << "lattice faster, l=" << LENGTH << ", rep=" << REPETITIONS << " : \t" << m << " +- " << sqrt(m2 - m * m) / sqrt(LOOPS - 1) << " ms" << std::endl;
	}
	{
		double m = 0., m2 = 0.;
		for (int iii = 0; iii < LOOPS; iii++) {
			std::vector<double> vec(LENGTH, 1.);
			std::vector<int> links;
			int links_per_element = 2;
			links.reserve(links_per_element * LENGTH);
			for(int i = 0; i < LENGTH; i++) {
				links.push_back(index1d(i - 1, LENGTH));
				links.push_back(index1d(i + 1, LENGTH));
			}
			ran2::RandomGenerator gen(SEED);
			const auto start = std::chrono::steady_clock::now();
			for(int i = 0; i < REPETITIONS; i++) {
				int r = gen.randL(0, LENGTH);
				vec[r] = vec[r] + (vec[links[tolink(r, 0, links_per_element)]] + vec[links[tolink(r, 1, links_per_element)]]) / 100.;
			}
			const auto end = std::chrono::steady_clock::now();
			double elapsed = (end - start).count() / 1000000.;
			m += elapsed / LOOPS;
			m2 += elapsed * elapsed / LOOPS;
		}
		std::cout << "vector + links, l=" << LENGTH << ", rep=" << REPETITIONS << " : \t" << m << " +- " << sqrt(m2 - m * m) / sqrt(LOOPS - 1) << " ms" << std::endl;
	}
	{
		double m = 0., m2 = 0.;
		for (int iii = 0; iii < LOOPS; iii++) {
			std::vector<double> vec(LENGTH, 1.);
			std::vector<double*> links;
			int links_per_element = 2;
			links.reserve(links_per_element * LENGTH);
			for(int i = 0; i < LENGTH; i++) {
				links.push_back(&(*(vec.begin() + index1d(i - 1, LENGTH))));
				links.push_back(&(*(vec.begin() + index1d(i + 1, LENGTH))));
			}
			ran2::RandomGenerator gen(SEED);
			const auto start = std::chrono::steady_clock::now();
			for(int i = 0; i < REPETITIONS; i++) {
				int r = gen.randL(0, LENGTH);
				vec[r] = vec[r] + (*links[tolink(r, 0, links_per_element)] + *links[tolink(r, 1, links_per_element)]) / 100.;
			}
			const auto end = std::chrono::steady_clock::now();
			double elapsed = (end - start).count() / 1000000.;
			m += elapsed / LOOPS;
			m2 += elapsed * elapsed / LOOPS;
		}
		std::cout << "pointer links, l=" << LENGTH << ", rep=" << REPETITIONS << " : \t" << m << " +- " << sqrt(m2 - m * m) / sqrt(LOOPS - 1) << " ms" << std::endl;
	}
	std::cout << "I choose lattice faster." << std::endl;
	std::cout << std::endl;
}