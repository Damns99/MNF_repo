#ifndef BOOTSTRAP_H
#define BOOTSTRAP_H

#include "ran2.h"
#include <math.h>
#include <vector>

namespace bootstrap {
	
	template <typename T, typename... Parameters>
	T bootstrapError(T (*func)(const std::vector<T>&, Parameters...), const std::vector<T>& data, int sample_number, int subsample_length, long seed, Parameters... pars) {
		T result;		
		T s1 = 0., s2 = 0.;
		int data_length = data.size();
		int n_draws = data_length / subsample_length;
		ran2::RandomGenerator gen(seed);
		
		for(int i = 0; i < sample_number; i++) {
			T fakemean = 0.;
			std::vector<T> fakedata;
			
			for(int j = 0; j < n_draws; j++) {
				int r = gen.randL(0, data_length);
				for(int k = 0; k < subsample_length; k++) fakedata.push_back(data[(r + k) % data_length]);
			}
			fakemean = func(fakedata, pars...);
			
			s1 += fakemean * fakemean / sample_number;
			s2 += fakemean / sample_number;
		}
		
		result = sqrt(s1 - s2 * s2);
		return result;
	}
	
	template <typename T, typename... Parameters>
	T bootstrapError(T (*func)(T, Parameters...), const std::vector<T>& data, int sample_number, int subsample_length, long seed, Parameters... pars) {
		T result;		
		T s1 = 0., s2 = 0.;
		int data_length = data.size();
		int n_draws = data_length / subsample_length;
		ran2::RandomGenerator gen(seed);
		
		for(int i = 0; i < sample_number; i++) {
			T fakemean = 0.;
			std::vector<T> fakedata;
			
			for(int j = 0; j < n_draws; j++) {
				int r = gen.randL(0, data_length);
				for(int k = 0; k < subsample_length; k++) fakemean += func(data[(r + k) % data_length], pars...);
			}
			
			s1 += fakemean * fakemean / sample_number;
			s2 += fakemean / sample_number;
		}
		
		result = sqrt(s1 - s2 * s2);
		return result;
	}
	
}

#endif