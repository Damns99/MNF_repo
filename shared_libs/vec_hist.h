#ifndef VEC_HIST_H
#define VEC_HIST_H

#include <vector>

namespace vecHist {
	
	template <typename T>
	std::vector<long> makeHist(const std::vector<T>& inputs, unsigned nbins, T a, T b) {
		std::vector<long> results(nbins, 0);
		for(auto element: inputs) {
			int hash = 1. * (element - a) / (b - a) * nbins;
			if(hash >= 0 && hash < nbins) results[hash] += 1;
		}
		return results;
	}
	
	template <typename T>
	std::vector<double> makeHist(const std::vector<T>& inputs, const std::vector<T>& weights, unsigned nbins, T a, T b) {
		std::vector<double> results;
		results.reserve(nbins);
		for(auto element = inputs.begin(), weight = weights.begin(); element != inputs.end() && weight != weights.end(); element++, weight++) {
			int hash = 1. * (*element - a) / (b - a) * nbins;
			if(hash >= 0 && hash < nbins) results[hash] += *weight;
		}
		return results;
	}
}

#endif