#ifndef BOX_MULLER_H
#define BOX_MULLER_H

#include <fstream>
#include <math.h>
#include <string>
#include <vector>

#include "ran2.h"

namespace boxMuller {
	
	class BoxMuller {
		private:
			ran2::RandomGenerator r = ran2::RandomGenerator(-42);
		public:
			BoxMuller(long seed);
			void reset(long seed);
			void rand(float mean, float stdev, float& r1, float& r2);
			std::vector<float> randToVec(float mean, float stdev, int n);
			void randToFile(float mean, float stdev, int n, const std::string filepath);
	};
	
}

#endif