// RAN2 from Numerical Recipes
#ifndef RAN2_H
#define RAN2_H

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

#include <fstream>
#include <string>

namespace ran2 {
	
	class RandomGenerator {
		private:
			long idum = -123456789;
			long idum2 = 123456789;
			long iy = 0;
			long iv[NTAB];
			long ran2Base();
			
		public:
			RandomGenerator();
			RandomGenerator(long seed);
			void reset(long seed);
			
			float randF();
			float randF(float a, float b);
			long randL();
			long randL(long a, long b);
			
			void fromFile(const std::string filepath);
			void toFile(const std::string filepath);
	};

}

#endif
