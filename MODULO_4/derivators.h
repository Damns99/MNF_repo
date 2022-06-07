#ifndef DERIVATORS_H
#define DERIVATORS_H

#include "bound_cond_vecs.h"

namespace derivators {
	typedef bound_cond_vecs::BoundCondVec<double> (*Derivator)(const bound_cond_vecs::BoundCondVec<double>&, const double);
	
	namespace utilities {
		bound_cond_vecs::BoundCondVec<double> fwd_der_function(const bound_cond_vecs::BoundCondVec<double>& input, const double dx) {
			bound_cond_vecs::BoundCondVec<double> output(input.len());
			for(int ii = 0; ii < input.len(); ii++) output[ii] = (input[ii + 1] - input[ii]) / dx;
			return output;
		}
	}
	Derivator fwd_derive = &(utilities::fwd_der_function);
}

#endif