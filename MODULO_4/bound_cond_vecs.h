#ifndef BOUND_COND_VECS
#define BOUND_COND_VECS

#include <cassert>

#define MAX_LENGTH 65536

#define PERIODIC_BC 101
#define ABSORBANT_BC 102
#define REFLECTIVE_BC 103

namespace bound_cond_vecs {
	template <typename T>
	class BoundCondVec{
		private:
			int length = 0;
			T* vec = nullptr;
			int mode = PERIODIC_BC;
			T placeholder = 0.;
		public:
			BoundCondVec(const int new_length) {
				assert(new_length <= MAX_LENGTH);
				vec = new T[new_length];
				for(int i = 0; i < new_length; i++) vec[i] = 0.;
				length = new_length;
			}
			BoundCondVec(const int new_length, const int new_mode) {
				assert(new_length <= MAX_LENGTH);
				vec = new T[new_length];
				for(int i = 0; i < new_length; i++) vec[i] = 0.;
				length = new_length;
				assert(new_mode == PERIODIC_BC || new_mode == ABSORBANT_BC || new_mode == REFLECTIVE_BC);
				mode = new_mode;
			}
			BoundCondVec(const int new_length, const T* new_vec) {
				assert(new_length <= MAX_LENGTH);
				vec = new T[new_length];
				for(int i = 0; i < new_length; i++) vec[i] = new_vec[i];
				length = new_length;
			}
			BoundCondVec(const int new_length, const T* new_vec, const int new_mode) {
				assert(new_length <= MAX_LENGTH);
				vec = new T[new_length];
				for(int i = 0; i < new_length; i++) vec[i] = new_vec[i];
				length = new_length;
				assert(new_mode == PERIODIC_BC || new_mode == ABSORBANT_BC || new_mode == REFLECTIVE_BC);
				mode = new_mode;
			}
			BoundCondVec(const BoundCondVec & other_bcv) {
				length = other_bcv.length;
				vec = new T[length];
				for(int i = 0; i < length; i++) vec[i] = other_bcv[i];
				mode = other_bcv.mode;
			}
			~BoundCondVec() {
				if(length > 0) delete[] vec;
			}
			
			T& operator[](int index) {
				if(mode == PERIODIC_BC) {
					return vec[(index - length * ((index / length) - (index < 0))) % length];
				}
				else if(mode == ABSORBANT_BC) {
					if(index >= 0 && index < length) return vec[index];
					else return placeholder;
				}
				else if(mode == REFLECTIVE_BC) {
					int m = (index + (index < 0)) / length;
					int r = (m + (index < 0)) % 2;
					int q = (index - length * (m - (index < 0))) % length;
					return vec[(1 - r) * q + r * (length - 1 - q)];
				}
				else return placeholder;
			}
			T operator [] (int index) const {
				if(mode == PERIODIC_BC) {
					return vec[(index - length * ((index / length) - (index < 0))) % length];
				}
				else if(mode == ABSORBANT_BC) {
					if(index >= 0 && index < length) return vec[index];
					else return placeholder;
				}
				else if(mode == REFLECTIVE_BC) {
					int m = (index + (index < 0)) / length;
					int r = (m + (index < 0)) % 2;
					int q = (index - length * (m - (index < 0))) % length;
					return vec[(1 - r) * q + r * (length - 1 - q)];
				}
				else return placeholder;
			}
			
			T* data() const {
				return vec;
			}
			int len() const {
				return length;
			}
			int getMode() const {
				return mode;
			}
	};
}

#endif