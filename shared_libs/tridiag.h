#ifndef TRIDIAG_H
#define TRIDIAG_H

namespace triDiag{
	
	void triDiagSolve(const int N, const double* a, const double* b, const double* c, const double* d, double* x);
	void cTriDiagSolve(const int N, const double* a, const double* b, const double* c, const double* d, double* x);
	
}

#endif