#include "tridiag.h"

void triDiag::triDiagSolve(const int N, const double* a, const double* b, const double* c, const double* d, double* x) {
	double* p = new double[N]; double* q = new double[N]; double* u = new double[N];
    p[0] = b[0];
    q[0] = - c[0] / p[0];
    u[0] = d[0] / p[0];
    for (int k = 1; k < N; k++) {
        p[k] = a[k] * q[k-1] + b[k];
        q[k] = - c[k] / p[k];
        u[k] = (d[k] - a[k] * u[k-1]) / p[k];
    }
    x[N-1] = u[N-1];
    for (int k = N-2; k >= 0; k--) {
        x[k] = q[k] * x[k+1] + u[k];
    }
	delete[] p; delete[] q; delete[] u;
}

void triDiag::cTriDiagSolve(const int N, const double* a, const double* b, const double* c, const double* d, double* x) {
    double* p = new double[N]; double* q = new double[N]; double* u = new double[N];
    double* s = new double[N]; double* t = new double[N]; double* v = new double[N];
    p[0] = b[0];
    q[0] = - c[0] / p[0];
    u[0] = d[0] / p[0];
    s[0] = - a[0] / p[0];
    for (int k = 1; k < N; k++) {
        p[k] = a[k] * q[k-1] + b[k];
        q[k] = - c[k] / p[k];
        u[k] = (d[k] - a[k] * u[k-1]) / p[k];
        s[k] = - a[k] * s[k-1] / p[k];
    }
    t[N-1] = 1; v[N-1] = 0;
    for (int k = N-2; k >= 0; k--) {
        t[k] = q[k] * t[k+1] + s[k];
        v[k] = q[k] * v[k+1] + u[k];
    }
    x[N-1] = (d[N-1] - c[N-1] * v[0] - a[N-1] * v[N-2]) / (c[N-1] * t[0] + a[N-1] * t[N-2] + b[N-1]);
    for (int k = N-2; k >= 0; k--) {
        x[k] = t[k] * x[N-1] + v[k];
    }
	delete[] p; delete[] q; delete[] u;
	delete[] s; delete[] t; delete[] v;
}
