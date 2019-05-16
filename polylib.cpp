#include <cmath>
#include "Polylib.h"

double JacP(double alpha, double beta, int n, double x) {
	if (n == 0) { return 1.0; }
	else if (n == 1) { return 0.5*(alpha - beta + (alpha + beta + 2)*x); }
	else if (n > 1) {
		double an1 = 2*((n-1)+1)*((n-1)+alpha+beta+1)*(2*(n-1)+alpha+beta);
		double an2 = (2*(n-1)+alpha+beta+1.)*(pow(alpha, 2) - pow(beta, 2));
		double an3 = (2*(n-1)+alpha+beta)*(2*(n-1)+alpha+beta+1)*(2*(n-1)+alpha+beta+2);
		double an4 = 2*((n-1)+alpha)*((n-1)+beta)*(2*(n-1)+alpha+beta+2);
		return (((an2 + an3*x)*JacP(alpha, beta, n-1, x) - an4*JacP(alpha, beta, n-2, x))/an1);
	}
}

double dJacP(double alpha, double beta, int n, double x) {
	return 0.5*(alpha+beta+n+1.0)*JacP(alpha+1, beta+1, n-1, x);
}

void JacPZ(double alpha, double beta, int n, double *x) {
	double r = 0., s = 0.0, delta = 0., tol = 0.000001;
	const long double PI = 3.14159265358979323846;
	for (int k = 0; k < n; ++k) {
		r = -std::cos(((2.*k+1.)/(2.*n))*PI);
		if (k > 0) r = (r + x[k-1])/2.;
		do {
			s = 0.0;
			for (int i = 0; i < k; ++i) s += 1./(r - x[i]);
			delta = -JacP(alpha, beta, n, r)/(dJacP(alpha, beta, n, r) - s*JacP(alpha, beta, n ,r));
			r += delta;
			if (std::fabs(delta) < tol) break;
		} while (true);
		x[k] = r;
	}
}