#include "OneDim.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>

void qsi2x(double *xg, double *loc, double *m, int Nel, int Nsp) {
	for (int i = 0; i != Nel; ++i) {
		for (int j = 0; j != Nsp; ++j) {
			xg[j + i*Nsp] = 0.5*(1.-loc[j])*m[i] + 0.5*(1.+loc[j])*m[i+1];
		}
	}
}

double Lagrange(double x, double* xg, int Nsp, int el, int k) {
	double prod = 1.0;
	for (int l = 0; l != Nsp; ++l) {
		if (l != k) {
			prod *= (x - xg[el*Nsp + l])/(xg[el*Nsp + k] - xg[el*Nsp + l]);
		}
	}
	return prod;
}

double LagrangeLoc(double x, double* loc, int Nsp, int k) {

}

void printSol(double* x, double* u, int tot, std::string fn) {
	static int cnt = 0;

	std::ofstream file("results/" + fn + std::to_string(cnt++) + ".plt");
	file << std::setprecision(6) << std::scientific;
	for (int i = 0; i != tot; ++i)
		file << x[i] << std::setw(15) << u[i] << "\n";
	file.close();
}