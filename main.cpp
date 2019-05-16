#include <iostream>
#include <string>
#include <cmath>
#include <iomanip>
#include "Polylib.h"
#include "Polykarn.h"
#include "OneDim.h"

#define PI 3.14159265358979323846

int main(int argc, char** argv) {

	if (argc < 2) { 
		std::cout << "One Dimensional solver using High-Order FR method\n\n"
				  << "Usage:\t./test\t<N. of elements>  [N. of solution points]  [Solution points (type)]\n\nwhere:\n"
				  << "  Solution point types:"
				  << "	\t0 - GAUSS (default)\n"
				  << "  			 \t1 - LOBATTO\n"
				  << "  			 \t2 - RADAU\n";
		return EXIT_FAILURE;
	}
	size_t Nel, Nsp = 4;
	SP spt = GAUSS;
	for (int i = 1; i != argc; ++i) {
		switch(i) {
			case 1:
				Nel = atoi(argv[i]);
			break;

			case 2:
				Nsp = atoi(argv[i]);
			break;

			case 3:
				spt = static_cast<SP>(atoi(argv[i]));
			break;

			default:
				std::cout << "NONE?\n";

		}
	}
	
	// Set advection speed and time-step
	double a = 1.0;
	double dt = 0.01;

	// Number of iterations
	int NMAX = 1000;


	/* (Uniform) Mesh generation */
	double *x = new double[Nel+1]();
	// x[0] = -PI;
	// x[Nel] = PI;
	x[0] = 0.;
	x[Nel] = 1.0;
	// double h = (2.*PI)/Nel;
	double h = 1./Nel;
	for (auto i = 1; i != Nel; ++i) {
		x[i] = x[i-1]+h;
	}

	/* Calculating solution points coordinates inside each element */
	// Local SPs
	double *loc = new double[Nsp]();
	switch (spt) {
		case GAUSS:
			JacPZ(0, 0, Nsp, loc);
		break;

		case LOBATTO:

		break;

		case RADAU:

		break;
	}

	// Mapping local SPs to global SPs within each element
	double *xg = new double[Nsp*Nel]();
	qsi2x(xg, loc, x, Nel, Nsp);

	/* Setting initial conditions on SPs */
	double *u = new double[Nsp*Nel]();
	for (int i = 0; i != Nsp*Nel; ++i)
		u[i] = std::exp(-40*std::pow(xg[i]-0.5, 2)); //u[i] = std::sin(xg[i]);

	printSol(xg, u, Nel*Nsp);

	double *f = new double [Nsp*Nel]();
	double **d = new double*[Nsp]();
		for (int i = 0; i != Nsp; ++i) d[i] = new double[Nsp]();
	double *fd = new double[Nsp*Nel]();
	double *fup = new double[Nel]();
	double **fjump = new double*[Nel]();
		for (int i = 0; i != Nel; ++i) fjump[i] = new double[2]();
	double *gL = new double [Nsp]();
	double *gR = new double [Nsp]();
	double *gdR = new double [Nsp*Nel]();
	double *gdL = new double [Nsp*Nel]();
	double *dF = new double [Nsp*Nel]();
	double uL, uR, fL, fR, aT;


	// Derivative matrix
	switch (spt) {
		case GAUSS:
		for (int i = 0; i != Nsp; ++i) {
			for (int j = 0; j != Nsp; ++j) {
				d[i][j] = dJacP(0, 0, Nsp, loc[i])/(dJacP(0, 0, Nsp, loc[j])*(loc[i] - loc[j]));
			}
		}
		for (int i = 0; i != Nsp; ++i) d[i][i] = loc[i]/(1.-std::pow(loc[i], 2));
			break;

		case LOBATTO:

		break;

		case RADAU:

		break;
	}

	/* ------- TIME MARCHING ------- */
	for (int n = 0; n != NMAX; ++n) {

		/* Calculate (discontinuous) flux and its derivative at SPs */
		for (int i = 0; i != Nsp*Nel; ++i)
			f[i] = a*u[i];


		// Derivative of disc. flux (local coordinates)
		for (int i = 0; i != Nel; ++i) {
			for (int p = 0; p != Nsp; ++p) {
				fd[i*Nsp + p] = 0;
				for (int j = 0; j != Nsp; ++j) {
					fd[i*Nsp + p] += d[p][j]*f[i*Nsp + j];
				}
				//fd[i*Nsp + p] *= 2./(x[i+1] - x[i]);
			}
		}

		/* Calculate upwind fluxes */

		// Interpolate u at SPs to find u(x) on adj. elements, then set uL and uR at interface
		for (int i = 0; i != Nel-1; ++i) {
			uR = 0;
			uL = 0;
			for (int j = 0; j != Nsp; ++j) {
				uL += u[i*Nsp + j]*Lagrange(x[i+1], xg, Nsp, i, j);
				uR += u[(i+1)*Nsp + j]*Lagrange(x[i+1], xg, Nsp, i+1, j);
			}
			if (std::fabs(uL-uR) < 0.001)
				aT = ((a*uR) - (a*uL))/(uR - uL);
			else
				aT = a;

			fup[i] = (aT >= 0) ? a*uL : a*uR;
		}
		// Right interface of last element
		uR = 0;
		uL = 0;
		for (int j = 0; j != Nsp; ++j) {
			uL += u[(Nel-1)*Nsp + j]*Lagrange(x[Nel], xg, Nsp, Nel-1, j);
			uR += u[0*Nsp + j]*Lagrange(x[Nel], xg, Nsp, 0, j);
		}
		if (std::fabs(uL-uR) < 0.001)
			aT = ((a*uR) - (a*uL))/(uR - uL);
		else
			aT = a;

		fup[Nel-1] = (aT >= 0) ? a*uL : a*uR;

		/* Calculate jumps at each interface */

		// Construct f(x) so we can calculate fL and fR
		for (int i = 0; i != Nel-1; ++i) {
			fR = 0;
			fL = 0;
			for (int j = 0; j != Nsp; ++j) {
				fL += f[i*Nsp + j]*LagrangeLoc(loc[0], xg, Nsp, i, j);
				fR += f[i*Nsp + j]*LagrangeLoc(loc[Nsp-1], xg, Nsp, i, j);
			}
			if (i == 0) {
				fjump[i][0] = fup[Nel-1] - fL;
			} else {
				fjump[i][0] = fup[i-1] - fL;
			}
			fjump[i][1] = fup[i] - fR;
		}
	
		// Jumps for last element
		fR = 0;
		fL = 0;
		for (int j = 0; j != Nsp; ++j) {
			fL += f[(Nel-1)*Nsp + j]*Lagrange(x[Nel-1], xg, Nsp, (Nel-1), j);
			fR += f[(Nel-1)*Nsp + j]*Lagrange(x[Nel], xg, Nsp, (Nel-1), j);
		}
		fjump[Nel-1][0] = fup[Nel-2] - fL;
		fjump[Nel-1][1] = fup[Nel-1] - fR;


		/* Calculage correction function, g */
		// g_DG:
		for (int i = 0; i != Nsp; ++i) {
			gL[i] = 0.5*std::pow(-1., Nsp)*(JacP(0, 0, Nsp, loc[i]) - JacP(0, 0, Nsp-1, loc[i]));
			gR[Nsp-1-i] = gL[i];
		}

		// Derivative of g
		for (int i = 0; i != Nel; ++i) {
			for (int p = 0; p != Nsp; ++p) {
				gdL[i*Nsp + p] = 0;
				gdR[i*Nsp + p] = 0;
				for (int j = 0; j != Nsp; ++j) {
					gdL[i*Nsp + p] = d[p][j]*gL[j];
					gdR[i*Nsp + p] = d[p][j]*gR[j];
				}
				// gdL[i*Nsp + p] *= 2./(x[i+1] - x[i]);
				// gdR[i*Nsp + p] *= 2./(x[i+1] - x[i]);
			}
		}

		/* Reconstructed (continuous) flux F */
		for (int i = 0; i != Nel; ++i) {
			for (int j = 0; j != Nsp; ++j) {
				dF[i*Nsp + j] = fd[i*Nsp +j] + fjump[i][0]*gdL[i*Nsp +j] + fjump[i][1]*gdR[i*Nsp + j];
				dF[i*Nsp + j] *= 2./(x[i+1] - x[i]);
			}
		}

		/* Advance solution in time */
		for (int i = 0; i != Nel*Nsp; ++i)
			u[i] = u[i] - dt*dF[i];

		printSol(xg, u, Nel*Nsp);

	}

	return 0;
}