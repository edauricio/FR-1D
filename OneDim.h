#include <string>

enum SP {GAUSS, LOBATTO, RADAU};
void qsi2x(double*, double*, double*, int, int);
double Lagrange(double, double*, int, int, int);
void printSol(double*, double*, int, std::string = "result");