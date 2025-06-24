#include <Navier_Stokes/linalg.h>

matrix *poisson(matrix *f, double dx, double dy, int itmax, double tol);

matrix *poisson_SOR(matrix *f, double dx, double dy, int itmax, double tol, double beta);
