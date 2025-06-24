#include <Navier_Stokes/linalg.h>

void voritcity_euler(matrix *w, matrix *u, matrix *v, matrix *dwdx, matrix *dwdy, matrix *d2wdx2, matrix *d2wdy2, double Re, double dt);

matrix *continuity_check(matrix *dudx, matrix *dvdy);