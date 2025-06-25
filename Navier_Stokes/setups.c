#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Navier_Stokes/linalg.h>
#define PI 3.14159

void taylor_green_vortex(matrix *u, matrix *v, double Lx, double Ly, double dx, double dy, double t, double nu){
    double a = (2 * PI) / Lx;
    double b = (2 * PI) / Ly;

    double A = 1;
    double B = - (A * a) / b;

    double k_sqr = a * a + b *b;

    double f = exp(-nu * k_sqr * t);
    
    for (int i = 0; i < u->nrow; i++){
        double x = i * dx;
        for (int j = 0; j < u->ncol; j++){
            double y = j * dy;
            // printf("%f,%f,%f,%f,%f\n",cos(a * x),sin(b * y),sin(a * x),cos(b * y),f);
            u->M[i][j] = A * cos(a * x) * sin(b * y) * f;
            v->M[i][j] = B * sin(a * x) * cos(b * y) * f;
        }
    }
}