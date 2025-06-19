#include <stdio.h>
#include <stdlib.h>
#include <Navier_Stokes/fluiddyn.h>
#include <Navier_Stokes/linalg.h>

void voritcity_euler(matrix w, matrix u, matrix v, matrix dwdx, matrix dwdy, matrix d2wdx2, matrix d2wdy2, double Re, double dt){
    int i, j;

    // General 3D equation https://uwaterloo.ca/applied-mathematics/current-undergraduates/continuum-and-fluid-mechanics-students/amath-463/vorticity
    // Since w is only in the z direction (w = del x \vec{u}) then first term is 0 as velocity is (u,v,0)
    // We assume incompressible flow such that \Del dot u = 0
    // We assume constant density so rho term goes to 0
    // nu is equivalent to 1/Re
    // Dw/Dt = material derivative or the chain rule applied of d/dt (w(x(t),y(t),t))
    // We calculate dw using this multiply dt to get how much it changs and add it to the original value

    for (i = 0; i < w.nrow; i++)
    {
        for (j = 0; j < w.ncol; j++)
        {
            w.M[i][j] = (-u.M[i][j] * dwdx.M[i][j] - v.M[i][j] * dwdy.M[i][j] + (1. / Re) * (d2wdx2.M[i][j] + d2wdy2.M[i][j])) * dt + w.M[i][j];
        }
    }
}

matrix continuity_check(matrix dudx, matrix dvdy){
    matrix cont_mtrx = init_matrix(dudx.nrow, dudx.ncol);

    for (int i = 0; i < dudx.nrow; i++){
        for (int j = 0; j < dudx.ncol; j++){
            cont_mtrx.M[i][j] = dudx.M[i][j] + dvdy.M[i][j];
        }
    }

    return cont_mtrx;
}
