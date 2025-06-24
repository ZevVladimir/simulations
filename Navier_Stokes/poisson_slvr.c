#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Navier_Stokes/linalg.h>

// Calculate the RMSE error between 2 matrices
double calc_error(matrix *u1, matrix *u2){
    if (u1->nrow != u2->nrow || u1->ncol != u2->ncol){
        printf("Error: matrix size mismatch\n");
        exit(1);
    }

    double e = 0;

    for (int i = 0; i < u1->nrow; i++){
        for (int j = 0; j < u1->ncol; j++){
            e += pow(u1->M[i][j] - u2->M[i][j],2);
        }
    }
    e = e/(u1->nrow*u1->ncol);
    return sqrt(e);
}

matrix *poisson(matrix *f, double dx, double dy, int itmax, double tol){
    double curr_error;

    matrix *u = init_matrix(f->nrow,f->ncol);
    matrix *u0 = init_matrix(f->nrow,f->ncol);

    for (int l = 0; l < itmax; l++){
        matrix_copy(u0,u);
        // Avoid edges 
        for (int i = 1; i < f->nrow-1; i++){
            for (int j = 1; j < f->ncol-1; j++){
                u->M[i][j] = ((dy * dy * (u->M[i + 1][j] + u->M[i - 1][j])) + (dx * dx * (u->M[i][j + 1] + u->M[i][j-1])) - (dx * dx * dy * dy * f->M[i][j])) / (2 * (dx * dx + dy * dy));
            }
        }
        // Determine current difference between the past two iterations
        curr_error = calc_error(u, u0);
        if (curr_error < tol)
        {
            printf("Poisson equation solved with %d iterations - root-sum-of-squares error: %E\n", l, curr_error);
            free_matrix(u0);
            return u;
        }
    }
    printf("Error: maximum number of iterations achieved for Poisson equation.\n");

    free_matrix(u);
    free_matrix(u0);
    exit(1);
}

//TODO understand this
matrix *poisson_SOR(matrix *f, double dx, double dy, int itmax, double tol, double beta)
{
    int i, j, k, nx, ny;
    double curr_error;
    //double u1max, u2max;
    matrix *u = init_matrix(f->nrow, f->ncol);
    matrix *u0 = init_matrix(f->nrow, f->ncol);
    nx = f->nrow;
    ny = f->ncol;

    for (k = 0; k < itmax; k++)
    {
        matrix_copy(u0, u);
        for (i = 1; i < nx - 1; i++)
        {
            for (j = 1; j < ny - 1; j++)
            {
                u->M[i][j] = beta * (dy * dy * (u->M[i + 1][j] + u->M[i - 1][j]) + dx * dx * (u->M[i][j + 1] + u->M[i][j - 1]) - dx * dx * dy * dy * f->M[i][j]) / (2 * (dx * dx + dy * dy)) + (1 - beta) * u0->M[i][j];
            }
        }
        curr_error = calc_error(u, u0);
        if (curr_error < tol)
        {
            printf("Poisson equation solved with %d iterations - root-sum-of-squares error: %E\n", k, curr_error);
            free_matrix(u0);
            return u;
        }
    }
    printf("Error: maximum number of iterations achieved for Poisson equation.\n");

    free_matrix(u);
    free_matrix(u0);
    exit(1);
}
