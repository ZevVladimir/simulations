#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <Navier_Stokes/linalg.h>
#include <Navier_Stokes/fdm.h>
#include <Navier_Stokes/fluiddyn.h>
#include <Navier_Stokes/poisson_slvr.h>
#include <Navier_Stokes/utils.h>
#include <Navier_Stokes/setups.h>

#define PI 3.14159

int main(int argc, char *argv[]){

    const char *outdir = "nav_stoke_output";
    mkdir(outdir, 0777);

    double Re = 1000;
    double Lx = 1;
    double Ly = 1;

    int nx = 64;
    int ny = 64;
    double dt = 0.005;
    double tf = 20.0;
    double max_courant = 1.0;
    int order = 6;
    
    //TODO understand this
    double beta = 0.5 * (2 / (1 + sin(PI / (nx + 1))) + 2 / (1 + sin(PI / (ny + 1)))); // SOR poisson parameter

    // The maximum number of iterations possible
    int  it_max = (int)((tf / dt) - 1);
    // Poisson equation criterion for convergence
    double poisson_tol = 1E-3;                                                         

    // Dirichlet Boundary Conditions
    // General initial velocity conditions
    double ui = 0.0; 
    double vi = 0.0;

    // Initial velocity conditions on the boundaries
    double u1 = 0.0; // right boundary condition
    double u2 = 0.0; // left boundary condition
    double u3 = 0.0; // bottom boundary condition
    double u4 = 1.0; // top boundary condition

    double v1 = 0.0;
    double v2 = 0.0;
    double v3 = 0.0;
    double v4 = 0.0;

    double dx = Lx / nx;
    double dy = Ly / ny;

    // Generates derivative operators
    matrix *d_x = fdm_1(nx, order, dx);
    matrix *d_y = fdm_1(ny, order, dy);
    matrix *d_x2 = fdm_2(nx, order, dx);
    matrix *d_y2 = fdm_2(ny, order, dy);
    
    // identity matrix
    matrix *Ix = iden_matrix(nx);
    matrix *Iy = iden_matrix(ny);

    // Using kronecker product turn the 1D derivative operators into 2D ones
    //TODO I think because they only operate in 1D not necessarily that are 1D matrices? 
    
    matrix *DX = kronecker_prod(Ix, d_x);   // kronecker product for x first derivative
    matrix *DY = kronecker_prod(d_y, Iy);   // kronecker product for y first derivative
    matrix *DX2 = kronecker_prod(Ix, d_x2); // kronecker product for x second derivative
    matrix *DY2 = kronecker_prod(d_y2, Iy); // kronecker product for y second derivative

    free_matrix(d_x);
    free_matrix(d_y);
    free_matrix(d_x2);
    free_matrix(d_y2);
    free_matrix(Ix);
    free_matrix(Iy);

    // Calculate the courant numbers for each dimension
    double C1 = u1 * dt / (dx);
    double C2 = u1 * dt / (dy);

    if ((C1 > max_courant) || (C2 > max_courant))
    {
        printf("Unstable Solution!\n");
        printf("r1: %lf\n", C1);
        printf("r2: %lf\n", C2);
        exit(1);
    }

    // Initialize velocity matrices on the inside 
    matrix *u = init_matrix(nx, ny);
    matrix *v = init_matrix(nx, ny);

    // vorticity matrix
    matrix *w = init_matrix(nx, ny);   

    // stream-function
    matrix *psi = init_matrix(nx, ny); 

    // taylor_green_vortex(u,v,Lx,Ly,dx,dy,0,1.0/Re);

    for (int i = 1; i < nx; i++){
        for (int j = 1; j < nx; j++){
            u->M[i][j] = ui;
            v->M[i][j] = vi;
        }
    }

    for (int t = 0; t <= it_max; t++){
        //Boundary conditions
        for (int i = 0; i < nx; i++){
            u->M[i][0] = u2; // Left
            v->M[i][0] = v2; // Left
            u->M[i][ny-1] = u1; // Right
            v->M[i][ny-1] = v1; // Right
        }

        for (int j = 0; j < ny; j++){
            u->M[0][j] = u4; // Bottom
            v->M[0][j] = v4; // Bottom
            u->M[nx-1][j] = u3; // Top
            v->M[nx-1][j] = v3; // Top
        }

        // Convert our 2D velocity arrays to be 1D since that is easier to work with and do the matrix math
        matrix *u0 = reshape_matrix(u,nx*ny,1);
        matrix *v0 = reshape_matrix(v,nx*ny,1);

        // Calculate derivatives needed to find vorticity

        // First calculate in our condensed matrix form
        matrix *dudy0 = matrix_multiplication(DY,u0);
        matrix *dvdx0 = matrix_multiplication(DX,v0);

        // Return to full matrices
        matrix *dudy = reshape_matrix(dudy0, nx, ny);
        matrix *dvdx = reshape_matrix(dvdx0, nx, ny);

        // No longer need temporary matrices
        free_matrix(u0);
        free_matrix(v0);
        free_matrix(dudy0);
        free_matrix(dvdx0);

        // Calculate vorticity [w = dvy/dx - dvx/dy in the z direction] only on the boundaries
        // Still setting up bndry conditions vorticity in the middle is calculated in the time loop
        for (int j = 0; j < ny; j++){
            w->M[0][j] = dvdx->M[0][j] - dudy->M[0][j];
            w->M[nx-1][j] = dvdx->M[nx-1][j] - dudy->M[nx-1][j];
        }
        for (int i = 0; i < nx; i++){
            w->M[i][0] = dvdx->M[i][0] - dudy->M[i][0];
            w->M[i][ny-1] = dvdx->M[i][ny-1] - dudy->M[i][ny-1];
        }

        // Output initial conditions
        if (t == 0){
            output_matrix(u, t, "u",outdir);
            output_matrix(v, t, "v",outdir);
            output_matrix(w,t,"omega",outdir);
        }

        // Now compute derivatives of vorticity needed for the time advancement
        matrix *w0 = reshape_matrix(w,nx*ny,1);
        matrix *dwdx0 = matrix_multiplication(DX,w0);
        matrix *dwdy0 = matrix_multiplication(DY,w0);

        matrix *dwdx = reshape_matrix(dwdx0,nx,ny);
        matrix *dwdy = reshape_matrix(dwdy0,nx,ny);

        matrix *d2wdx20 = matrix_multiplication(DX2,w0);
        matrix *d2wdy20 = matrix_multiplication(DY2,w0);

        matrix *d2wdx2 = reshape_matrix(d2wdx20,nx,ny);
        matrix *d2wdy2 = reshape_matrix(d2wdy20,nx,ny);

        free_matrix(dwdx0);
        free_matrix(dwdy0);
        free_matrix(d2wdx20);
        free_matrix(d2wdy20);

        voritcity_euler(w,u,v,dwdx,dwdy,d2wdx2,d2wdy2,Re,dt);

        // Reset psi
        free_matrix(psi);

        // Invert the sign so it correctly matches the poisson equation we are solving 
        invert_sign(w);

        // TODO add poisson Successive Over-Relaxation

        // psi = poisson(w,dx,dy,it_max,poisson_tol);

        psi = poisson_SOR(w,dx,dy,it_max,poisson_tol,beta);

        // Un-invert the sign for use afterwards
        invert_sign(w);

        // Convert psi to 1D array for multiplication with derivatives
        // derivatives of psi are the velocities

        // Computes velocities using the fact that the stream function is built for incompressible flows
        matrix *psi0 = reshape_matrix(psi, nx * ny, 1);
        matrix *dpsidx0 = matrix_multiplication(DX, psi0);
        matrix *dpsidy0 = matrix_multiplication(DY, psi0);

        matrix *dpsidx = reshape_matrix(dpsidx0, nx, ny);
        matrix *dpsidy = reshape_matrix(dpsidy0, nx, ny);

        free_matrix(psi0);
        free_matrix(dpsidx0);
        free_matrix(dpsidy0);

        free_matrix(u);
        free_matrix(v);
        u = init_matrix(nx, ny);
        v = init_matrix(nx, ny);
        matrix_copy(u, dpsidy);
        invert_sign(dpsidx);
        matrix_copy(v, dpsidx);

        // Check that the new values satisfy the continuity equation that \del * u is close to 0 or that the fluid is still incompressible
        u0 = reshape_matrix(u, nx*ny, 1);
        v0 = reshape_matrix(v, nx*ny, 1);

        matrix *dudx0 = matrix_multiplication(DX,u0);
        matrix *dvdy0 = matrix_multiplication(DY,v0);

        matrix *dudx = reshape_matrix(dudx0,nx,ny);
        matrix *dvdy = reshape_matrix(dvdy0,nx,ny);

        matrix *continuity_ret = continuity_check(dudx,dvdy);
        
        if ((t % 10 == 0) && (t != 0)) {
            output_matrix(u, t, "u",outdir);
            output_matrix(v, t, "v",outdir);
            output_matrix(w, t, "omega",outdir);
            printf("Iteration: %d | ", t);
            printf("Time: %lf | ", (double)t * dt);
            printf("Progress: %.2lf%%\n", (double)100 * t / it_max);
            printf("Continuity max: %E | ", max_val(continuity_ret));
            printf("Continuity min: %E\n", min_val(continuity_ret));
        }       

        free_matrix(u0);
        free_matrix(v0);
        free_matrix(dudx0);
        free_matrix(dvdy0);

        //TODO output values

        free_matrix(dwdx);
        free_matrix(dwdy);
        free_matrix(d2wdx2);
        free_matrix(d2wdy2);
        free_matrix(dpsidx);
        free_matrix(dpsidy);
        free_matrix(dudx);
        free_matrix(dudy);
        free_matrix(dvdx);
        free_matrix(dvdy);
        free_matrix(continuity_ret);
    }

    free_matrix(u);
    free_matrix(v);
    free_matrix(w);
    free_matrix(psi);
    free_matrix(DX);
    free_matrix(DY);
    free_matrix(DX2);
    free_matrix(DY2);

    printf("Simulation complete!\n");

    return 0;
}