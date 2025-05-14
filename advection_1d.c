#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>  // For creating directories (Linux/Mac)
#include <sys/types.h> // For directory types

#define N 100 // number of cells
#define xl 0.0 // left side
#define xr 10.0 // right side
#define ti 0.0
#define tf 10.0
#define V 0.5 // advection velocity
#define CFL 0.9 // Courant number
#define NU 0.05 // diffusion constant


void initialize(double *u, double *x, double dx){
    for(int i = 0; i < N-1; i++){
        x[i] = xl + i * dx;
        u[i] = exp(-pow(x[i] - 5.0, 2) / 0.5);
    }
}

void calc_flux(double *u, double *F){
    for (int i = 0; i < N - 1; i++){
        F[i] = (V/2) * (u[i+1] + u[i]) - (abs(V)/2) * (u[i+1] - u[i]);
    }
}

void update_adv_flux(double *u, double *u_new, double *F, double dx, double dt){
    for (int i = 1; i < N - 1; i++){
        u_new[i] = u[i] - dt/dx * (F[i] - F[i-1]);
    }
    u_new[0] = u_new[1];
    u_new[N-1] = u_new[N-2];

}

void update_adv(double *u, double *u_new, double dx, double dt) {
    double lambda = V * dt / dx;

    // Apply upwind scheme
    if(V > 0){
        for (int i = 1; i < N; i++) {  // Skip i=0 (left boundary)
            u_new[i] = u[i] - lambda * (u[i] - u[i-1]);
        }
        u_new[0] = u[N-1];
    }
    else{
        for (int i = 0; i < N-1; i++){
            u_new[i] = u[i] - lambda * (u[i+1] - u[i]);
        }
        u_new[N-1] = u_new[0];
    }
}

void update(double *u, double *u_new, double dx, double dt){
    double lambda = V * dt / dx;
    if (lambda >= 1.0) {
        printf("Warning: CFL number is too large, possible instability!\n");
    }
    for (int i = 1; i < N - 1; i++){
        double advection = - V * (u[i] - u[i-1])/dx;
        double diffusion = NU * (u[i + 1] - 2 * u[i] + u[i - 1]) / (dx * dx);
        u_new[i] = u[i] + dt * (advection + diffusion);
    }
    //Periodic bndry conds
    u_new[0] = u_new[N-2];
    u_new[N-1] = u_new[1];

    //Reflective bndry conds
    // u[0] = u[1];      
    // u[N+1] = u[N];    
}

void export_to_csv(double *u, double *x, int step) {
    char directory_name[] = "data_adv_1d"; // Subdirectory name
    char filename[1000];             // File path for the CSV file

    // Create the subdirectory if it doesn't exist
    struct stat st = {0};
    if (stat(directory_name, &st) == -1) {
        mkdir(directory_name, 0700); // Create directory with read/write/execute permissions
    }

    // Construct the full file path in the subdirectory
    sprintf(filename, "%s/adv1d_step_%d.csv", directory_name, step);

    // Open the file for writing
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        printf("Error: Could not write to file %s\n", filename);
        return;
    }

    // Write the CSV header
    fprintf(file, "x,u\n");

    // Write particle data
    for (int i = 0; i < N; i++) {
        fprintf(file, "%f,%f\n", 
                x[i],u[i]
        );
    }

    fclose(file);
    // printf("Exported data for step %d to %s\n", step, filename);
}

int main(){
    double x[N], u[N], u_new[N], F[N-1];
    double dx = (xr - xl) / (N);
    // double dt = fmin(0.5 * dx / V, 0.5 * dx * dx / NU); // Stable time step 
    double dt = CFL * dx / V;
    int n_step = 500;

    initialize(u,x,dx);

    for(int i = 0; i < n_step; i++){
        // calc_flux(u,F);
        update_adv(u, u_new, dx, dt);
        for (int i = 0; i < N; i++) {
            u[i] = u_new[i];
        }
        if(i % 5 == 0 || i == 0){
            export_to_csv(u,x,i);
        }
    }

}