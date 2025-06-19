#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>  // For creating directories (Linux/Mac)
#include <sys/types.h> // For directory types

#define NX 100 // number of x cells
#define NY 100 // number of y cells
#define xl 0.0 // left side
#define xr 10.0 // right side
#define yl 0.0 // left side
#define yr 10.0 // right side
#define ti 0.0
#define tf 10.0
#define VX 0.5 // advection velocity
#define VY 0.5
#define CFL 0.9 // Courant number
#define NU 0.05 // diffusion constant

// Use this function to return where in the 1D array the corresponding 2D array cell is
#define MOD(a, b) (((a) + (b)) % (b))  // Ensures positive result even for a < 0
#define IDX(i, j) ((MOD(i, NX)) * NY + MOD(j, NY))


void initialize(double *u, double dx, double dy){
    for(int i = 0; i < NX-1; i++){
        for(int j = 0; j < NY-1; j++){
            u[IDX(i,j)] = exp(-pow((dx*i + dy*j) - 5.0, 2) / 0.5);
        }
    }
}

void update_adv(double *u_old, double *u_new, double dx, double dy, double dt) {
    // Apply upwind scheme
    for (int i = 0; i < NX; i++) {  // Skip i=0 (left boundary)
        for (int j = 0; j < NY; j++){
            u_new[IDX(i,j)] = u_old[IDX(i,j)]
                - dt * VX * (u_old[IDX(i,j)] - u_old[IDX(MOD(i-1,NX),j)]) / dx
                - dt * VY * (u_old[IDX(i,j)] - u_old[IDX(i,MOD(j-1,NY))]) / dy;
        }
    }
}

void update_diff(double *u_old, double *u_new, double dx, double dy, double dt){
    
}

void export_to_csv(double *u, int step) {
    char directory_name[] = "data_adv_2d"; // Subdirectory name
    char filename[1000];             // File path for the CSV file

    // Create the subdirectory if it doesn't exist
    struct stat st = {0};
    if (stat(directory_name, &st) == -1) {
        mkdir(directory_name, 0700); // Create directory with read/write/execute permissions
    }

    // Construct the full file path in the subdirectory
    sprintf(filename, "%s/adv2d_step_%d.csv", directory_name, step);

    // Open the file for writing
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        printf("Error: Could not write to file %s\n", filename);
        return;
    }

    // Write the CSV header
    // fprintf(file, "x\n");

    // Write particle data
    for (int i = 0; i < NX*NY; i++) {
        fprintf(file, "%f\n", u[i]);
    }

    fclose(file);
    // printf("Exported data for step %d to %s\n", step, filename);
}

int main(){
    double *u = malloc(NX * NY * sizeof(double));
    double *u_new = malloc(NX * NY * sizeof(double));

    double dx = (xr - xl) / (NX);
    double dy = (yr - yl) / (NY);
    // double dt = fmin(0.5 * dx / V, 0.5 * dx * dx / NU); // Stable time step 
    double dt = CFL * dx / VX;
    printf("%f",dt);
    dt = 0.01;
    int n_step = 500;

    initialize(u,dx,dy);

    for(int i = 0; i < n_step; i++){
        // calc_flux(u,F);
        update_adv(u, u_new, dx, dy, dt);
        double *tmp = u;
        u = u_new;
        u_new = tmp;
        if(i % 5 == 0 || i == 0){
            export_to_csv(u,i);
        }
    }

    free(u);
    free(u_new);
    return 0;
}