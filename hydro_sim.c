#include <stdio.h>
#include <math.h>
#include <sys/stat.h>  // For creating directories (Linux/Mac)
#include <sys/types.h> // For directory types

#define NX 100    // Number of grid points
#define GAMMA 1.4 // Specific heat ratio
#define DT 0.001  // Time step
#define DX 0.01   // Grid spacing
#define T_MAX 0.001 // Simulation time

typedef struct {
    double rho; // Density
    double u;   // Velocity
    double P;   // Pressure
    double e;   // Energy
} State;

void grid_init(State grid[]) {
    double dx = 1.0 / NX;
    for (int i = 0; i < NX; i++) {
        grid[i].e = 1.0 + 2.0 * exp(- (((dx * i)-0.5) * ((dx * i)-0.5)) / 0.006666); // Gaussian energy distrubtion
        grid[i].u = 0.0;
        grid[i].rho = 1.0;
        grid[i].P = grid[i].rho * grid[i].e * (GAMMA - 1.0);
    }
}

void compute_flux(State grid[], double flux[][3]) {
    for (int i = 0; i < NX; i++) {
        double rho = grid[i].rho;
        double u = grid[i].u;
        double P = grid[i].P;
        double E = P / (GAMMA - 1.0) + 0.5 * rho * u * u;

        flux[i][0] = rho * u;         // Mass flux
        flux[i][1] = rho * u * u + P; // Momentum flux
        flux[i][2] = u * (E + P);     // Energy flux
    }
}

void update(State grid[], double flux[][3]) {
    State new_grid[NX];
    for (int i = 1; i < NX - 1; i++) {
        double dU[3];
        dU[0] = -(flux[i][0] - flux[i - 1][0]) * DT / DX;
        dU[1] = -(flux[i][1] - flux[i - 1][1]) * DT / DX;
        dU[2] = -(flux[i][2] - flux[i - 1][2]) * DT / DX;

        new_grid[i].rho = grid[i].rho + dU[0];
        new_grid[i].u = (grid[i].rho * grid[i].u + dU[1]) / new_grid[i].rho;
        new_grid[i].P = (grid[i].P / (GAMMA - 1.0) + dU[2] - 0.5 * dU[0] * new_grid[i].u * new_grid[i].u) * (GAMMA - 1.0);

        new_grid[i].u = grid[i].u - (DT/2) * ((new_grid[i+1].P - new_grid[i].P) / flux[i][0]);
        new_grid[i].x = grid[i].x + (DT/2) * grid[i].u;
        new_grid[i].e = grid[i].e - (DT/2) * grid[i].P * ((grid[i].u - grid[i-1].u) / flux[i][0]);
        
    }

    // Update the grid with new values
    for (int i = 1; i < NX - 1; i++) {
        grid[i] = new_grid[i];
    }
}

void save_to_file(State grid[], int time_step) {
    char filename[50];
    char directory_name[] = "hydro_data"; // Subdirectory name

    // Create the subdirectory if it doesn't exist
    struct stat st = {0};
    if (stat(directory_name, &st) == -1) {
        mkdir(directory_name, 0700); // Create directory with read/write/execute permissions
    }
    sprintf(filename, "%s/hydro_output_%04d.txt", directory_name,time_step); // Generate filenames like output_0000.txt
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }

    for (int i = 0; i < NX; i++) {
        double x = i * DX;
        fprintf(file, "%.5f %.5f %.5f %.5f %.5f\n", x, grid[i].e, grid[i].rho, grid[i].u, grid[i].P);
    }

    fclose(file);
}


void simulate() {
    State grid[NX];
    double flux[NX][3];
    grid_init(grid);

    int time_step = 0;
    for (double t = 0; t < T_MAX; t += DT) {
        save_to_file(grid, time_step);
        compute_flux(grid, flux);
        update(grid, flux);
        time_step++;
    }
}

int main() {
    simulate();
    return 0;
}
