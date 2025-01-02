#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <dirent.h>
#include <string.h>
#include <sys/stat.h>  // For creating directories (Linux/Mac)
#include <sys/types.h> // For directory types
#include <string.h>
#include "structs.h"       // For Particle structure and constants
#include "integrators.h"    // For rk4 integrator
#include "barnes_hut.h"

// Define the Plummer scale radius as a constant
#define SCALE_RADIUS 20.0 // Plummer scale radius
#define M0 1e4 //mass of BH 
#define EPSILON 0.1
#define G 1
#define BOX_SIZE 5000

// Compute derivatives (gravitational forces) for the n-body system
void derivs(double t, double y[], double dydx[], int N) {
    // Zero out forces and velocities
    for (int i = 0; i < N; i++) {
        dydx[4 * i]     = y[4 * i + 2]; // dx/dt = vx
        dydx[4 * i + 1] = y[4 * i + 3]; // dy/dt = vy
        dydx[4 * i + 2] = 0.0;          // dvx/dt = Fx / m
        dydx[4 * i + 3] = 0.0;          // dvy/dt = Fy / m
    }
    // Compute pairwise forces
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i != j) {
                double dx = y[4 * j] - y[4 * i];         // xj - xi
                double dy = y[4 * j + 1] - y[4 * i + 1]; // yj - yi
                double r2 = dx * dx + dy * dy + EPSILON * EPSILON; // Softened distance squared
                double r = sqrt(r2);
                double F;

                F = G / r2; 

                // Force components
                double Fx = F * dx / r;
                double Fy = F * dy / r;

                // Update accelerations (action-reaction pairs)
                dydx[4 * i + 2] += Fx;  // dvx/dt for particle i
                dydx[4 * i + 3] += Fy;  // dvy/dt for particle i
            }
        }
        
    }
}

// Function to compute circular velocity for a Plummer sphere
double compute_circular_velocity(double r, double total_mass, double scale_radius) {
    double enclosed_mass = total_mass * pow(r, 3) / pow(r * r + scale_radius * scale_radius, 1.5);
    return sqrt(GRAVITATIONAL_CONSTANT * enclosed_mass / r);
}

// Assign proper velocities for circular orbits around the galaxy COM
void assign_orbital_velocities(Particle *particles, int N, double total_mass, double scale_radius, double ofst_x, double ofst_y, double ofst_vx, double ofst_vy) {
    //Black hole
    particles[0].x  = ofst_x;
    particles[0].y  = ofst_y;
    particles[0].vx  = ofst_vx;
    particles[0].vy  = ofst_vy;
    for (int i = 1; i < N; i++) {
        // Compute distance from the origin
        double r = sqrt(particles[i].x * particles[i].x + particles[i].y * particles[i].y);
        if (r == 0.0) continue; // Avoid division by zero

        // Compute circular velocity
        double v_circular = compute_circular_velocity(r, total_mass, scale_radius);

        // Compute tangential velocity components
        double theta = atan2(particles[i].y, particles[i].x); // Angle in polar coordinates
        particles[i].vx = -v_circular * sin(theta) + ofst_vx; // Perpendicular to radius
        particles[i].vy = v_circular * cos(theta) + ofst_vy;
    
        particles[i].x += ofst_x;
        particles[i].y += ofst_y;
    }
}

// Generate a Plummer sphere in 2D
void generate_plummer_2d(Particle *particles, int N, double total_mass) {
    double scale_mass = total_mass / N; // Equal mass per particle

    //Black hole
    particles[0].x  = 0.0;
    particles[0].y  = 0.0;
    particles[0].vx  = 0.0;
    particles[0].vy  = 0.0;
    particles[0].mass = M0;
    
    for (int i = 1; i < N; i++) {
        // Sample radial distance using the Plummer profile cumulative distribution
        double r = SCALE_RADIUS / sqrt(pow((double)rand() / RAND_MAX, -2.0 / 3.0) - 1);

        // Sample angular position randomly
        double theta = 2 * M_PI * ((double)rand() / RAND_MAX);

        // Convert polar to Cartesian coordinates
        particles[i].x = r * cos(theta);
        particles[i].y = r * sin(theta);

        // Assign particle mass
        particles[i].mass = scale_mass;

        // Initialize velocities as zero (to be assigned later)
        particles[i].vx = 0.0;
        particles[i].vy = 0.0;
    }
}

void export_to_csv(Particle *particles, int N, int step) {
    char directory_name[] = "data"; // Subdirectory name
    char filename[1000];             // File path for the CSV file

    // Create the subdirectory if it doesn't exist
    struct stat st = {0};
    if (stat(directory_name, &st) == -1) {
        mkdir(directory_name, 0700); // Create directory with read/write/execute permissions
    }

    // Construct the full file path in the subdirectory
    sprintf(filename, "%s/galaxy_step_%d.csv", directory_name, step);

    // Open the file for writing
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        printf("Error: Could not write to file %s\n", filename);
        return;
    }

    // Write the CSV header
    fprintf(file, "mass,x,y,vx,vy\n");

    // Write particle data
    for (int i = 0; i < N; i++) {
        fprintf(file, "%f,%f,%f,%f,%f\n", 
                particles[i].mass, // mass
                particles[i].x,    // x position
                particles[i].y,    // y position
                particles[i].vx,   // x velocity
                particles[i].vy    // y velocity
        );
    }

    fclose(file);
    printf("Exported data for step %d to %s\n", step, filename);
}

int load_from_csv(Particle *particles, int N, const char *directory_name, int *latest_step) {
    DIR *dir;
    struct dirent *entry;
    int max_step = -1;
    char latest_file[1000] = {0};

    // Open the directory
    if ((dir = opendir(directory_name)) == NULL) {
        printf("Error: Could not open directory %s\n", directory_name);
        return -1;
    }

    // Find the file with the largest step number
    while ((entry = readdir(dir)) != NULL) {
        if (strstr(entry->d_name, "galaxy_step_") != NULL) {
            int step;
            // Extract the step number from the filename (e.g., galaxy_step_123.csv)
            if (sscanf(entry->d_name, "galaxy_step_%d.csv", &step) == 1) {
                if (step > max_step) {
                    max_step = step;
                    // Construct the full path to the file
                    sprintf(latest_file, "%s/%s", directory_name, entry->d_name);
                }
            }
        }
    }
    closedir(dir);

    // Check if no valid snapshot files were found
    if (max_step == -1) {
        printf("Error: No snapshot files found in %s\n", directory_name);
        return -1;
    }

    // Load data from the file with the largest step number
    FILE *file = fopen(latest_file, "r");
    if (file == NULL) {
        printf("Error: Could not open file %s\n", latest_file);
        return -1;
    }

    printf("Resuming simulation from %s (step %d)\n", latest_file, max_step);

    // Store the largest step number found
    *latest_step = max_step;

    // Skip the header
    char line[256];
    fgets(line, sizeof(line), file);

    // Read particle data
    for (int i = 0; i < N; i++) {
        if (fgets(line, sizeof(line), file) != NULL) {
            sscanf(line, "%lf,%lf,%lf,%lf,%lf", 
                   &particles[i].mass, // mass
                   &particles[i].x,    // x position
                   &particles[i].y,    // y position
                   &particles[i].vx,   // x velocity
                   &particles[i].vy    // y velocity
            );
        }
    }
    fclose(file);
    return 0;
}


int main(int argc,char *argv[]) {
    int resume_simulation = 0;
    int latest_step = 0;  // To store the step number of the most recent snapshot

    // if (argc != 14) {
	// printf("Usage: %s Number of particles (galaxy 1) %s Number of particles (galaxy 2) %s Number of snapshots %s Time step %s Step output %s \nGlxy 1 x %s Glxy 1 y %s Glxy 1 vx %s Glxy 1 vy %s \nGlxy 2 x %s Glxy 2 y %s Glxy 2 vx %s Glxy 2 vy %s\n",
	// argv[0],argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],argv[9],argv[10],argv[11],argv[12],argv[13]);
	// return 1;
	// }
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--resume") == 0) {
            resume_simulation = 1;
            break;
        }
    }

	printf("%s Number of particles (galaxy 1) %s Number of particles (galaxy 2) %s Number of snapshots %s Time step %s Step output %s \nGlxy 1: x=%s, y=%s, vx=%s, vy=%s\nGlxy 2: x=%s, y=%s, vx=%s, vy=%s\n",
	argv[0],argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],argv[9],argv[10],argv[11],argv[12],argv[13]);

    // Simulation parameters
    int glxy1_N = atoi(argv[1]) + 1;  // Number of particles galaxy 1 + BH
    int glxy2_N = atoi(argv[2]) + 1;  // Number of particles galaxy 2 + BH
    double mass = 1.0;                // Mass of each particle
    int steps = atoi(argv[3]);        // Number of simulation steps
    double dt = atof(argv[4]);        // Time step
    int tot_nptls = glxy1_N + glxy2_N;
    int step_output = atoi(argv[5]); // After x many time steps output the data
    Particle *combined_galaxy = malloc(tot_nptls * sizeof(Particle));

    if (resume_simulation) {
        // Resume simulation from the most recent snapshot
        if (load_from_csv(combined_galaxy, tot_nptls, "data", &latest_step) != 0) {
            printf("Error: Failed to load snapshot. Exiting.\n");
            return 1;
        }
    } else {
        // Initialize galaxies as before (your original code)
        Particle *galaxy1 = malloc(glxy1_N * sizeof(Particle));
        Particle *galaxy2 = malloc(glxy2_N * sizeof(Particle));

        generate_plummer_2d(galaxy1, glxy1_N, mass * glxy1_N);
        generate_plummer_2d(galaxy2, glxy2_N, mass * glxy2_N);

        // Assign orbital velocities and combine into the state vector
        double glxy1_x = atof(argv[6]), glxy1_y = atof(argv[7]);
        double glxy1_vx = atof(argv[8]), glxy1_vy = atof(argv[9]);
        double glxy2_x = atof(argv[10]), glxy2_y = atof(argv[11]);
        double glxy2_vx = atof(argv[12]), glxy2_vy = atof(argv[13]);

        assign_orbital_velocities(galaxy1, glxy1_N, mass * glxy1_N + M0, SCALE_RADIUS, glxy1_x, glxy1_y, glxy1_vx, glxy1_vy);
        assign_orbital_velocities(galaxy2, glxy2_N, mass * glxy2_N + M0, SCALE_RADIUS, glxy2_x, glxy2_y, glxy2_vx, glxy2_vy);

        
        if (combined_galaxy == NULL) {
            fprintf(stderr, "Memory allocation for galaxy combination failed.\n");
            exit(EXIT_FAILURE);
        }

        memcpy(combined_galaxy, galaxy1, glxy1_N * sizeof(Particle));
        memcpy(combined_galaxy + glxy1_N, galaxy2, glxy2_N * sizeof(Particle));

        free(galaxy1);
        free(galaxy2);
    }

    // Simulation loop
    for (int step = latest_step; step < (steps+latest_step); step++) {
        Node *tree = build_tree(combined_galaxy,tot_nptls,-BOX_SIZE,BOX_SIZE,-BOX_SIZE,BOX_SIZE);
        if (step == 0){
        visualize_tree(tree,"tree_struct.txt");
        }
        // Export data every 10 steps and at the beginning 
        if (step % step_output == 0 | step == 0) {
            export_to_csv(combined_galaxy, tot_nptls, step);
        }

        lf2_BH(tree,combined_galaxy,tot_nptls,0.5,EPSILON,G,dt);

        // Compute derivatives
        // derivs(0.0, y, dydx, tot_nptls, glxy1_N);
        // // Integrate using Runge-Kutta
        // double *yout = malloc(4 * tot_nptls * sizeof(double)); // Temporary storage for updated state
        // rk4(y, dydx, 4 * tot_nptls, 0.0, dt, yout, tot_nptls, glxy1_N, derivs);

        // Update state vector
        // memcpy(y, yout, 4 * tot_nptls * sizeof(double));
        // free(yout);
        free_tree(tree);
    }

    // Free memory
    free(combined_galaxy);

    return 0;
}
