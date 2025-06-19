#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Navier_Stokes/linalg.h>
#include <sys/stat.h>  // For creating directories (Linux/Mac)
#include <sys/types.h> // For directory types

#define NX 11 // number of x cells (includes 1 extra for wall cells)
#define NY 11 // number of y cells (includes 1 extra for wall cells)
#define xl 0.0 // left side
#define xr 10.0 // right side
#define yl 0.0 // left side
#define yr 10.0 // right side
#define ti 0.0
#define tf 6.0
#define alpha 19.0 // mm^2/s
#define AC_temp 20
#define AC_speed 2
#define cx  NX / 2.0
#define cy  NY / 2.0
#define strength  10.0

void initialize(matrix temp, double dx, double dy){
    for(int i = 0; i < temp.nrow; i++){
        for(int j = 0; j < temp.ncol; j++){
            temp.M[i][j] = 26.67;
        }
    }
    temp.M[1][1] = AC_temp;
    temp.M[1][2] = AC_temp;
}

void update_vfield(matrix vx_old, matrix vx_new, matrix vy_old, matrix vy_new, double dx, double dy, double dt) {
    for (int i = 0; i < vx_old.nrow; i++) {  
        for (int j = 0; j < vx_old.ncol; j++){
            double vx_fan = 0.0;
            double vy_fan = 0.0;     
            if ((i > (cx - 2) && i < (cx + 2)) && (j > (cy - 2) && j < (cy + 2))){
                double x = j;
                double y = i;
                vx_fan = -(y - cy);
                vy_fan = (x - cx);
                double norm = sqrt(vx_fan * vx_fan + vy_fan * vy_fan) + 1e-8;
                vx_fan = strength * vx_fan / norm;
                vy_fan = strength * vy_fan / norm;
            }
            
            double viscosity = 0.1;

            double adv_vx = 0.0;
            double adv_vy = 0.0;
            double lap_vx = 0.0;
            double lap_vy = 0.0;
            if (i > 0 && i < vx_old.nrow - 1 && j > 0 && j < vx_old.ncol - 1) {
                adv_vx = -vx_old.M[i][j] * (vx_old.M[i+1][j] - vx_old.M[i-1][j]) / (2.0 * dx);
                adv_vy = -vy_old.M[i][j] * (vy_old.M[i][j+1] - vy_old.M[i][j-1]) / (2.0 * dy);
                
                lap_vx = (vx_old.M[i+1][j] + vx_old.M[i-1][j] - 2 * vx_old.M[i][j]) / (dx*dx)
                        + (vx_old.M[i][j+1] + vx_old.M[i][j-1] - 2 * vx_old.M[i][j]) / (dy*dy);
                lap_vx += viscosity * lap_vx;

                lap_vy = (vy_old.M[i+1][j] + vy_old.M[i-1][j] - 2 * vy_old.M[i][j]) / (dx*dx)
                        + (vy_old.M[i][j+1] + vy_old.M[i][j-1] - 2 * vy_old.M[i][j]) / (dy*dy);
                lap_vy += viscosity * lap_vy;
            }

            vx_new.M[i][j] = vx_old.M[i][j] + dt * (vx_fan + adv_vx + lap_vx);
            vy_new.M[i][j] = vy_old.M[i][j] + dt * (vy_fan + adv_vy + lap_vx);

            if ((i == 1 && j == 1) || (i == 1 && j==2)){
                vy_new.M[i][j] = AC_speed;
            }
        }
    }
}

void update_temp(matrix temp_old, matrix temp_new, double dx, double dy, double dt) {
    for (int i = 0; i < temp_old.ncol; i++) {  
        for (int j = 0; j < temp_old.nrow; j++){
            if ((i == 1 && j == 1) || (i == 1 && j == 2)){
                temp_new.M[i][j] = AC_temp;
                continue;
            }
            // Boundary conditions keep the walls at the original temperature
            if (i == 0 || i == NX - 1 || j == 0 || j == NY - 1) {
                temp_new.M[i][j] = temp_old.M[i][j];
                continue;
            }
            // Calculate diffusion
            double diffus = alpha * (
                    ((temp_old.M[i+1][j] - 2 * temp_old.M[i][j] + temp_old.M[i-1][j])/(dx*dx)) +
                    ((temp_old.M[i][j+1] - 2 * temp_old.M[i][j] + temp_old.M[i][j-1])/(dy*dy))
            );
            
            // Calculate advection
            double x = j;
            double y = i;
            double vx = -(y - cy);
            double vy = (x - cx);
            double norm = sqrt(vx*vx + vy*vy) + 1e-8;
            vx = strength * vx / norm;
            vy = strength * vy / norm;

            double advec = -1 * (vx * (temp_old.M[i][j+1] - temp_old.M[i][j-1]) / (2.0 * dx)) + (vy * (temp_old.M[i+1][j] - temp_old.M[i-1][j]) / (2.0 * dy));
            
            temp_new.M[i][j] = temp_old.M[i][j] + dt * (diffus + advec);
            }

        }
    }

void export_vel_field(int step, matrix vx, matrix vy) {
    char filename[256];
    sprintf(filename, "data_room_2d/velocity_step_%d.csv", step);

    FILE *file = fopen(filename, "w");
    if (!file) {
        printf("Error opening velocity output file\n");
        return;
    }

    for (int i = 0; i < vx.nrow; i++) {
        for (int j = 0; j < vx.ncol; j++) {
            fprintf(file, "%f,%f", vx.M[i][j], vy.M[i][j]);
            if (j < vx.nrow - 1) fprintf(file, ",");
        }
        fprintf(file, "\n");
    }

    fclose(file);
}

void export_to_csv(matrix temp, int step) {
    char directory_name[] = "data_room_2d"; // Subdirectory name
    char filename[1000];             // File path for the CSV file

    // Create the subdirectory if it doesn't exist
    struct stat st = {0};
    if (stat(directory_name, &st) == -1) {
        mkdir(directory_name, 0700); // Create directory with read/write/execute permissions
    }

    // Construct the full file path in the subdirectory
    sprintf(filename, "%s/room_2d_step_%d.csv", directory_name, step);

    // Open the file for writing
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        printf("Error: Could not write to file %s\n", filename);
        return;
    }

    // Write particle data
    for (int i = 0; i < temp.nrow; i++) {
        for (int j = 0; j < temp.ncol; j++) {
            fprintf(file, "%f", temp.M[i][j]);
            if (j < NY - 1) fprintf(file, ",");
        }
    fprintf(file, "\n");
}

    fclose(file);
    // printf("Exported data for step %d to %s\n", step, filename);
}

int main(){
    matrix temp = init_matrix(NX, NY);
    matrix temp_new = init_matrix(NX, NY);
    matrix vx = init_matrix(NX, NY);
    matrix vx_new = init_matrix(NX, NY);
    matrix vy = init_matrix(NX, NY);
    matrix vy_new = init_matrix(NX, NY);

    double dx = (xr - xl) / (NX);
    double dy = (yr - yl) / (NY);
    double dt = 0.01;
    
    int n_step = (tf - ti) / dt;

    initialize(temp,dx,dy);

    for(int i = 0; i < n_step; i++){
        if(i == 0){
            export_to_csv(temp,i);
            export_vel_field(i,vx,vy);
        }
        update_vfield(vx,vx_new,vy,vy_new,dx,dy,dt);

        matrix tmp_vx = vx;
        vx = vx_new;
        vx_new = tmp_vx;

        matrix tmp_vy = vy;
        vy = vy_new;
        vy_new = tmp_vy;

        update_temp(temp, temp_new, dx, dy, dt);

        matrix tmp_temp = temp;
        temp = temp_new;
        temp_new = tmp_temp;

        if(i % 5 == 0 && i != 0){
            export_to_csv(temp,i);
            export_vel_field(i,vx,vy);
        }
    }

    free_matrix(&temp);
    free_matrix(&temp_new);
    free_matrix(&vx);
    free_matrix(&vx_new);
    free_matrix(&vy);
    free_matrix(&vy_new);

    return 0;
}