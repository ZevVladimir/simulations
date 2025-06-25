#include <stdio.h>
#include <Navier_Stokes/linalg.h>

void output_matrix(matrix *A, int step, const char *name, const char *directory) {
    char filename[200];
    sprintf(filename, "%s/%s_step_%04d.csv", directory, name, step);

    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
        printf("Error: cannot open file %s for writing\n", filename);
        return;
    }

    for (int i = 0; i < A->nrow; i++) {
        for (int j = 0; j < A->ncol; j++) {
            fprintf(fp, "%.6f", A->M[i][j]);
            if (j < A->ncol - 1) {
                fprintf(fp, ",");
            }
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
}
