#include <stdio.h>
#include <stdlib.h>
#include <Navier_Stokes/linalg.h>

void zero_matrix(matrix A){
    for (int i = 0; i < A.nrow; i++){
        for (int j = 0; j < A.ncol; j++){
            A.M[i][j] = 0.0;
        }
    }
}

double **alloc_matrix(int nrow, int ncol){
    double **M;

    if ((nrow < 1) || (ncol < 1)){
        printf("Matrix must have at least 1 row and one column\n");
        exit(1);
    }

    M = (double **) malloc(nrow * sizeof(double *));

    if (M == NULL){
        printf("Insufficient memory to construct matrix\n");
        exit(1);
    }

    for (int i = 0; i < nrow; i++){
        M[i] = (double *)malloc(ncol * sizeof(double));

        if (M[i] == NULL){
            printf("Insufficient memory to construct matrix row\n");
            exit(1);
        }
    }

    return M;
}

matrix init_matrix(int nrow, int ncol){
    matrix A;
    A.M = alloc_matrix(nrow, ncol);
    A.nrow = nrow;
    A.ncol = ncol;
    zero_matrix(A);

    return A;
}

// We pass a pointer to the matrix instead of the matrix itself as we don't want to create a copy and leave the original still existing
void free_matrix(matrix *A){
    if (A->M == NULL)
        return;

    if ((A->nrow < 1) || (A->ncol < 1))
    {
        printf("Matrix must have at least 1 row and one column\n");
        exit(1);
    }

    for (int i = 0; i < A->nrow; i++) {
        free(A->M[i]);
    }

    free(A->M);
}

matrix kronecker_prod(matrix *A, matrix *B){
    int a = A->nrow;
    int b = A->ncol;
    int p = B->nrow;
    int q = B->ncol;

    matrix C;
    C.nrow = a * p;
    C.ncol = b * q;
    C.M = alloc_matrix(C.nrow, C.ncol);

    for (int i = 0; i < C.nrow; i++){
        for (int j = 0; j < C.ncol; j++){
            C.M[i][j] = A->M[i/p][j/q] * B->M[i%p][j%q];
        }
    }
    return C;
}

matrix iden_matrix(int n){
    matrix I = init_matrix(n,n);

    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (i == j){
                I.M[i][j] = 1.0;
            }
        }
    }

    return I;
}

matrix reshape_matrix(matrix A, int new_nrow, int new_ncol){
    if ((A.nrow * A.ncol) != (new_nrow * new_ncol))
    {
        printf("Error: the reshaped matrix must have the same number of elements\n");
        printf("Number of elements of input matrix: %d\n", A.nrow * A.ncol);
        printf("Number of elements of output matrix: %d\n", new_nrow * new_ncol);
        exit(1);
    }

    matrix B = init_matrix(new_nrow,new_ncol);

    // Indices for the inputted matrix
    int k,l;
    k = 0;
    l = 0;

    for (int i = 0; i < new_nrow; i++){
        for (int j = 0; j < new_ncol; j++){
            B.M[i][j] = A.M[k][l];
            //Also increment k and l according to how many columns there are in A
            if (l < A.ncol - 1){
                l++;
            }
            else{
                k++;
                l = 0;
            }
        }
    }

    return B;
}

matrix matrix_transpose(matrix A){
    matrix A_T = init_matrix(A.ncol,A.nrow);

    for (int i = 0; i < A.nrow; i++){
        for (int j = 0; j < A.ncol; j++){
            A_T.M[i][j] = A.M[j][i];
        }
    }

    return A_T;
}

matrix matrix_multiplication(matrix A, matrix B){
    if (A.ncol != B.nrow){
        printf("Error: the number of columns of A should be equal to the number of rows of B");
        exit(1);
    }
    
    // We can transpose B to take advantage of the fact that C is row major based
    matrix B_T = matrix_transpose(B);

    matrix C = init_matrix(A.nrow,B.ncol);

    for (int i = 0; i < A.nrow; i++){
        for (int j = 0; j < B.ncol; j++){
            for (int k = 0; k < A.ncol; k++){
                C.M[i][j] += A.M[i][k] * B_T.M[j][k];
            }
        }
    }

    return C;
}

void invert_sign(matrix A){
    for (int i = 0; i < A.nrow; i++){
        for (int j = 0; j < A.ncol; j++){
            A.M[i][j] = -A.M[i][j];
        }
    }
}

void matrix_copy(matrix A, matrix B){
    for (int i = 0; i < A.nrow; i++){
        for (int j = 0; j < A.ncol; j++){
            A.M[i][j] = B.M[i][j];
        }
    }
}

double max_val(matrix A){
    double max_val = -__DBL_MAX__;

    for (int i = 0; i < A.nrow; i++){
        for (int j = 0; j < A.ncol; j++){
            if (A.M[i][j] > max_val){
                max_val = A.M[i][j];
            }
        }
    }

    return max_val;
}

double min_val(matrix A){
    double min_val = __DBL_MAX__;

    for (int i = 0; i < A.nrow; i++){
        for (int j = 0; j < A.ncol; j++){
            if (A.M[i][j] < min_val){
                min_val = A.M[i][j];
            }
        }
    }

    return min_val;
}
