typedef struct matrix
{
    double **M;
    int nrow;
    int ncol;
} matrix;

void zero_matrix(matrix A);

double **alloc_matrix(int nrow, int ncol);

matrix init_matrix(int nrow, int ncol);

void free_matrix(matrix *A);

matrix kronecker_prod(matrix *A, matrix *B);

matrix iden_matrix(int n);

matrix matrix_transpose(matrix A);

matrix reshape_matrix(matrix A, int new_nrow, int new_ncol);

matrix matrix_multiplication(matrix A, matrix B);

void invert_sign(matrix A);

void matrix_copy(matrix A, matrix B);

double max_val(matrix A);

double min_val(matrix A);
