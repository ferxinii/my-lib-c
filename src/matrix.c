#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "matrix.h"


double **malloc_matrix(int N1, int N2)
{
    double **array = malloc(sizeof(double*) * N1);
    for (int ii=0; ii<N1; ii++) {
        array[ii] = malloc(sizeof(double) * N2);
    }
    return array;
}


void free_matrix(double **array, int N1)
{
    for (int ii=0; ii<N1; ii++) {
        free(array[ii]);
    }
    free(array);
}


void print_matrix(double **matrix, int N_rows, int N_cols)
{   
    for (int ii=0; ii<N_rows; ii++) {
        for (int jj=0; jj<N_cols; jj++) {
            printf("%f  ", matrix[ii][jj]);
        }
        printf("\n");
    }
}


double determinant(double **matrix, int N)
{   
    double det = 0;
    if (N == 2) {
        det = matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1];   
    } else {
        int sign = 1;
        for (int irow=0; irow<N; irow++) {
            double **reduced = new_matrix_removed_row_column(matrix, N, irow, 0);
            det += sign * matrix[irow][0] * determinant(reduced, N-1);
            sign *= -1;
            free_matrix(reduced, N-1);
        }
    }
    return det;   
}


double **new_matrix_removed_row_column(double **matrix, int N, int row, int col)
{
    int ii_eff = 0, jj_eff = 0;
    double **out = malloc_matrix(N-1, N-1);
    for (int ii=0; ii<N; ii++) {
        if (ii != row) {
            jj_eff = 0;
            for (int jj=0; jj<N; jj++) {
                if (jj != col) {
                    out[ii_eff][jj_eff] = matrix[ii][jj];
                    jj_eff++;
                }
            }
            ii_eff++;
        }
    }
    return out;
}


double **new_arr_selected_columns(double **array, int N_rows, int N_cols, int *id_keep, int N_keep)
{
    double **out = malloc_matrix(N_rows, N_keep);
    for (int kk=0; kk<N_keep; kk++) {
        for (int jj=0; jj<N_cols; jj++) {
            if (jj == id_keep[kk]) {
                for(int ii=0; ii<N_rows; ii++) {
                    out[ii][kk] = array[ii][jj];
                }
            }
        }
    }
    return out;
}


double dot_product(double *v1, double *v2, double N)
{
    double out = 0;
    for (int ii=0; ii<N; ii++) {
        out += v1[ii] * v2[ii];
    }
    return out;
}


void equal_vectors(double *in, double *out, int N)
{
    for (int ii=0; ii<N; ii++) {
        out[ii] = in[ii];
    }
}



void multiply_matrices_2(double **A, double **B, double **out)
{
    out[0][0] = A[0][0]*B[0][0] + A[0][1]*B[1][0];
    out[0][1] = A[0][0]*B[0][1] + A[0][1]*B[1][1];
    out[1][0] = A[1][0]*B[0][0] + A[1][1]*B[1][0];
    out[1][1] = A[1][0]*B[0][1] + A[1][1]*B[1][1];
}


void equal_matrices_2(double **in, double **out)
{
    out[0][0] = in[0][0];
    out[0][1] = in[0][1];
    out[1][0] = in[1][0];
    out[1][1] = in[1][1];
}


void subtract_matrices_2(double **A, double **B, double **out)
{
    out[0][0] = A[0][0] - B[0][0];
    out[0][1] = A[0][1] - B[0][1];
    out[1][0] = A[1][0] - B[1][0];
    out[1][1] = A[1][1] - B[1][1];
}


void identity_matrix_2(double **A)
{
    A[0][0] = 1;
    A[0][1] = 0;
    A[1][0] = 0;
    A[1][1] = 1;
}


void multiply_matrix_vector_2(double **A, double *x, double *out)
{
    out[0] = A[0][0] * x[0] + A[0][1] * x[1];
    out[1] = A[1][0] * x[0] + A[1][1] * x[1];
}


void inverse_matrix_2(double **in, double **out)
{
    double det;
    det = in[0][0] * in[1][1] - in[1][0] * in[0][1];
    
    if (fabs(det) < 1e-16 ) {
        printf("Dividing by a very small number!!\n");
    }
    if (det == 0) {
        printf("Actually, dividing by 0...\n STOPPING\n");
        exit(1);
    }

    out[0][0] = 1 / det * in[1][1];
    out[0][1] = - 1 / det * in[0][1];
    out[1][0] = - 1 / det * in[1][0];
    out[1][1] = 1 / det * in[1][1];
}


int vector_components_are_lt_abs(double *vec, int N, double value)
{
    int indicator = 1;
    for (int ii=0; ii<N; ii++) {
        if (fabs(vec[ii]) >= value) indicator = 0;
    }
    return indicator;
}


void sum_vectors(double *in, double *value, double *out, int N)
{
    for (int ii=0; ii<N; ii++) {
        out[ii] = in[ii] + value[ii];
    }
}


void subtract_vectors(double *in1, double *in2, double *out, int N)
{
    for (int ii=0; ii<N; ii++) {
        out[ii] = in1[ii] - in2[ii];
    }
}


void multiply_matrix_vector(double **matrix, double *vec, double *out, int N)
{
    for (int ii=0; ii<N; ii++) {
        out[ii] = 0;
        for (int jj=0; jj<N; jj++) {
            out[ii] += matrix[ii][jj] * vec[jj];
        }
    }
}


void change_sign_vector_inplace(double *inout, int N)
{
    for (int ii=0; ii<N; ii++) {
        inout[ii] = -1 * inout[ii];
    }
}


void inverse_matrix(double **in, double **out, int N)
{
    double **L = malloc_matrix(N, N);
    double **U = malloc_matrix(N, N);

    lu_decomposition(in, L, U, N);

    double *B = (double *)malloc(N * sizeof(double));
    double *Y = (double *)malloc(N * sizeof(double));
    double *X = (double *)malloc(N * sizeof(double));

    // Solve for each column of the inverse matrix
    for (int jj = 0; jj < N; jj++) {
        // Set B as the j-th column of the identity matrix
        for (int ii = 0; ii < N; ii++) {
            B[ii] = (ii == jj) ? 1.0 : 0.0;
        }
        forward_substitution(L, B, Y, N);  // Solve L * Y = B (forward substitution)
        backward_substitution(U, Y, X, N);  // Solve U * X = Y (backward substitution)

        // Set the j-th column of the inverse matrix
        for (int ii = 0; ii < N; ii++) {
            out[ii][jj] = X[ii];
        }
    }

    free(B); free(Y); free(X);
    free_matrix(L, N); free_matrix(U, N);
}



void lu_decomposition(double **A, double **L, double **U, int N)
{
    for (int i = 0; i < N; i++) {
        // Upper triangular matrix U
        for (int k = i; k < N; k++) {
            double sum = 0.0;
            for (int j = 0; j < i; j++) {
                sum += (L[i][j] * U[j][k]);
            }
            U[i][k] = A[i][k] - sum;
        }

        // Lower triangular matrix L
        for (int k = i; k < N; k++) {
            if (i == k) {
                L[i][i] = 1.0;  // Diagonal as 1
            } else {
                double sum = 0.0;
                for (int j = 0; j < i; j++) {
                    sum += (L[k][j] * U[j][i]);
                }
                L[k][i] = (A[k][i] - sum) / U[i][i];
            }
        }
    }
}


void forward_substitution(double **L, double *B, double *Y, int N)
{  // Function to solve a system of linear equations L * Y = B using forward substitution
    for (int i = 0; i < N; i++) {
        double sum = 0.0;
        for (int j = 0; j < i; j++) {
            sum += L[i][j] * Y[j];
        }
        Y[i] = (B[i] - sum);
    }
}


void backward_substitution(double **U, double *Y, double *X, int N)
{  // Function to solve a system of linear equations U * X = Y using backward substitution
    for (int i = N - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < N; j++) {
            sum += U[i][j] * X[j];
        }
        X[i] = (Y[i] - sum) / U[i][i];
    }
}


void normalize_vector(double *inout, int dim)
{
    double norm = 0;
    for (int ii=0; ii<dim; ii++) {
        norm += inout[ii] * inout[ii];
    }
    norm = sqrt(norm);
    
    for (int ii=0; ii<dim; ii++) {
        inout[ii] = inout[ii] / norm;
    }
}


void find_eigenvector(double **mat, int dim, double shift, double *out, int max_it, double tol)
{
    for (int ii=0; ii<dim; ii++) {  // Initial guess is vector (1, 1, ... 1)
        out[ii] = 1;  
    }
    normalize_vector(out, dim);
    
    double **aux1 = malloc_matrix(dim, dim);  //  A - s * Id
    for (int ii=0; ii<dim; ii++) {
        for (int jj=0; jj<dim; jj++) {
            if (ii == jj) {
                aux1[ii][jj] = mat[ii][jj] - shift;
            } else {
                aux1[ii][jj] = mat[ii][jj];
            }
        }
    }
    
    double **aux2 = malloc_matrix(dim, dim);  // inverse of aux1
    inverse_matrix(aux1, aux2, dim);
    
    double *prev_vec = malloc(sizeof(double) * dim);
    double *diff_vec = malloc(sizeof(double) * dim);
    equal_vectors(out, prev_vec, dim);
    int it = 0;
    for (it=0; it<max_it; it++) {
        multiply_matrix_vector(aux2, prev_vec, out, dim);
        normalize_vector(out, dim);

       // Ensure consistent direction w dot product w prev_vec
        double dot_product = 0;
        for (int jj = 0; jj < dim; jj++) {
            dot_product += out[jj] * prev_vec[jj];
        }
        if (dot_product < 0) {
            for (int jj = 0; jj < dim; jj++) {
                out[jj] = -out[jj];   
            }
        }

        subtract_vectors(out, prev_vec, diff_vec, dim);
        equal_vectors(out, prev_vec, dim);
        if (vector_components_are_lt_abs(diff_vec, dim, tol)) break;
    }
    
    if (it >= max_it - 1) {
        printf("Reached maximum iters approximating eigenvector...");
        exit(1);
    }
}


void matrix_multiplication(double **A, double **B, double **C, int m, int n, int p)
{  // m: rows of A, n: cols of A, p: cols of B
    for (int ii = 0; ii < m; ii++) {
        for (int jj = 0; jj < p; jj++) {
            C[ii][jj] = 0.0;
            for (int kk = 0; kk < n; kk++) {
                C[ii][jj] += A[ii][kk] * B[kk][jj];
            }
        }
    }
}


void transpose(double **A, double **At, int m, int n)
{   // m: rows of A, n: cols of A
    for (int ii = 0; ii < m; ii++) {
        for (int jj = 0; jj < n; jj++) {
            At[jj][ii] = A[ii][jj];
        }
    }
}


void least_squares(double **A, double *b, double *x, int m, int n)
{   // m: rows of A, n: columns of a
    double **At = malloc_matrix(m, n);
    double **AtA = malloc_matrix(n, n);
    double *Atb = malloc(sizeof(double) * n);

    transpose(A, At, m, n);
    // Compute A^T * A
    matrix_multiplication(At, A, AtA, n, m, n);
    // Compute A^T * b
    multiply_matrix_vector(At, b, Atb, n);

    // Solve (A^T * A) * x = A^T * b
    double **L = malloc_matrix(n, n);
    double **U = malloc_matrix(n, n);
    double *Y = malloc(sizeof(double) * n);
    lu_decomposition(AtA, L, U, n);
    forward_substitution(L, Atb, Y, n);
    backward_substitution(U, Y, x, n);

    // Free allocated memory
    free_matrix(At, m);
    free_matrix(AtA, n);
    free(Atb);
    free_matrix(L, n);
    free_matrix(U, n);
    free(Y);
}
