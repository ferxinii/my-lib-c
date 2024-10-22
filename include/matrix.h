#ifndef MATRIX_H
#define MATRIX_H

double **malloc_matrix(int N1, int N2);

void free_matrix(double **array, int N1);

void print_matrix(double **matrix, int N_rows, int N_cols);

double determinant(double **matrix, int N);

double **new_matrix_removed_row_column(double **matrix, int N, int row, int col);

double dot_product(double *v1, double *v2, double N);

void equal_vectors(double *in, double *out, int N);

void multiply_matrices_2(double **A, double **B, double **out);

void equal_matrices_2(double **in, double **out);

void subtract_matrices_2(double **A, double **B, double **out);

void identity_matrix_2(double **A);

void multiply_matrix_vector_2(double **mat, double *vec, double *out);

void inverse_matrix_2(double **in, double **out);

int vector_components_are_lt_abs(double *vec, int N, double value);

void sum_vectors(double *in, double *value, double *out, int N);

void multiply_matrix_vector(double **matrix, double *vec, double *out, int N);

void change_sign_vector_inplace(double *inout, int N);

void inverse_matrix(double **in, double **out, int N);

void lu_decomposition(double **A, double **L, double **U, int N);

void forward_substitution(double **L, double *B, double *Y, int N);

void backward_substitution(double **U, double *Y, double *X, int N);

#endif
