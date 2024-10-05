#ifndef MATRIX_H
#define MATRIX_H

double **malloc_matrix(int N1, int N2);

void free_matrix(double **array, int N1);

void print_matrix(double **matrix, int N_rows, int N_cols);

double determinant(double **matrix, int N);

double **new_matrix_removed_row_column(double **matrix, int N, int row, int col);

double dot_product(double *v1, double *v2, double N);

void multiply_matrices_2(double **A, double **B, double **out);

void equal_matrices_2(double **in, double **out);

void subtract_matrices_2(double **A, double **B, double **out);

void identity_matrix_2(double **A);

void multiply_matrix_vector_2(double **mat, double *vec, double *out);

void inverse_matrix_2(double **in, double **out);

#endif
