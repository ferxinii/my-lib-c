#ifndef MATRIX_H
#define MATRIX_H

double **malloc_matrix(int N1, int N2);

void free_matrix(double **array, int N1);

void print_matrix(double **matrix, int N_rows, int N_cols);

double determinant(double **matrix, int N);

double **new_matrix_removed_row_column(double **matrix, int N, int row, int col);

double **new_arr_selected_columns(double **array, int N_rows, int N_cols, int *id_keep, int N_keep);

double dot_product(double *v1, double *v2, double N);

void equal_vectors(double *in, double *out, int N);

void multiply_matrices_2(double **A, double **B, double **out);

void equal_matrices_2(double **in, double **out);

void subtract_matrices_2(double **A, double **B, double **out);

void identity_matrix_2(double **A);

void multiply_matrix_vector_2(double **mat, double *vec, double *out);

void inverse_matrix_2(double **in, double **out);

int vector_components_are_lt_abs(double *vec, int N, double value);

void subtract_vectors(double *in1, double *in2, double *out, int N);

void sum_vectors(double *in, double *value, double *out, int N);

void multiply_matrix_vector(double **matrix, double *vec, double *out, int N);

void change_sign_vector_inplace(double *inout, int N);

void inverse_matrix(double **in, double **out, int N);

void lu_decomposition(double **A, double **L, double **U, int N);

void forward_substitution(double **L, double *B, double *Y, int N);

void backward_substitution(double **U, double *Y, double *X, int N);

void normalize_vector(double *inout, int dim);

void find_eigenvector(double **mat, int dim, double shift, double *out, int max_it, double tol);

void matrix_multiplication(double **A, double **B, double **C, int m, int n, int p);

void transpose(double **A, double **At, int m, int n);

void least_squares(double **A, double *b, double *x, int m, int n);

#endif
