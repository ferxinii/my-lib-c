#ifndef MATRIX_H
#define MATRIX_H

double **malloc_matrix(int N1, int N2);

void free_matrix(double **array, int N1);

void multiply_matrices_2(double **A, double **B, double **out);

void equal_matrices_2(double **in, double **out);

void subtract_matrices_2(double **A, double **B, double **out);

void identity_matrix_2(double **A);

void multiply_matrix_vector_2(double **mat, double *vec, double *out);

void inverse_matrix_2(double **in, double **out);

#endif
