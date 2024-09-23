#include <stdlib.h>
#include <stdio.h>
#include <math.h>


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

