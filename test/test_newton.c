#include <stdio.h>
#include <math.h>
#include "nmds.h"


// Test function definitions
void fun(double *in, double *out) {
    out[0] = in[0]*in[0] + in[1]*in[1] - 1;  // x1^2 + x2^2 - 1
    out[1] = in[0] - in[1];  // x1 - x2
}

void jacobian(double *in, double **out) {
    out[0][0] = 2 * in[0];  // df1/dx1 = 2*x1
    out[0][1] = 2 * in[1];  // df1/dx2 = 2*x2
    out[1][0] = 1;          // df2/dx1 = 1
    out[1][1] = -1;         // df2/dx2 = -1
}

double fun_1D(double x) {
    return x * x - 2;  // f(x) = x^2 - 2
}

double dfun_1D(double x) {
    return 2 * x;  // f'(x) = 2x
}


int main() {
    int dim = 2;
    double x0[2] = {0.5, 0.5};  // Initial guess
    double x_out[2];  // Solution to be found

    // Perform the Newton method
    newton_method_vectorial(fun, jacobian, dim, x0, x_out, 1e-15, 100);
    double root = newton_method_1D(fun_1D, dfun_1D, 1, 1e-15, 100);

    // Output result
    printf("Solution: x1 = %.6f, x2 = %.6f\n", x_out[0], x_out[1]);
    printf("\n\n1D:\n");
    printf("Root found: x = %.10f\n", root);
    printf("Expected root: x = %.10f\n", sqrt(2));
    
    return 0;
}
