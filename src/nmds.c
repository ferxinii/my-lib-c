#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
#include "plot.h"
#include "nmds.h"


double **integrate_orbit(int (*taylor_uniform_step__ODE_NAME__tag)(MY_FLOAT *, MY_FLOAT *, int, int, double, double, MY_FLOAT *, MY_FLOAT *, int *, MY_JET *, int),
double *x0, int dim, double abs_total_time, int N_steps_forward, int N_steps_backward)
{
    static int tag = 0; tag++;
    double **orbit = malloc_matrix(N_steps_forward + N_steps_backward, dim);
    if (!orbit) {
        printf("Error allocating array for orbit!\n");  exit(1);
    }
    double t = 0;
    double delta_t = abs_total_time / (N_steps_forward + N_steps_backward);  
    double log10abserr = log10(1e-16), log10relerr = log10(1e-16);
    double x[dim];
    equal_vectors(x0, x, dim);

    // Forward
    for (int ii=N_steps_backward-1; ii>=0; ii--) {
        int out = taylor_uniform_step__ODE_NAME__tag(&t, x, -1, 1, log10abserr, log10relerr, 
                                                   NULL, &delta_t, NULL, NULL, tag);
        if (out != 0) {
            printf("Error in Taylor!\n");
            exit(1);
        }
        equal_vectors(x, orbit[ii], dim);
    }
    // Backward
    t = 0; 
    tag++;
    equal_vectors(x0, x, dim);
    for (int ii=N_steps_backward; ii<N_steps_backward + N_steps_forward; ii++) {
        int out = taylor_uniform_step__ODE_NAME__tag(&t, x, 1, 1, log10abserr, log10relerr, 
                                                   NULL, &delta_t, NULL, NULL, tag);
        if (out != 0) {
            printf("Error in Taylor!\n");
            exit(1);
        }
        equal_vectors(x, orbit[ii], dim);
    }

    return orbit;
}


double ***sample_orbits(int (*taylor_uniform_step__ODE_NAME__tag)(MY_FLOAT *, MY_FLOAT *, int, int, double, double, MY_FLOAT *, MY_FLOAT *, int *, MY_JET *, int),
                        double **array_IC, int N_IC, int dim, double abs_total_time, int N_steps_forward, int N_steps_backward)
{
    double ***orbits = malloc(sizeof(double **) * N_IC);
    for (int ii=0; ii<N_IC; ii++) {
        orbits[ii] = integrate_orbit(taylor_uniform_step__ODE_NAME__tag, array_IC[ii], 
                                     dim, abs_total_time, N_steps_forward, N_steps_backward);
    }
    return orbits;
}


double **initial_conditions_circle(double x, double y, double r, int N)
{
    double dth = 2 * M_PI / (N - 1);
    double th = 0;
    double **out = malloc_matrix(N, 2);
    for (int ii=0; ii<N; ii++) {
        out[ii][0] = x + r * cos(th);
        out[ii][1] = y + r * sin(th);
        th += dth;
    }
    return out;
}


double **initial_conditions_grid(double xmin, double xmax, double ymin, double ymax, double Nx, double Ny)
{
    double **out = malloc_matrix(Nx * Ny, 2);
    int kk = 0;
    for (int ii=0; ii<Nx; ii++) {
        for (int jj=0; jj<Ny; jj++) {
            out[kk][0] = xmin + ii * (xmax - xmin)/(Nx - 1);
            out[kk][1] = ymin + jj * (ymax - ymin)/(Ny - 1);
            kk++;
        }
    }
    return out;
}


double **initial_conditions_line(double xmin, double xmax, double ymin, double ymax, double N)
{
    double **out = malloc_matrix(N, 2);
    for (int ii=0; ii<N; ii++) {
        out[ii][0] = xmin + ii * (xmax - xmin)/(N - 1);
        out[ii][1] = ymin + ii * (ymax - ymin)/(N - 1);
    }
    return out;
}


void plot_orbits_2D(double ***orbits_xy, int N_orbits, int N_steps, const char *title, const char *file_name, 
                    int mark_IC, double *plotDimensions_x0_xf_y0_yf, char *config, double *arrows_size_freq_offset)
{
    FILE *pipe = popen_gnuplot();
    char buffer[256];
    if (plotDimensions_x0_xf_y0_yf) {
        snprintf(buffer, 256, "set xrange[%f:%f]", plotDimensions_x0_xf_y0_yf[0], plotDimensions_x0_xf_y0_yf[1]);
        set_config(pipe, buffer);
        snprintf(buffer, 256, "set yrange[%f:%f]", plotDimensions_x0_xf_y0_yf[2], plotDimensions_x0_xf_y0_yf[3]);
        set_config(pipe, buffer);
    }
    set_config(pipe, "set xlabel 'x'");
    set_config(pipe, "set ylabel 'y'");
    set_config(pipe, "set grid");
    set_config(pipe, "set key box opaque");
    snprintf(buffer, 256, "set title \"%s\" font ',24'", title);
    set_config(pipe, buffer);

    start_plot(pipe, file_name);
    for (int ii=0; ii<N_orbits; ii++) {
        if (mark_IC == 1) {
            add_point(pipe, orbits_xy[ii][0][0], orbits_xy[ii][0][1], "ps 2 pt 7 lc 8 notitle");
        }
        add_array_points(pipe, orbits_xy[ii], N_steps, config);
        if (arrows_size_freq_offset) {
            snprintf(buffer, 256, "lc 2 lw 3 size %f,20 fixed notitle", arrows_size_freq_offset[0]);
            add_arrows_from_array_points(pipe, orbits_xy[ii], N_steps, arrows_size_freq_offset[1], ii*arrows_size_freq_offset[2],
                                         buffer);
        }
    }
    end_plot(pipe);
    pclose(pipe); 
}


double newton_method_1D(double (*fun)(double), double (*dfun)(double), double x0, double eps, int itmax, int verbose)
{
    int it = 0;
    double x, x_prev = x0;
    while (it < itmax) {
        double den = dfun(x_prev);
        if (fabs(den) < 1e-15) {
            printf("Error in newton_method_1D: Dividing by 0...\n");
            exit(1);
        }
        double delta = - fun(x_prev) / den;
        x = x_prev + delta;
        x_prev = x;
        it++;
        if (verbose == 1) printf("it: %d, delta: %f \n", it, delta);
        if (fabs(delta) < eps) break;
    }
    if (it >= itmax) {
        printf("Reached maximum iters in newton_method_1D!\n");
        exit(1);
    }
    return x;
}


void newton_method_vectorial(void (*fun)(double *in, double *out), void (*jacobian)(double *in, double **out), 
                               int dim, double *x0, double *x_out, double eps, int itmax, int verbose)
{
    double *x = malloc(sizeof(double) * dim);
    double *x_prev = malloc(sizeof(double) * dim);
    double *fx_prev = malloc(sizeof(double) * dim);
    double *delta = malloc(sizeof(double) * dim);
    double **J = malloc_matrix(dim, dim);
    double **J_inv = malloc_matrix(dim, dim);
    
    int it = 0;
    equal_vectors(x0, x_prev, dim);
    while (it < itmax) {
        fun(x_prev, fx_prev);

        jacobian(x_prev, J);
        inverse_matrix(J, J_inv, dim);
        
        multiply_matrix_vector(J_inv, fx_prev, delta, dim);  // delta = - J_inv * f(x)
        change_sign_vector_inplace(delta, dim);
        sum_vectors(x_prev, delta, x, dim);  // x = x_prev + delta
        equal_vectors(x, x_prev, dim);
        it++;
        if (verbose == 1) printf("it: %d, delta[0]: %f \n", it, delta[0]);
        if (vector_components_are_lt_abs(delta, dim, eps)) break;
    }
    if (it >= itmax) {
        printf("Reached maximum iters in newton_method_vectorial!\n");
        exit(1);
    }
    equal_vectors(x, x_out, dim);

    free(x); free(x_prev); free(fx_prev); free(delta);
    free_matrix(J, dim); free_matrix(J_inv, dim);
}


int sign(double a)
{
    if (a > 0) return 1;
    else if (a < 0) return -1;
    else {
        printf("Error! 0 does not have a sign...?\n");
        exit(1);
    }
}


double poincare_t(int (*taylor_uniform_step__ODE_NAME__tag)(MY_FLOAT *, MY_FLOAT *, int, int, double, double, MY_FLOAT *, MY_FLOAT *, int *, MY_JET *, int),
        double (*g)(double *x), void (*grad_g)(double *x, double *out), 
        void (*vector_field)(double *x, double *out),
        int dim, 
        double *x0, 
        int dir, 
        double tol,
        int itmax, 
        int n_crossings)
{
    // We consider the poincare section g(x,y) = 0
    // We want to find t such that g(phi(t,x0)) = 0
    // That is, the zero of the function G(t) = g(phi(t,x0))
    
    // First, find the initial guess taking into account n_crossings
    double t;
    double log10abserr = log10(1e-16), log10relerr = log10(1e-16);
    double delta_t_aux = 1e-1, delta_t;  // ATTENTION This is arbitrary...
    int counter_crossings;  
    double x[dim]; 
    double t_prev, x_prev[dim];
    int out;
    
    int max_tries = 5;
    for (int ii=0; ii<max_tries; ii++) {
        t = 0;
        counter_crossings = 0;
        equal_vectors(x0, x, dim);

        delta_t = delta_t_aux;
        while (counter_crossings != n_crossings) {
            t_prev = t;
            equal_vectors(x, x_prev, dim);
            out = taylor_uniform_step__ODE_NAME__tag(&t, x, dir, 1, log10abserr, log10relerr, 
                                                     NULL, &delta_t, NULL, NULL, 0);
            if (out != 0) {
                printf("Something went wrong in Taylor step...\n");
                exit(1);
            }
            if (g(x_prev) * g(x) < 0) counter_crossings++;
        }
        // initial guess is t_prev, and t_prev < tau < t    
        
        // Now, iterate using Newton's method:
        t = t_prev;
        equal_vectors(x_prev, x, dim);
        int it = 0;
        while (fabs(g(x)) > tol && it < itmax) {
            it++;
            double Dg[dim], vf[dim];
            double denom;    

            grad_g(x, Dg);
            vector_field(x, vf);
            denom = dot_product(Dg, vf, dim);  // != 0 by definition (g transversal)

            delta_t = - g(x) / denom;
            
            // Update position and time
            dir = sign(delta_t);
            out = taylor_uniform_step__ODE_NAME__tag(&t, x, dir, 1, log10abserr, log10relerr, 
                                                     NULL, &delta_t, NULL, NULL, it);
            if (out != 0) {
                printf("Something went wrong in Taylor step...\n");
                exit(1);
            }
        }

        if (it < itmax) {
            return t;
        }

        delta_t_aux *= 0.1;
    }


    printf("Try %d in poincaré_t, delta_t = %f...\n", max_tries, delta_t_aux);
    printf("Reached maximum iterations in Newton's...\n");
    exit(1);
}


double max_abs_diff_w_initial(double (*fun)(double *x), double **orbit, int N_steps)
{

    double val0 = fun(orbit[0]);

    double max = 0;
    for (int ii=0; ii<N_steps; ii++) {
        double val = fun(orbit[ii]);
        if (fabs(val - val0) > max) max = fabs(val - val0);
    }
    return max;
}
