
typedef void* MY_JET;

typedef double MY_FLOAT;

double **integrate_orbit(int (*taylor_uniform_step__ODE_NAME__tag)(MY_FLOAT *, MY_FLOAT *, int, int, double, double, MY_FLOAT *, MY_FLOAT *, int *, MY_JET *, int),
                         double *x0, int dim, double abs_total_time, int N_steps_forward, int N_steps_backward);

void plot_orbits_2D(double ***orbits_xy, int N_orbits, int N_steps, const char *title, const char *file_name, 
                    int mark_IC, double *plotDimensions_x0_xf_y0_yf, int *arrows_freq_offset);

double newton_method_1D(double (*fun)(double), double (*dfun)(double), double x0, double eps, int itmax, int verbose);

void newton_method_vectorial(void (*fun)(double *in, double *out), void (*jacobian)(double *in, double **out), 
                               int dim, double *x0, double *x_out, double eps, int itmax, int verbose);

double poincare_t(int (*taylor_uniform_step__ODE_NAME__tag)(MY_FLOAT *, MY_FLOAT *, int, int, double, double, MY_FLOAT *, MY_FLOAT *, int *, MY_JET *, int),
        double (*g)(double *x), void (*grad_g)(double *x, double *out), 
        void (*vector_field)(double *x, double *out),
        int dim, double *x0, int dir, double tol, int itmax, int n_crossings);

