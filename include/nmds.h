
typedef void* MY_JET;

typedef double MY_FLOAT;

int sign(double a);

double **integrate_orbit(int (*taylor_uniform_step__ODE_NAME__tag)(MY_FLOAT *, MY_FLOAT *, int, int, double, double, MY_FLOAT *, MY_FLOAT *, int *, MY_JET *, int),
                         double *x0, int dim, double abs_total_time, int N_steps_forward, int N_steps_backward);

double ***sample_orbits(int (*taylor_uniform_step__ODE_NAME__tag)(MY_FLOAT *, MY_FLOAT *, int, int, double, double, MY_FLOAT *, MY_FLOAT *, int *, MY_JET *, int),
                        double **array_IC, int N_IC, int dim, double abs_total_time, int N_steps_forward, int N_steps_backward);

double **initial_conditions_circle(double x, double y, double r, int N);

double **initial_conditions_grid(double xmin, double xmax, double ymin, double ymax, double Nx, double Ny);

double **initial_conditions_line(double xmin, double xmax, double ymin, double ymax, double N);

void plot_orbits_2D(double ***orbits_xy, int N_orbits, int N_steps, const char *title, const char *file_name, 
                    double *IC, double *plotDimensions_x0_xf_y0_yf, char *config, double *arrows_size_freq_offset);

double bisection_method_1D(double (*fun)(double), double x1, double x2, double tol, int it_max);

double newton_method_1D(double (*fun)(double), double (*dfun)(double), double x0, double eps, int itmax, int verbose);

void newton_method_vectorial(void (*fun)(double *in, double *out), void (*jacobian)(double *in, double **out), 
                               int dim, double *x0, double *x_out, double eps, int itmax, int verbose);

double poincare_t(int (*taylor_uniform_step__ODE_NAME__tag)(MY_FLOAT *, MY_FLOAT *, int, int, double, double, MY_FLOAT *, MY_FLOAT *, int *, MY_JET *, int),
        double (*g)(double *x), void (*grad_g)(double *x, double *out), 
        void (*vector_field)(double *x, double *out),
        int dim, double *x0, int dir, double t_steps_0, double tol, int itmax, int n_crossings);

double max_abs_diff_w_initial(double (*fun)(double *x), double **orbit, int N_steps);

