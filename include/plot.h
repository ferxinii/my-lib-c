#ifndef PLOT_H
#define PLOT_H

#include <stdio.h>

FILE *popen_gnuplot(void);

FILE *popen_gnuplot_video(const char *video_output, int framerate);

void set_config(FILE *pipe, const char *config);

void start_plot(FILE *pipe, const char *file_output);

void end_plot(FILE *pipe);

void video_to_gif(const char *name_video, const char *name_gif);

void add_point(FILE *pipe, double x, double y, const char *config);

void add_array_points(FILE *pipe, double **points, int N, const char *config);

void add_arrows_from_array_points(FILE* pipe, double **points, int N, int spacing, int offset, const char *config);

void add_yline(FILE *pipe, double y, const char *config);

void add_line(FILE *pipe, double x0, double xf, double y0, double yf, const char *config);

void add_function(FILE *pipe, double x0, double xf, int N, double (*fun)(double), const char *config);

void plot_array_2D(double **array, int N, const char *title, char *xlabel, char *ylabel, const char *file_name, 
                   double *IC, double *plotDimensions_x0_xf_y0_yf, char *config, double *arrows_size_freq_offset);

#endif
