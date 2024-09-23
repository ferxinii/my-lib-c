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

void add_yline(FILE *pipe, double y, const char *config);

void add_line(FILE *pipe, double x0, double xf, double y0, double yf, const char *config);

void add_function(FILE *pipe, double x0, double xf, int N, double (*fun)(double), const char *config);

#endif