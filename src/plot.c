#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>


FILE *popen_gnuplot(void)
{
    FILE *pipe = popen("gnuplot -persistent 2>&1", "w");
    if (!pipe) {
        printf("ERROR! Cannot open gnuplot pipe...\n");
        exit(1);
    }

    fprintf(pipe, "set terminal pngcairo enhanced font 'Arial,18' size 1080,1080 enhanced \n");

    return pipe;
}


FILE *popen_gnuplot_video(const char *video_output, int framerate)
{   
    char buffer[256];
    snprintf(buffer, 256, "rm -f %s", video_output);
    system(buffer);

    FILE *pipe = popen_gnuplot();
    fprintf(pipe, "set output '| ffmpeg -loglevel error -f png_pipe -r %d -s:v 1920x1080 -i pipe: -pix_fmt yuv420p -c:v libx264 -crf 18 %s'\n", framerate, video_output);
    
    return pipe;
}


void set_config(FILE *pipe, const char *config)
{
    fprintf(pipe, "%s\n", config);
}


void start_plot(FILE *pipe, const char *file_output)
{
    fprintf(pipe, "");
    if (file_output) {
        fprintf(pipe, "set output '%s'\n", file_output);
    }

    fflush(pipe);
    fprintf(pipe, "plot ");
}


void end_plot(FILE *pipe)
{   
    fflush(pipe);
    fprintf(pipe, "\n");
}


void video_to_gif(const char *name_video, const char *name_gif)
{
      sleep(1);  //Delay necessary to ensure mp4 exists?
      char command_mp4_to_gif[1024];
      char cwd[256];
      getcwd(cwd, 256);
      snprintf(command_mp4_to_gif, 1024, "ffmpeg -loglevel error -i %s/%s -vf \"fps=12,scale=640:-1:flags=lanczos\" -gifflags +transdiff -y %s 2>&1", cwd, name_video, name_gif);
      system(command_mp4_to_gif);
}


void add_point(FILE *pipe, double x, double y, const char *config)
{
    fprintf(pipe, "\"<echo \'%f %f\'\" w p %s, ", x, y, config);
} 


void add_array_points(FILE *pipe, double **points, int N, const char *config)
{ 
    fprintf(pipe, "\"<echo \'");
    for (int ii=0; ii<N; ii++) {
        fprintf(pipe, "%f %f", points[ii][0], points[ii][1]);
        if (ii < N-1) {
            fprintf(pipe, "\\n");
        }
    }
    fprintf(pipe, "\'\" %s, ", config);
}


void add_arrows_from_array_points(FILE* pipe, double **points, int N, int spacing, int offset, const char *config)
{   
    fprintf(pipe, "\"<echo \'");
    for (int ii=0; ii<N; ii++) {
        if ((offset +ii)%spacing == 0 && ii+1 < N) {
            fprintf(pipe, "%f %f %f %f", points[ii][0], points[ii][1], points[ii+1][0]-points[ii][0], points[ii+1][1]-points[ii][1]);
        }

        if (ii < N-1) {
            fprintf(pipe, "\\n");
        }
    }
    fprintf(pipe, "\'\" w vectors %s, ", config);

}


void add_yline(FILE *pipe, double y, const char *config) 
{
    fprintf(pipe, "\"<echo \'%f %f\\n%f %f'\" w l %s, ", 0.0, y, 1.0, y, config);
}


void add_line(FILE *pipe, double x0, double y0, double xf, double yf, const char *config)
{
    fprintf(pipe, "\"<echo \'%f %f %f %f\' \" w vectors nohead %s, ", x0, y0, xf-x0, yf-y0, config);
}


void add_function(FILE *pipe, double x0, double xf, int N, double (*fun)(double), const char *config)
{
    fprintf(pipe, "\"<echo \'");
    for (int ii=0; ii<N; ii++) {
        double x = x0 + ii * (xf - x0) / (N - 1);
        fprintf(pipe, "%f %f", x, fun(x));
        if (ii < N-1) {
            fprintf(pipe, "\\n");
        }
    }
    fprintf(pipe, "\'\" %s, ", config);
}


void plot_array_2D(double **array, int N, const char *title, char *xlabel, char *ylabel, const char *file_name, 
                   double *IC, double *plotDimensions_x0_xf_y0_yf, char *config, double *arrows_size_freq_offset)
{
    FILE *pipe = popen_gnuplot();
    char buffer[256];
    if (plotDimensions_x0_xf_y0_yf) {
        snprintf(buffer, 256, "set xrange[%f:%f]", plotDimensions_x0_xf_y0_yf[0], plotDimensions_x0_xf_y0_yf[1]);
        set_config(pipe, buffer);
        snprintf(buffer, 256, "set yrange[%f:%f]", plotDimensions_x0_xf_y0_yf[2], plotDimensions_x0_xf_y0_yf[3]);
        set_config(pipe, buffer);
    }
    if (xlabel) {
        snprintf(buffer, 256, "set xlabel '%s'", xlabel);
        set_config(pipe, buffer);
    }
    if (ylabel) {
        snprintf(buffer, 256, "set ylabel '%s'", ylabel);
        set_config(pipe, buffer);
    }
    set_config(pipe, "set grid");
    set_config(pipe, "set key box opaque");
    snprintf(buffer, 256, "set title \"%s\" font ',24'", title);
    set_config(pipe, buffer);

    start_plot(pipe, file_name);
    if (IC) {
        add_point(pipe, IC[0], IC[1], "ps 2 pt 7 lc 8 notitle");
    }
    add_array_points(pipe, array, N, config);
    if (arrows_size_freq_offset) {
        snprintf(buffer, 256, "lc 2 lw 3 size %f,20 fixed notitle", arrows_size_freq_offset[0]);
        add_arrows_from_array_points(pipe, array, N, arrows_size_freq_offset[1], arrows_size_freq_offset[2],
                                     buffer);
    }
    end_plot(pipe);
    pclose(pipe); 
}


