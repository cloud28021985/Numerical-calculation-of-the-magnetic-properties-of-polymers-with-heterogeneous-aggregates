// dimer


#include "header.h"


void collect_y_vs_x(int size, const char *str_1[], const char *str_2) {
    int i, n;
    double x, y, array_x[N_plot_num + 1], array_y[N_plot_num + 1];
    n = 0;
    FILE *myfile_1 = NULL;
    for (i = 0; i <= size - 1; ++i) {
        myfile_1 = fopen(str_1[i], "r");
        while(fscanf(myfile_1, "%10lf %20lf\n", &x, &y) != EOF) {
            array_x[n] = x;
            array_y[n] = y;
            ++n;
        }
        fclose(myfile_1);
    }
    FILE *myfile_2 = NULL;
    myfile_2 = fopen(str_2, "w");
    for(i = 0; i <= n - 1; ++i) {
        fprintf(myfile_2, "%10.6lf %20.10lf\n", array_x[i], array_y[i]);
    }
    fclose(myfile_2);
}


void MPI_write_mu_1_mean_num_vs_h(int size, int rank, double h_max, double sigma, const char *str[]) {
    int i, j;
    double h, h_plot, mu;
    h_plot = h_max / N_plot_num;
    for (i = 0; i <= size - 1; ++i) {
        if (rank == i) {
            h = i * h_max / size;
            FILE *myfile = NULL;
            myfile = fopen(str[i], "w");
            if (i == 0) {
                fprintf(myfile, "%10.6lf %20.10lf\n", 0.0, 0.0);
            }
            for(j = 0; j < N_plot_num / size; ++j) {
                h += h_plot;
                worker_mu_1_mean_num(h, &mu, sigma);
                fprintf(myfile, "%10.6lf %20.10lf\n", h, mu);
            }
            fclose(myfile);
        }
    }
}


void MPI_write_mu_2_inf_num_vs_h(int size, int rank, double h_max, double lambda, double sigma, const char *str[]) {
    int i, j;
    double h, h_plot, mu;
    h_plot = h_max / N_plot_num;
    for (i = 0; i <= size - 1; ++i) {
        if (rank == i) {
            h = i * h_max / size;
            FILE *myfile = NULL;
            myfile = fopen(str[i], "w");
            if (i == 0) {
                fprintf(myfile, "%10.6lf %20.10lf\n", 0.0, 0.0);
            }
            for(j = 0; j < N_plot_num / size; ++j) {
                h += h_plot;
                worker_mu_2_inf_num(h, lambda, &mu, sigma);
                fprintf(myfile, "%10.6lf %20.10lf\n", h, mu);
            }
            fclose(myfile);
        }
    }
}


void MPI_write_mu_2_mean_ac_vs_h(int size, int rank, double h_max, double lambda, double sigma, const char *str[]) {
    int i, j;
    double h, h_plot, mu, Z_0;
    h_plot = h_max / N_plot_num;
    worker_Z_0_ac(lambda, sigma, &Z_0);
    for (i = 0; i <= size - 1; ++i) {
        if (rank == i) {
            h = i * h_max / size;
            FILE *myfile = NULL;
            myfile = fopen(str[i], "w");
            if (i == 0) {
                fprintf(myfile, "%10.6lf %20.10lf\n", 0.0, 0.0);
            }
            for(j = 0; j < N_plot_num / size; ++j) {
                h += h_plot;
                worker_mu_2_mean_ac(h, &mu, lambda, sigma, Z_0);
                fprintf(myfile, "%10.6lf %20.10lf\n", h, mu);
            }
            fclose(myfile);
        }
    }
}


void MPI_write_mu_2_mean_num_vs_h(int size, int rank, double h_max, double lambda, double sigma, const char *str[]) {
    int i, j;
    double h, h_plot, mu;
    h_plot = h_max / N_plot_num;
    for (i = 0; i <= size - 1; ++i) {
        if (rank == i) {
            h = i * h_max / size;
            FILE *myfile = NULL;
            myfile = fopen(str[i], "w");
            if (i == 0) {
                fprintf(myfile, "%10.6lf %20.10lf\n", 0.0, 0.0);
            }
            for(j = 0; j < N_plot_num / size; ++j) {
                h += h_plot;
                worker_mu_2_mean_num(h, &mu, lambda, sigma);
                fprintf(myfile, "%10.6lf %20.10lf\n", h, mu);
            }
            fclose(myfile);
        }
    }
}


void write_mu_1_inf_num_vs_h(double h_max, double sigma, const char *str) {
    int i;
    double h, h_plot, mu_inf_1;
    h_plot = h_max / N_plot_num;
    h = 0.0;
    FILE *myfile = NULL;
    myfile = fopen(str, "w");
    fprintf(myfile, "%10s %20s\n", "# h", "mu_inf_1");
    fprintf(myfile, "%10.6lf %20.6lf\n", 0.0, 0.0);
    for(i = 0; i < N_plot_num; ++i) {
        h += h_plot;
        worker_mu_1_inf_num(h, &mu_inf_1, sigma);
        fprintf(myfile, "%10.6lf %20.6lf\n", h, mu_inf_1);
    }
    fclose(myfile);
}


void write_mu_inf_approx_vs_h(double h_max) {
    int i;
    double h, h_plot, mu_inf_1, mu_inf_2;
    h_plot = h_max / N_plot_theory;
    h = 0.0;
    FILE *myfile = NULL;
    myfile = fopen("data/mu_inf_approx_vs_h.txt", "w");
    fprintf(myfile, "%10s %20s %20s\n", "# h", "mu_inf_1", "mu_inf_2");
    fprintf(myfile, "%10.6lf %20.6lf %20.6lf\n", 0.0, 0.0, 0.0);
    for(i = 0; i < N_plot_theory; ++i) {
        h += h_plot;
        worker_mu_inf_approx(1, h, &mu_inf_1);
        worker_mu_inf_approx(2, h, &mu_inf_2);
        fprintf(myfile, "%10.6lf %20.6lf %20.6lf\n", h, mu_inf_1, mu_inf_2);
    }
    fclose(myfile);
}


void write_mu_mean_approx_vs_h(double h_max) {
    int i;
    double h, h_plot, mu_mean_1, mu_mean_2;
    h_plot = h_max / N_plot_theory;
    h = 0.0;
    FILE *myfile = NULL;
    myfile = fopen("data/mu_mean_approx_vs_h.txt", "w");
    fprintf(myfile, "%10s %20s %20s %20s\n", "# h", "L(h)", "mu_mean_1", "mu_mean_2");
    fprintf(myfile, "%10.6lf %20.6lf %20.6lf %20.6lf\n", 0.0, 0.0, 0.0, 0.0);
    for(i = 0; i < N_plot_theory; ++i) {
        h += h_plot;
        worker_mu_mean_approx(1, h, &mu_mean_1);
        worker_mu_mean_approx(2, h, &mu_mean_2);
        fprintf(myfile, "%10.6lf %20.6lf %20.6lf %20.6lf\n", h, L(h), mu_mean_1, mu_mean_2);
    }
    fclose(myfile);
}
