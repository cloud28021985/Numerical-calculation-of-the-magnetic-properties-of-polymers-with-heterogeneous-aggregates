// dimer


#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


#define N_MAX_PROC 6 // max number of processors
#define N_plot_num 6 // number of points for the plot
#define N_plot_theory 200 // number of points for the plot
#define N_Riemann_1 100 // number of nodes in numerical integration
#define N_Riemann_3 40 // number of nodes in numerical integration
#define N_Riemann 4 // number of nodes in numerical integration
#define pi M_PI // the number pi


double int_pow(double x, int y);
double L(double x);
void collect_y_vs_x(int size, const char *str_1[], const char *str_2);
void MPI_write_mu_1_mean_num_vs_h(int size, int rank, double h_max, double sigma, const char *str[]);
void MPI_write_mu_2_inf_num_vs_h(int size, int rank, double h_max, double lambda, double sigma, const char *str[]);
void MPI_write_mu_2_mean_ac_vs_h(int size, int rank, double h_max, double lambda, double sigma, const char *str[]);
void MPI_write_mu_2_mean_num_vs_h(int size, int rank, double h_max, double lambda, double sigma, const char *str[]);
void worker_f_0_ac(double *f_0, double lambda, double phi_1, double phi_2, double sigma, double tau_1, double tau_2, double Z_0);
void worker_f_f_0_ac(double *f, double lambda, double phi_1, double phi_2, double sigma, double tau_1, double tau_2,
    double theta_1, double theta_2, double varphi_1, double varphi_2, double Z_0);
void worker_f_mu_1(double *f_mu, double h, double psi, double sigma, double theta, double varphi);
void worker_f_mu_1_inf(double *f_mu, double h, double sigma, double theta);
void worker_f_mu_1_mean(double *f_mu, double h, double psi, double sigma);
void worker_f_mu_2(double *f_mu, double h, double lambda, double psi, double sigma, double theta_1, double theta_2,
    double varphi_1, double varphi_2);
void worker_f_mu_2_ac(double *f, double f_0, double h, double lambda, double phi_1, double phi_2, double psi, double sigma,
    double tau_1, double tau_2, double theta_1, double theta_2, double varphi_1, double varphi_2);
void worker_f_mu_2_inf(double delta_varphi, double *f_mu, double h, double lambda, double sigma,
    double theta_1, double theta_2);
void worker_f_mu_2_mean(double *f_mu, double h, double lambda, double psi, double sigma);
void worker_f_mu_2_mean_ac(double *f, double h, double lambda, double psi, double sigma, double Z_0);
void worker_f_mu_mean_approx(int n, double *f, double h, double x);
void worker_f_Z_0_ac(double *f, double lambda, double phi_1, double phi_2, double sigma, double tau_1, double tau_2,
    double theta_1, double theta_2, double varphi_1, double varphi_2);
void worker_h_cdot_mu(double h, double *h_cdot_mu, double psi, double theta, double varphi);
void worker_mu_1_num(double h, double *mu, double psi, double sigma);
void worker_mu_1_inf_num(double h, double *mu, double sigma);
void worker_mu_1_mean_num(double h, double *mu, double sigma);
void worker_mu_2_ac(double h, double lambda, double *mu, double psi, double sigma, double Z_0);
void worker_mu_2_num(double h, double lambda, double *mu, double psi, double sigma);
void worker_mu_2_mean_ac(double h, double *mu, double lambda, double sigma, double Z_0);
void worker_mu_2_inf_num(double h, double lambda, double *mu, double sigma);
void worker_mu_2_mean_num(double h, double *mu, double lambda, double sigma);
void worker_mu_cdot_nu(double *mu_cdot_nu, double phi, double tau, double theta, double varphi);
void worker_mu_inf_approx(int n, double h, double *mu);
void worker_mu_mean_approx(int n, double h, double *mu);
void worker_u_1(double h_cdot_mu, double sigma, double theta, double *u);
void worker_u_1_inf(double h, double sigma, double theta, double *u);
void worker_u_2(double h_cdot_mu_1, double h_cdot_mu_2, double lambda, double sigma, double theta_1,
    double theta_2, double *u, double varphi_1, double varphi_2);
void worker_u_2_ac(double h_cdot_mu_1, double h_cdot_mu_2, double lambda, double phi_1, double phi_2, double sigma,
    double tau_1, double tau_2, double theta_1, double theta_2, double *u, double varphi_1, double varphi_2);
void worker_u_2_0_ac(double lambda, double phi_1, double phi_2, double sigma, double tau_1, double tau_2,
    double theta_1, double theta_2, double *u, double varphi_1, double varphi_2);
void worker_u_2_inf(double delta_varphi, double h, double lambda, double sigma, double theta_1, double theta_2, double *u);
void worker_u_a(double sigma, double theta, double *u_a);
void worker_u_a_ac(double phi, double sigma, double tau, double theta, double *u, double varphi);
void worker_u_dd(double lambda, double theta_1, double theta_2, double *u_dd, double varphi_1, double varphi_2);
void worker_u_dd_inf(double delta_varphi, double lambda, double theta_1, double theta_2, double *u_dd);
void worker_u_Z(double h_cdot_mu, double *u_Z);
void worker_Z_0_ac(double lambda, double sigma, double *Z_0);
void worker_u_Z_inf(double h, double theta, double *u_Z);
void write_mu_1_inf_num_vs_h(double h_max, double sigma, const char *str);
void write_mu_inf_approx_vs_h(double h_max);
void write_mu_mean_approx_vs_h(double h_max);
