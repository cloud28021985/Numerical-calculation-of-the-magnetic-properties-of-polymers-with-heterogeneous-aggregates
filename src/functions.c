// dimer


#include "header.h"


double coth(double x) {
    return (exp(x) + exp(- x)) / (exp(x) - exp(- x));
}


double L(double x) {
    return coth(x) - 1.0 / x;
}


void worker_f_0_ac(double *f_0, double lambda, double phi_1, double phi_2, double sigma, double tau_1, double tau_2,
    double Z_0) {
    int i_1, i_2, j_1, j_2;
    double f, h_theta, h_varphi, sum, theta_1, theta_2, theta_max, varphi_1, varphi_2, varphi_max;
    theta_max = pi;
    varphi_max = 2.0 * pi;
    h_theta = theta_max / N_Riemann;
    h_varphi = varphi_max / N_Riemann;
    sum = 0.0;
    for(i_1 = 0; i_1 < N_Riemann; ++i_1) {
        theta_1 = (0.5 + i_1) * h_theta;
        for(i_2 = 0; i_2 < N_Riemann; ++i_2) {
            theta_2 = (0.5 + i_2) * h_theta;
            for(j_1 = 0; j_1 < N_Riemann; ++j_1) {
                varphi_1 = (0.5 + j_1) * h_varphi;
                for(j_2 = 0; j_2 < N_Riemann; ++j_2) {
                    varphi_2 = (0.5 + j_2) * h_varphi;
                    worker_f_f_0_ac(&f, lambda, phi_1, phi_2, sigma, tau_1, tau_2, theta_1, theta_2, varphi_1, varphi_2, Z_0);
                    sum += int_pow(h_theta * h_varphi, 2) * f;
                }
            }
        }
    }
    *f_0 = sum;
}


void worker_f_f_0_ac(double *f, double lambda, double phi_1, double phi_2, double sigma, double tau_1, double tau_2,
    double theta_1, double theta_2, double varphi_1, double varphi_2, double Z_0) {
    double u_0;
    worker_u_2_0_ac(lambda, phi_1, phi_2, sigma, tau_1, tau_2, theta_1, theta_2, &u_0, varphi_1, varphi_2);
    *f = exp(- u_0) * sin(theta_1) * sin(theta_2) / Z_0;
}


void worker_f_mu_1(double *f_mu, double h, double psi, double sigma, double theta, double varphi) {
    double h_cdot_mu, u;
    worker_h_cdot_mu(h, &h_cdot_mu, psi, theta, varphi);
    worker_u_1(h_cdot_mu, sigma, theta, &u);
    f_mu[0] = exp(- u) * sin(theta);
    f_mu[1] = h_cdot_mu / h * f_mu[0];
}


void worker_f_mu_1_inf(double *f_mu, double h, double sigma, double theta) {
    double u;
    worker_u_1_inf(h, sigma, theta, &u);
    f_mu[0] = exp(- u) * sin(theta);
    f_mu[1] = cos(theta) * f_mu[0];
}


void worker_f_mu_1_mean(double *f_mu, double h, double psi, double sigma) {
    double mu;
    worker_mu_1_num(h, &mu, psi, sigma);
    *f_mu = mu * sin(psi);
}


void worker_f_mu_2(double *f_mu, double h, double lambda, double psi, double sigma, double theta_1, double theta_2,
    double varphi_1, double varphi_2) {
    double h_cdot_mu_1, h_cdot_mu_2, u;
    worker_h_cdot_mu(h, &h_cdot_mu_1, psi, theta_1, varphi_1);
    worker_h_cdot_mu(h, &h_cdot_mu_2, psi, theta_2, varphi_2);
    worker_u_2(h_cdot_mu_1, h_cdot_mu_2, lambda, sigma, theta_1, theta_2, &u, varphi_1, varphi_2);
    f_mu[0] = exp(- u) * sin(theta_1) * sin(theta_2);
    f_mu[1] = h_cdot_mu_1 / h * f_mu[0];
}


void worker_f_mu_2_ac(double *f, double f_0, double h, double lambda, double phi_1, double phi_2, double psi, double sigma,
    double tau_1, double tau_2, double theta_1, double theta_2, double varphi_1, double varphi_2) {
    double h_cdot_mu_1, h_cdot_mu_2, u;
    worker_h_cdot_mu(h, &h_cdot_mu_1, psi, theta_1, varphi_1);
    worker_h_cdot_mu(h, &h_cdot_mu_2, psi, theta_2, varphi_2);
    worker_u_2_ac(h_cdot_mu_1, h_cdot_mu_2, lambda, phi_1, phi_2, sigma, tau_1, tau_2, theta_1, theta_2, &u, varphi_1,
        varphi_2);
    f[0] = exp(- u) * f_0 * sin(theta_1) * sin(theta_2) * sin(tau_1) * sin(tau_2);
    f[1] = h_cdot_mu_1 / h * f[0];
}


void worker_f_mu_2_inf(double delta_varphi, double *f_mu, double h, double lambda, double sigma,
    double theta_1, double theta_2) {
    double u;
    worker_u_2_inf(delta_varphi, h, lambda, sigma, theta_1, theta_2, &u);
    f_mu[0] = exp(- u) * sin(theta_1) * sin(theta_2);
    f_mu[1] = cos(theta_1) * f_mu[0];
}


void worker_f_mu_2_mean(double *f_mu, double h, double lambda, double psi, double sigma) {
    double mu;
    worker_mu_2_num(h, lambda, &mu, psi, sigma);
    *f_mu = mu * sin(psi);
}


void worker_f_mu_2_mean_ac(double *f, double h, double lambda, double psi, double sigma, double Z_0) {
    double mu;
    worker_mu_2_ac(h, lambda, &mu, psi, sigma, Z_0);
    *f = mu * sin(psi);
}


void worker_f_mu_mean_approx(int n, double *f, double h, double x) {
    *f = x * tanh(n * h * x);
}


void worker_f_Z_0_ac(double *f, double lambda, double phi_1, double phi_2, double sigma, double tau_1, double tau_2,
    double theta_1, double theta_2, double varphi_1, double varphi_2) {
    double u_0;
    worker_u_2_0_ac(lambda, phi_1, phi_2, sigma, tau_1, tau_2, theta_1, theta_2, &u_0, varphi_1, varphi_2);
    *f = exp(- u_0) * sin(theta_1) * sin(theta_2) * sin(tau_1) * sin(tau_2);
}


void worker_h_cdot_mu(double h, double *h_cdot_mu, double psi, double theta, double varphi) {
    *h_cdot_mu = h * (sin(psi) * sin(theta) * cos(varphi) + cos(psi) * cos(theta));
}


void worker_mu_cdot_nu(double *mu_cdot_nu, double phi, double tau, double theta, double varphi) {
    *mu_cdot_nu = sin(theta) * sin(tau) * cos(varphi - phi) + cos(theta) * cos(tau);
}


void worker_mu_inf_approx(int n, double h, double *mu) {
    *mu = tanh(n * h);
}


void worker_u_1(double h_cdot_mu, double sigma, double theta, double *u) {
    double u_a, u_Z;
    worker_u_Z(h_cdot_mu, &u_Z);
    worker_u_a(sigma, theta, &u_a);
    *u = u_Z + u_a;
}


void worker_u_1_ac(double h_cdot_mu, double phi, double sigma, double tau, double theta, double *u, double varphi) {
    double u_a, u_Z;
    worker_u_Z(h_cdot_mu, &u_Z);
    worker_u_a_ac(phi, sigma, tau, theta, &u_a, varphi);
    *u = u_Z + u_a;
}


void worker_u_1_inf(double h, double sigma, double theta, double *u) {
    double u_a, u_Z;
    worker_u_Z_inf(h, theta, &u_Z);
    worker_u_a(sigma, theta, &u_a);
    *u = u_Z + u_a;
}


void worker_u_2(double h_cdot_mu_1, double h_cdot_mu_2, double lambda, double sigma, double theta_1,
    double theta_2, double *u, double varphi_1, double varphi_2) {
    double u_1, u_2, u_dd;
    worker_u_1(h_cdot_mu_1, sigma, theta_1, &u_1);
    worker_u_1(h_cdot_mu_2, sigma, theta_2, &u_2);
    worker_u_dd(lambda, theta_1, theta_2, &u_dd, varphi_1, varphi_2);
    *u = u_1 + u_2 + u_dd;
}


void worker_u_2_0_ac(double lambda, double phi_1, double phi_2, double sigma, double tau_1, double tau_2,
    double theta_1, double theta_2, double *u, double varphi_1, double varphi_2) {
    double u_a_1, u_a_2, u_dd;
    worker_u_a_ac(phi_1, sigma, tau_1, theta_1, &u_a_1, varphi_1);
    worker_u_a_ac(phi_2, sigma, tau_2, theta_2, &u_a_2, varphi_2);
    worker_u_dd(lambda, theta_1, theta_2, &u_dd, varphi_1, varphi_2);
    *u = u_a_1 + u_a_2 + u_dd;
}


void worker_u_2_ac(double h_cdot_mu_1, double h_cdot_mu_2, double lambda, double phi_1, double phi_2, double sigma,
    double tau_1, double tau_2, double theta_1, double theta_2, double *u, double varphi_1, double varphi_2) {
    double u_1, u_2, u_dd;
    worker_u_1_ac(h_cdot_mu_1, phi_1, sigma, tau_1, theta_1, &u_1, varphi_1);
    worker_u_1_ac(h_cdot_mu_2, phi_2, sigma, tau_2, theta_2, &u_2, varphi_2);
    worker_u_dd(lambda, theta_1, theta_2, &u_dd, varphi_1, varphi_2);
    *u = u_1 + u_2 + u_dd;
}


void worker_u_2_inf(double delta_varphi, double h, double lambda, double sigma, double theta_1, double theta_2, double *u) {
    double u_1, u_2, u_dd;
    worker_u_1_inf(h, sigma, theta_1, &u_1);
    worker_u_1_inf(h, sigma, theta_2, &u_2);
    worker_u_dd_inf(delta_varphi, lambda, theta_1, theta_2, &u_dd);
    *u = u_1 + u_2 + u_dd;
}


void worker_u_a(double sigma, double theta, double *u_a) {
    *u_a = - sigma * int_pow(cos(theta), 2);
}


void worker_u_a_ac(double phi, double sigma, double tau, double theta, double *u, double varphi) {
    double mu_cdot_nu;
    worker_mu_cdot_nu(&mu_cdot_nu, phi, tau, theta, varphi);
    *u = - sigma * int_pow(mu_cdot_nu, 2);
}


void worker_u_dd(double lambda, double theta_1, double theta_2, double *u_dd, double varphi_1, double varphi_2) {
    *u_dd = lambda * (sin(theta_1) * sin(theta_2) * cos(varphi_1 - varphi_2) - 2.0 * cos(theta_1) * cos(theta_2));
}


void worker_u_dd_inf(double delta_varphi, double lambda, double theta_1, double theta_2, double *u_dd) {
    *u_dd = lambda * (sin(theta_1) * sin(theta_2) * cos(delta_varphi) - 2.0 * cos(theta_1) * cos(theta_2));
}


void worker_u_Z(double h_cdot_mu, double *u_Z) {
    *u_Z = - h_cdot_mu;
}


void worker_u_Z_inf(double h, double theta, double *u_Z) {
    *u_Z = - h * cos(theta);
}
