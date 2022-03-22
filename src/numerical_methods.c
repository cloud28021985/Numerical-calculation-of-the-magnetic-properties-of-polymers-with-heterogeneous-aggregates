// dimer


#include "header.h"


// exponentiation by squaring
double int_pow(double x, int y) {
    double res = 1.0;
    while (y) {
        if (y & 1) {
            res *= x;
        }
        y >>= 1;
        x *= x;
    }
    return res;
}


void worker_mu_1_num(double h, double *mu, double psi, double sigma) {
    int i, j;
    double h_theta, h_varphi, theta, theta_max, varphi, varphi_max, f_mu[2], Z[2];
    theta_max = pi;
    varphi_max = 2.0 * pi;
    h_theta = theta_max / N_Riemann_3;
    h_varphi = varphi_max / N_Riemann_3;
    Z[0] = 0.0;
    Z[1] = 0.0;
    for (i = 0; i < N_Riemann_3; ++i) {
        for (j = 0; j < N_Riemann_3; ++j) {
            theta = (0.5 + i) * h_theta;
            varphi = (0.5 + j) * h_varphi;
            worker_f_mu_1(f_mu, h, psi, sigma, theta, varphi);
            Z[0] += f_mu[0];
            Z[1] += f_mu[1];
        }
    }
    *mu = Z[1] / Z[0];
}


void worker_mu_1_inf_num(double h, double *mu, double sigma) {
    int i;
    double h_theta, theta, theta_max, f_mu[2], Z[2];
    theta_max = pi;
    h_theta = theta_max / N_Riemann_1;
    Z[0] = 0.0;
    Z[1] = 0.0;
    for (i = 0; i < N_Riemann_1; ++i) {
        theta = (0.5 + i) * h_theta;
        worker_f_mu_1_inf(f_mu, h, sigma, theta);
        Z[0] += f_mu[0];
        Z[1] += f_mu[1];
    }
    *mu = Z[1] / Z[0];
}


void worker_mu_1_mean_num(double h, double *mu, double sigma) {
    int i;
    double f_mu, h_psi, psi, psi_max, sum;
    psi_max = pi / 2.0;
    h_psi = psi_max / N_Riemann_3;
    sum = 0.0;
    for (i = 0; i < N_Riemann_3; ++i) {
        psi = (0.5 + i) * h_psi;
        worker_f_mu_1_mean(&f_mu, h, psi, sigma);
        sum += h_psi * f_mu;
    }
    *mu = sum;
}


void worker_mu_2_ac(double h, double lambda, double *mu, double psi, double sigma, double Z_0) {
    int i_1, i_2, j_1, j_2, k_1, k_2, n_1, n_2;
    double f_0, h_phi, h_tau, h_theta, h_varphi, phi_1, phi_2, phi_max, tau_1, tau_2, tau_max, theta_1, theta_2, theta_max,
        varphi_1, varphi_2, varphi_max, f[2], Z[2];
    theta_max = pi;
    varphi_max = 2.0 * pi;
    tau_max = pi / 2.0;
    phi_max = 2.0 * pi;
    h_theta = theta_max / N_Riemann;
    h_varphi = varphi_max / N_Riemann;
    h_tau = tau_max / N_Riemann;
    h_phi = phi_max / N_Riemann;
    Z[0] = 0.0;
    Z[1] = 0.0;
    for(k_1 = 0; k_1 < N_Riemann; ++k_1) {
        tau_1 = (0.5 + k_1) * h_tau;
        for(k_2 = 0; k_2 < N_Riemann; ++k_2) {
            tau_2 = (0.5 + k_2) * h_tau;
            for(n_1 = 0; n_1 < N_Riemann; ++n_1) {
                phi_1 = (0.5 + n_1) * h_phi;
                for(n_2 = 0; n_2 < N_Riemann; ++n_2) {
                    phi_2 = (0.5 + n_2) * h_phi;
                    worker_f_0_ac(&f_0, lambda, phi_1, phi_2, sigma, tau_1, tau_2, Z_0);
                    for(i_1 = 0; i_1 < N_Riemann; ++i_1) {
                        theta_1 = (0.5 + i_1) * h_theta;
                        for(i_2 = 0; i_2 < N_Riemann; ++i_2) {
                            theta_2 = (0.5 + i_2) * h_theta;
                            for(j_1 = 0; j_1 < N_Riemann; ++j_1) {
                                varphi_1 = (0.5 + j_1) * h_varphi;
                                for(j_2 = 0; j_2 < N_Riemann; ++j_2) {
                                    varphi_2 = (0.5 + j_2) * h_varphi;
                                    worker_f_mu_2_ac(f, f_0, h, lambda, phi_1, phi_2, psi, sigma, tau_1, tau_2,
                                        theta_1, theta_2, varphi_1, varphi_2);
                                    Z[0] += f[0];
                                    Z[1] += f[1];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    *mu = Z[1] / Z[0];
}


void worker_mu_2_inf_num(double h, double lambda, double *mu, double sigma) {
    int i_1, i_2, j;
    double delta_varphi, h_theta, h_varphi, theta_1, theta_2, theta_max, varphi_max, f_mu[2], Z[2];
    theta_max = pi;
    varphi_max = 2.0 * pi;
    h_theta = theta_max / N_Riemann_3;
    h_varphi = varphi_max / N_Riemann_3;
    Z[0] = 0.0;
    Z[1] = 0.0;
    for (i_1 = 0; i_1 < N_Riemann_3; ++i_1) {
        for (i_2 = 0; i_2 < N_Riemann_3; ++i_2) {
            for (j = 0; j < N_Riemann_3; ++j) {
                theta_1 = (0.5 + i_1) * h_theta;
                theta_2 = (0.5 + i_2) * h_theta;
                delta_varphi = (0.5 + j) * h_varphi;
                worker_f_mu_2_inf(delta_varphi, f_mu, h, lambda, sigma, theta_1, theta_2);
                Z[0] += f_mu[0];
                Z[1] += f_mu[1];
            }
        }
    }
    *mu = Z[1] / Z[0];
}


void worker_mu_2_mean_ac(double h, double *mu, double lambda, double sigma, double Z_0) {
    int i;
    double f, h_psi, psi, psi_max, sum;
    psi_max = pi / 2.0;
    h_psi = psi_max / N_Riemann;
    sum = 0.0;
    for (i = 0; i < N_Riemann; ++i) {
        psi = (0.5 + i) * h_psi;
        worker_f_mu_2_mean_ac(&f, h, lambda, psi, sigma, Z_0);
        sum += h_psi * f;
    }
    *mu = sum;
}


void worker_mu_2_mean_num(double h, double *mu, double lambda, double sigma) {
    int i;
    double f_mu, h_psi, psi, psi_max, sum;
    psi_max = pi / 2.0;
    h_psi = psi_max / N_Riemann;
    sum = 0.0;
    for (i = 0; i < N_Riemann; ++i) {
        psi = (0.5 + i) * h_psi;
        worker_f_mu_2_mean(&f_mu, h, lambda, psi, sigma);
        sum += h_psi * f_mu;
    }
    *mu = sum;
}


void worker_mu_mean_approx(int n, double h, double *mu) {
    int i;
    double f, h_int, sum, x;
    h_int = 1.0 / N_Riemann_1;
    sum = 0.0;
    for (i = 0; i < N_Riemann_1; ++i) {
        x = (0.5 + i) * h_int;
        worker_f_mu_mean_approx(n, &f, h, x);
        sum += h_int * f;
    }
    *mu = sum;
}


void worker_mu_2_num(double h, double lambda, double *mu, double psi, double sigma) {
    int i_1, i_2, j_1, j_2;
    double h_theta, h_varphi, theta_1, theta_2, theta_max, varphi_1, varphi_2, varphi_max, f_mu[2], Z[2];
    theta_max = pi;
    varphi_max = 2.0 * pi;
    h_theta = theta_max / N_Riemann;
    h_varphi = varphi_max / N_Riemann;
    Z[0] = 0.0;
    Z[1] = 0.0;
    for (i_1 = 0; i_1 < N_Riemann; ++i_1) {
        for (i_2 = 0; i_2 < N_Riemann; ++i_2) {
            for (j_1 = 0; j_1 < N_Riemann; ++j_1) {
                for (j_2 = 0; j_2 < N_Riemann; ++j_2) {
                    theta_1 = (0.5 + i_1) * h_theta;
                    theta_2 = (0.5 + i_2) * h_theta;
                    varphi_1 = (0.5 + j_1) * h_varphi;
                    varphi_2 = (0.5 + j_1) * h_varphi;
                    worker_f_mu_2(f_mu, h, lambda, psi, sigma, theta_1, theta_2, varphi_1, varphi_2);
                    Z[0] += f_mu[0];
                    Z[1] += f_mu[1];
                }
            }
        }
    }
    *mu = Z[1] / Z[0];
}


void worker_Z_0_ac(double lambda, double sigma, double *Z_0) {
    int i_1, i_2, j_1, j_2, k_1, k_2, n_1, n_2;
    double f, h_phi, h_tau, h_theta, h_varphi, phi_1, phi_2, phi_max, sum, tau_1, tau_2, tau_max, theta_1, theta_2, theta_max,
        varphi_1, varphi_2, varphi_max;
    theta_max = pi;
    varphi_max = 2.0 * pi;
    tau_max = pi / 2.0;
    phi_max = 2.0 * pi;
    h_theta = theta_max / N_Riemann;
    h_varphi = varphi_max / N_Riemann;
    h_tau = tau_max / N_Riemann;
    h_phi = phi_max / N_Riemann;
    sum = 0.0;
    for(k_1 = 0; k_1 < N_Riemann; ++k_1) {
        tau_1 = (0.5 + k_1) * h_tau;
        for(k_2 = 0; k_2 < N_Riemann; ++k_2) {
            tau_2 = (0.5 + k_2) * h_tau;
            for(n_1 = 0; n_1 < N_Riemann; ++n_1) {
                phi_1 = (0.5 + n_1) * h_phi;
                for(n_2 = 0; n_2 < N_Riemann; ++n_2) {
                    phi_2 = (0.5 + n_2) * h_phi;
                    for(i_1 = 0; i_1 < N_Riemann; ++i_1) {
                        theta_1 = (0.5 + i_1) * h_theta;
                        for(i_2 = 0; i_2 < N_Riemann; ++i_2) {
                            theta_2 = (0.5 + i_2) * h_theta;
                            for(j_1 = 0; j_1 < N_Riemann; ++j_1) {
                                varphi_1 = (0.5 + j_1) * h_varphi;
                                for(j_2 = 0; j_2 < N_Riemann; ++j_2) {
                                    varphi_2 = (0.5 + j_2) * h_varphi;
                                    worker_f_Z_0_ac(&f, lambda, phi_1, phi_2, sigma, tau_1, tau_2, theta_1, theta_2,
                                        varphi_1, varphi_2);
                                    sum += int_pow(h_theta * h_varphi * h_tau * h_phi, 2) * f;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    *Z_0 = sum;
}
