// dimer


#include "header.h"


int main(int argc, char **argv) {
    int size, rank;
    double delta_time, finish_time, h_max, lambda, sigma_1, sigma_2, start_time;
    const char *MPI_data_mu_1_mean_sigma_1_num_vs_h[N_MAX_PROC] = {
        "MPI_data/mu_1_mean_sigma_1_num_vs_h_part_1.txt",
        "MPI_data/mu_1_mean_sigma_1_num_vs_h_part_2.txt",
        "MPI_data/mu_1_mean_sigma_1_num_vs_h_part_3.txt",
        "MPI_data/mu_1_mean_sigma_1_num_vs_h_part_4.txt",
        "MPI_data/mu_1_mean_sigma_1_num_vs_h_part_5.txt",
        "MPI_data/mu_1_mean_sigma_1_num_vs_h_part_6.txt"
    };
    const char *MPI_data_mu_1_mean_sigma_2_num_vs_h[N_MAX_PROC] = {
        "MPI_data/mu_1_mean_sigma_2_num_vs_h_part_1.txt",
        "MPI_data/mu_1_mean_sigma_2_num_vs_h_part_2.txt",
        "MPI_data/mu_1_mean_sigma_2_num_vs_h_part_3.txt",
        "MPI_data/mu_1_mean_sigma_2_num_vs_h_part_4.txt",
        "MPI_data/mu_1_mean_sigma_2_num_vs_h_part_5.txt",
        "MPI_data/mu_1_mean_sigma_2_num_vs_h_part_6.txt"
    };
    const char *MPI_data_mu_2_inf_sigma_1_num_vs_h[N_MAX_PROC] = {
        "MPI_data/mu_2_inf_sigma_1_num_vs_h_part_1.txt",
        "MPI_data/mu_2_inf_sigma_1_num_vs_h_part_2.txt",
        "MPI_data/mu_2_inf_sigma_1_num_vs_h_part_3.txt",
        "MPI_data/mu_2_inf_sigma_1_num_vs_h_part_4.txt",
        "MPI_data/mu_2_inf_sigma_1_num_vs_h_part_5.txt",
        "MPI_data/mu_2_inf_sigma_1_num_vs_h_part_6.txt"
    };
    const char *MPI_data_mu_2_inf_sigma_2_num_vs_h[N_MAX_PROC] = {
        "MPI_data/mu_2_inf_sigma_2_num_vs_h_part_1.txt",
        "MPI_data/mu_2_inf_sigma_2_num_vs_h_part_2.txt",
        "MPI_data/mu_2_inf_sigma_2_num_vs_h_part_3.txt",
        "MPI_data/mu_2_inf_sigma_2_num_vs_h_part_4.txt",
        "MPI_data/mu_2_inf_sigma_2_num_vs_h_part_5.txt",
        "MPI_data/mu_2_inf_sigma_2_num_vs_h_part_6.txt"
    };
    const char *MPI_data_mu_2_mean_sigma_1_num_vs_h[N_MAX_PROC] = {
        "MPI_data/mu_2_mean_sigma_1_num_vs_h_part_1.txt",
        "MPI_data/mu_2_mean_sigma_1_num_vs_h_part_2.txt",
        "MPI_data/mu_2_mean_sigma_1_num_vs_h_part_3.txt",
        "MPI_data/mu_2_mean_sigma_1_num_vs_h_part_4.txt",
        "MPI_data/mu_2_mean_sigma_1_num_vs_h_part_5.txt",
        "MPI_data/mu_2_mean_sigma_1_num_vs_h_part_6.txt"
    };
    const char *MPI_data_mu_2_mean_sigma_2_num_vs_h[N_MAX_PROC] = {
        "MPI_data/mu_2_mean_sigma_2_num_vs_h_part_1.txt",
        "MPI_data/mu_2_mean_sigma_2_num_vs_h_part_2.txt",
        "MPI_data/mu_2_mean_sigma_2_num_vs_h_part_3.txt",
        "MPI_data/mu_2_mean_sigma_2_num_vs_h_part_4.txt",
        "MPI_data/mu_2_mean_sigma_2_num_vs_h_part_5.txt",
        "MPI_data/mu_2_mean_sigma_2_num_vs_h_part_6.txt"
    };
    const char *MPI_data_mu_2_mean_ac_sigma_1_num_vs_h[N_MAX_PROC] = {
        "MPI_data/mu_2_mean_ac_sigma_1_num_vs_h_part_1.txt",
        "MPI_data/mu_2_mean_ac_sigma_1_num_vs_h_part_2.txt",
        "MPI_data/mu_2_mean_ac_sigma_1_num_vs_h_part_3.txt",
        "MPI_data/mu_2_mean_ac_sigma_1_num_vs_h_part_4.txt",
        "MPI_data/mu_2_mean_ac_sigma_1_num_vs_h_part_5.txt",
        "MPI_data/mu_2_mean_ac_sigma_1_num_vs_h_part_6.txt"
    };
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    lambda = 8.0; // dipole-dipole interaction parameter
    sigma_1 = 10.0; // anisotropy parameter
    sigma_2 = 20.0; // anisotropy parameter
    h_max = 2.0;
    start_time = MPI_Wtime();
    if (rank == 0) {
        printf("---------------------------\n\n\n");
        write_mu_mean_approx_vs_h(h_max);
        write_mu_inf_approx_vs_h(h_max);
        write_mu_1_inf_num_vs_h(h_max, sigma_1, "data/mu_1_inf_sigma_1_num_vs_h.txt");
        write_mu_1_inf_num_vs_h(h_max, sigma_2, "data/mu_1_inf_sigma_2_num_vs_h.txt");
    }
    MPI_write_mu_1_mean_num_vs_h(size, rank, h_max, sigma_1, MPI_data_mu_1_mean_sigma_1_num_vs_h);
    MPI_write_mu_1_mean_num_vs_h(size, rank, h_max, sigma_2, MPI_data_mu_1_mean_sigma_2_num_vs_h);
    MPI_write_mu_2_inf_num_vs_h(size, rank, h_max, lambda, sigma_1, MPI_data_mu_2_inf_sigma_1_num_vs_h);
    MPI_write_mu_2_inf_num_vs_h(size, rank, h_max, lambda, sigma_2, MPI_data_mu_2_inf_sigma_2_num_vs_h);
    MPI_write_mu_2_mean_num_vs_h(size, rank, h_max, lambda, sigma_1, MPI_data_mu_2_mean_sigma_1_num_vs_h);
    MPI_write_mu_2_mean_num_vs_h(size, rank, h_max, lambda, sigma_2, MPI_data_mu_2_mean_sigma_2_num_vs_h);
    MPI_write_mu_2_mean_ac_vs_h(size, rank, h_max, lambda, sigma_1, MPI_data_mu_2_mean_ac_sigma_1_num_vs_h);
    MPI_Barrier(MPI_COMM_WORLD);
    finish_time = MPI_Wtime();
    delta_time = finish_time - start_time;
    if (rank == 0) {
        collect_y_vs_x(size, MPI_data_mu_1_mean_sigma_1_num_vs_h, "data/mu_1_mean_sigma_1_num_vs_h.txt");
        collect_y_vs_x(size, MPI_data_mu_1_mean_sigma_2_num_vs_h, "data/mu_1_mean_sigma_2_num_vs_h.txt");
        collect_y_vs_x(size, MPI_data_mu_2_inf_sigma_1_num_vs_h, "data/mu_2_inf_sigma_1_num_vs_h.txt");
        collect_y_vs_x(size, MPI_data_mu_2_inf_sigma_2_num_vs_h, "data/mu_2_inf_sigma_2_num_vs_h.txt");
        collect_y_vs_x(size, MPI_data_mu_2_mean_sigma_1_num_vs_h, "data/mu_2_mean_sigma_1_num_vs_h.txt");
        collect_y_vs_x(size, MPI_data_mu_2_mean_sigma_2_num_vs_h, "data/mu_2_mean_sigma_2_num_vs_h.txt");
        collect_y_vs_x(size, MPI_data_mu_2_mean_ac_sigma_1_num_vs_h, "data/mu_2_mean_ac_sigma_1_num_vs_h.txt");
        printf("C runtime = %.2lf sec = %.2lf min = %.2lf hours\n\n", delta_time, delta_time / 60.0, delta_time / 3600.0);
    }
    MPI_Finalize();
    return 0;
}
