# dimer


import matplotlib.pyplot as plt
import numpy as np
import timeit
start = timeit.default_timer()


h_apr, L, mu_1_apr, mu_2_apr = np.loadtxt('data/mu_mean_approx_vs_h.txt', comments = '#', unpack = True)
h_2_num, mu_2_num = np.loadtxt('data/mu_2_mean_sigma_1_num_vs_h.txt', comments = '#', unpack = True)
h_2_ac, mu_2_ac = np.loadtxt('data/mu_2_mean_ac_sigma_1_num_vs_h.txt', comments = '#', unpack = True)
fig = plt.figure(figsize = (6.0, 4.5))
plt.plot(h_apr, mu_2_apr, color = 'black', linewidth = 2.0)
plt.plot(h_2_ac, mu_2_ac, linestyle = '--', marker = 'o', color = 'black', linewidth = 1.0, markersize = 6.0)
plt.plot(h_2_num, mu_2_num, linestyle = '--', marker = 's', color = 'black', linewidth = 1.0, markersize = 6.0)
plt.xlabel('$h$', fontsize = 16)
plt.ylabel('$<μ_h>$', fontsize = 16)
plt.tight_layout()
plt.grid(True)
plt.savefig('figs/comparison.pdf')
plt.close('all')


h_apr, L, mu_1_apr, mu_2_apr = np.loadtxt('data/mu_mean_approx_vs_h.txt', comments = '#', unpack = True)
fig = plt.figure(figsize = (6.0, 4.5))
plt.plot(h_apr, mu_2_apr, color = 'black', linewidth = 2.0)
plt.plot(h_apr, mu_1_apr, linestyle = '--', color = 'black', linewidth = 2.0)
plt.plot(h_apr, L, linestyle = '-.', color = 'black', linewidth = 2.0)
plt.xlabel('$h$', fontsize = 16)
plt.ylabel('$<μ_h>$', fontsize = 16)
plt.tight_layout()
plt.grid(True)
plt.savefig('figs/asymptotics.pdf')
plt.close('all')


h_1_num, mu_1_num = np.loadtxt('data/mu_1_mean_sigma_1_num_vs_h.txt', comments = '#', unpack = True)
h_2_num, mu_2_num = np.loadtxt('data/mu_2_mean_sigma_1_num_vs_h.txt', comments = '#', unpack = True)
fig = plt.figure(figsize = (6.0, 4.5))
plt.plot(h_apr, mu_1_apr, linestyle = '--', color = 'black', linewidth = 2.0)
plt.plot(h_apr, mu_2_apr, color = 'black', linewidth = 2.0)
plt.plot(h_1_num, mu_1_num, linestyle = '--', marker = 'o', color = 'black', linewidth = 1.0, markersize = 6.0)
plt.plot(h_2_num, mu_2_num, linestyle = '--', marker = 's', color = 'black', linewidth = 1.0, markersize = 6.0)
plt.xlabel('$h$', fontsize = 16)
plt.ylabel('$<μ_h>$', fontsize = 16)
plt.tight_layout()
plt.grid(True)
plt.savefig('figs/mean_10.pdf')
plt.close('all')


h_1_num, mu_1_num = np.loadtxt('data/mu_1_mean_sigma_2_num_vs_h.txt', comments = '#', unpack = True)
h_2_num, mu_2_num = np.loadtxt('data/mu_2_mean_sigma_2_num_vs_h.txt', comments = '#', unpack = True)
fig = plt.figure(figsize = (6.0, 4.5))
plt.plot(h_apr, mu_1_apr, linestyle = '--', color = 'black', linewidth = 2.0)
plt.plot(h_apr, mu_2_apr, color = 'black', linewidth = 2.0)
plt.plot(h_1_num, mu_1_num, linestyle = '--', marker = 'o', color = 'black', linewidth = 1.0, markersize = 6.0)
plt.plot(h_2_num, mu_2_num, linestyle = '--', marker = 's', color = 'black', linewidth = 1.0, markersize = 6.0)
plt.xlabel('$h$', fontsize = 16)
plt.ylabel('$<μ_h>$', fontsize = 16)
plt.tight_layout()
plt.grid(True)
plt.savefig('figs/mean_20.pdf')
plt.close('all')


h_apr, mu_1_apr, mu_2_apr = np.loadtxt('data/mu_inf_approx_vs_h.txt', comments = '#', unpack = True)
h_1_num, mu_1_num = np.loadtxt('data/mu_1_inf_sigma_1_num_vs_h.txt', comments = '#', unpack = True)
h_2_num, mu_2_num = np.loadtxt('data/mu_2_inf_sigma_1_num_vs_h.txt', comments = '#', unpack = True)
fig = plt.figure(figsize = (6.0, 4.5))
plt.plot(h_apr, mu_1_apr, linestyle = '--', color = 'black', linewidth = 2.0)
plt.plot(h_apr, mu_2_apr, color = 'black', linewidth = 2.0)
plt.plot(h_1_num, mu_1_num, linestyle = '--', marker = 'o', color = 'black', linewidth = 1.0, markersize = 6.0)
plt.plot(h_2_num, mu_2_num, linestyle = '--', marker = 's', color = 'black', linewidth = 1.0, markersize = 6.0)
plt.xlabel('$h$', fontsize = 16)
plt.ylabel('$μ_h$', fontsize = 16)
plt.tight_layout()
plt.grid(True)
plt.savefig('figs/inf_10.pdf')
plt.close('all')


h_1_num, mu_1_num = np.loadtxt('data/mu_1_inf_sigma_2_num_vs_h.txt', comments = '#', unpack = True)
h_2_num, mu_2_num = np.loadtxt('data/mu_2_inf_sigma_2_num_vs_h.txt', comments = '#', unpack = True)
fig = plt.figure(figsize = (6.0, 4.5))
plt.plot(h_apr, mu_1_apr, linestyle = '--', color = 'black', linewidth = 2.0)
plt.plot(h_apr, mu_2_apr, color = 'black', linewidth = 2.0)
plt.plot(h_1_num, mu_1_num, linestyle = '--', marker = 'o', color = 'black', linewidth = 1.0, markersize = 6.0)
plt.plot(h_2_num, mu_2_num, linestyle = '--', marker = 's', color = 'black', linewidth = 1.0, markersize = 6.0)
plt.xlabel('$h$', fontsize = 16)
plt.ylabel('$μ_h$', fontsize = 16)
plt.tight_layout()
plt.grid(True)
plt.savefig('figs/inf_20.pdf')
plt.close('all')


stop = timeit.default_timer()
runtime = stop - start
print('plot runtime = {:4.2f} sec\n' .format(runtime))
