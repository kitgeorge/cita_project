import numpy as np
import matplotlib.pyplot as plt
import math
import sys
sys.path.append('/home/kit/Documents/cita_project/')
from py_src.disc_plotter import DiscPlotter
import math

N_particles = 10000
N_timesteps = 1001

# n_max = 64
# l_max = 32

# bfe_coefficients = np.loadtxt("../data/n_body/test_bfe_coefficients.csv")
# bfe_coefficients = bfe_coefficients.reshape(N_timesteps, n_max + 1, l_max + 1, 2)

# for i in range(N_timesteps):
#     print(i)
#     plt.imshow(np.log(np.power(np.power(bfe_coefficients[i, :, :, 0], 2) + np.power(bfe_coefficients[i, :, :, 1], 2), 0.5)), origin="lower")
#     plt.savefig("../plots/n_body/bfe_coefficients/t={}.png".format(i), dpi=500)

# exit()

# v_0 = 0.22

# bin_interval = 0.005
# bin_max = 0.05
# N_bins = int(bin_max/bin_interval)

# theta_R_binned = [[] for k in range(N_bins)]

# init_theta_R_values = np.loadtxt("../data/n_body/init_theta_R_values.csv")
# plt.hist(init_theta_R_values, bins=1000)
# plt.savefig("../plots/init_theta_R_values.png")
# plt.close()

# theta_R_values = np.loadtxt("../data/n_body/theta_R_values.csv")
# theta_R_values = theta_R_values.reshape(N_particles, N_timesteps)

# def Phi_eff_min(v_0, L):
#     return v_0**2*(1/6 + math.log(math.sqrt(3)*abs(L)/v_0))


# E_L_values = np.loadtxt("../data/n_body/E_L_values.csv")
# E_L_values = E_L_values.reshape(N_particles, N_timesteps, 2)


# oscillator_E_values = np.array([[E_L[0] - Phi_eff_min(v_0, E_L[1])
#                                  for E_L in row] for row in E_L_values])


# for i in range(N_bins):
#     for j in range(N_particles):
#         if(oscillator_E_values[j, 0] >= i*bin_interval
#            and oscillator_E_values[j, 0] < (i + 1)*bin_interval):
#             theta_R_binned[i].append(theta_R_values[j, 0])

# for i in range(N_bins):
#     plt.hist(theta_R_binned[i], bins=10)
#     plt.savefig("../plots/n_body/sub_hists/E<{}.png".format((i + 1)*bin_interval),
#                 dpi=500)
#     plt.close()

# for i in range(N_timesteps):
    # plt.hist(E_L_values[:, i, 0], bins=10)
    # plt.savefig("../plots/n_body/E_distribution/t={}.png".format(i),
    #             dpi = 500)
    # plt.close()
    # plt.hist(E_L_values[:, i, 1], bins=10)
    # plt.savefig("../plots/n_body/L_distribution/t={}.png".format(i),
    #             dpi = 500)
    # plt.close()

# for i in range(N_particles):
#     print(i)
#     plt.plot(E_L_values[i, :, 0])
#     plt.savefig("../plots/n_body/E_evolution/particle_{}.png".format(i))
#     plt.close()
#     plt.plot(E_L_values[i, :, 1])
#     plt.savefig("../plots/n_body/L_evolution/particle_{}.png".format(i))
#     plt.close()

# for i in range(N_timesteps):
#     plt.hist(theta_R_values[:, i], bins=10)
#     plt.ylim(0, 1500)
#     plt.savefig("../plots/n_body/theta_R_distribution/t={}.png".format(i),
#                 dpi=500)
#     plt.close()


# coefficient_norms = np.loadtxt("../data/n_body/test_bfe_coefficient_norms.csv")
# plt.plot(coefficient_norms)
# plt.savefig("../plots/n_body_trial_norms.png")
# # exit()

data = np.loadtxt("../data/n_body/test_trajectories.csv")
data = data.reshape(N_timesteps, N_particles, 2, 2)

for i in range(N_timesteps):
    fig = plt.figure()
    ax = fig.add_subplot(projection='polar')
    c = ax.scatter(data[i, :, 0, 1], data[i, :, 0, 0], s=0.1)
    ax.set_ylim([0, 10])
    plt.savefig("../plots/n_body_trial/t={:04d}Myr".format(i*10))
    plt.close()

